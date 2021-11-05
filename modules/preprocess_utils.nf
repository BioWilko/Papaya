nextflow.enable.dsl=2

process check_barcodes {
    conda "environments/base.yaml"
    input:
        path(runpath)
        val(barcodes)
    output:
        file("barcodes.txt")
    
    script:
        if(barcodes == 0){
            """
            python_funcs.py -a -r ${runpath} > barcodes.txt
            """
        } else {
            """
            IFS=',' read -r -a barcode_array <<< "${barcodes}"
            for barcode in \${barcode_array[@]}; do
                echo \$barcode >> barcodes.txt
            done
            """
        }
}

process grab_fastq {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/", mode: 'copy', overwrite: true
    input:
        path(runpath)
        val(barcode)
    output:
        path "${barcode}_all_reads.fastq", emit: all_reads
        val barcode, emit: barcode

    """
    python_funcs.py -b ${barcode} -r ${runpath}
    cat *.fastq > ${barcode}_all_reads.fastq
    """
}

process expunge {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/", mode: 'copy', overwrite: true
    input:
        val(barcode)
        file(fastq)
    output:
        path "${barcode}_expunged_reads.fastq", emit: expunged_reads
        val barcode, emit: barcode
    
    """
    minimap2 -x map-ont ${params.expunge_ref_idx} ${fastq} > aligned_reads.paf
    python_funcs.py -x aligned_reads.paf > to_expunge
    grep -E "^@" ${fastq} | cut -c2- | sed 's/ .*//' > read.ids
    grep -vxFf to_expunge read.ids > to_save
    seqtk subseq ${fastq} to_save > ${barcode}_expunged_reads.fastq
    """
}

process racon {
    conda "environments/base.yaml"
    input:
        val(barcode)
        val(reads)
        val(contigs)
    output:
        path "${barcode}_polished_contigs.fasta", emit: polished_contigs
    
    """
    minimap2 -x map-ont ${contigs} ${reads} > overlaps_1.paf
    racon -t ${params.threads} ${reads} overlaps_1.paf ${contigs} > temp_contigs.fasta
    minimap2 -x map-ont temp_contigs.fasta ${reads} > overlaps_2.paf
    racon -t ${params.threads} ${reads} overlaps_2.paf temp_contigs.fasta > ${barcode}_polished_contigs.fasta
    """
}