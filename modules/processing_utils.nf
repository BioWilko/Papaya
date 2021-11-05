nextflow.enable.dsl=2

process flye_assemble {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/", mode: 'copy', overwrite: true
    input:
        val(barcode)
        file(reads)
    output:
        val barcode, emit: barcode
        path "${barcode}_all_contigs.fasta", emit: contigs
        path "${barcode}_assembly_report.txt"
        path "${reads}", emit: reads

    """
    flye --nano-hq ${reads} --meta --out-dir . -t ${params.threads}
    mv assembly.fasta ${barcode}_all_contigs.fasta
    mv assembly_info.txt ${barcode}_assembly_report.txt
    """
}

process kraken2_contigs {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/", mode: 'copy', overwrite: true
    input:
        val(barcode)
        file(contigs)
    output:
        file("*")
        
 
    """
    kraken2 --db ${params.kraken_database} --threads ${params.threads} --classified-out ${barcode}_classified_contigs.fasta --unclassified-out ${barcode}_unclassified_contigs.fasta --report ${barcode}_contig_kraken.report ${contigs} > kraken_contig.report
    """
    
}

process kraken2_reads {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/", mode: 'copy', overwrite: true
    input:
        val(barcode)
        file(reads)

    output:
        path "*.fasta"
        path "${barcode}_read_kraken.report", emit: kraken_report
  
 
    """
    kraken2 --db ${params.kraken_database} --threads ${params.threads} --classified-out ${barcode}_classified_reads.fasta --unclassified-out ${barcode}_unclassified_reads.fasta --report ${barcode}_read_kraken.report ${reads}
    """
    
}

process bracken {
    conda "environments/base.yaml"
    publishDir "${params.outdir}/${barcode}/" , mode: 'copy', overwrite: true
    input:
        val(barcode)
        file(kraken_report)
    output:
        path "${barcode}_bracken.report"

    """
    bracken -d ${params.kraken_database} -i ${kraken_report} -o ${barcode}_bracken.report -r 300 -l ${params.bracken_tax_level}
    """
    
}