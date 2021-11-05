nextflow.enable.dsl=2

include {check_barcodes; grab_fastq; expunge; racon} from "./modules/preprocess_utils.nf"
include {flye_assemble; kraken2_reads; kraken2_contigs; bracken} from "./modules/processing_utils.nf"

workflow preprocess_expunge{
    take:
        runpath
        barcodes

    main:       
        grab_fastq(runpath, barcodes)
        
        expunge(grab_fastq.out.barcode, grab_fastq.out.all_reads)

        
    emit:
        reads = expunge.out.expunged_reads
        barcode = expunge.out.barcode
}

workflow preprocess {
    take:
        runpath
        barcodes

    main:
        barcodes
            .count()
            .set{count}
        

        grab_fastq(runpath, barcodes)

    emit:
        reads = grab_fastq.out.all_reads
        barcode = grab_fastq.out.barcode
}

workflow assemble {
    take:
        barcode
        reads

    main:
        flye_assemble(barcode, reads)
    
    emit:
        contigs = flye_assemble.out.contigs
        barcode = flye_assemble.out.barcode
        reads = flye_assemble.out.reads
    
}

workflow classify_contigs {
    take:
        barcode
        contigs
    
    main:
        kraken2_contigs(barcode, contigs)
    
    emit:
        contig_outfiles = kraken2_contigs.out
}

workflow polish {
    take:
        barcode
        reads
        contigs
    
    main:
        racon(barcode, reads, contigs)
    
    emit:
        polished_contigs = racon.out.polished_contigs
        barcode = barcode
}

workflow classify_and_estimate_abundance {
    take:
        barcode
        reads
    
    main:
        kraken2_reads(barcode, reads)
        bracken(barcode, kraken2_reads.out.kraken_report)
}