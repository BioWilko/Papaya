nextflow.enable.dsl=2



def printHelp() {
    log.info"""
    
    """
}


if(!params.runpath) error "No filepath to run directory supplied with --runpath"

if(!params.outdir) error "No output directory supplied with --outdir"

if(params.help) {
    printHelp()
    exit 0
}

include {check_barcodes; grab_fastq; expunge; racon} from "./modules/preprocess_utils.nf"

include {preprocess; preprocess_expunge; assemble; classify_contigs; classify_and_estimate_abundance; polish} from "./workflows.nf"

workflow {
    barcodes = Channel.of(params.barcodes)
    runpath = params.runpath
    // genome = params.genome
    check_barcodes(runpath, barcodes)

    check_barcodes.out
        .splitText()
        .map{it -> it.trim()}
        .set{barcodes}
    
    main:
        if(params.expunge){
            preprocess_expunge(runpath, barcodes)
            preprocess_expunge.out.barcode
                .set{barcode}
            preprocess_expunge.out.reads
                .set{reads}
        } else {
            preprocess(runpath, barcodes)
            preprocess.out.barcode
                .set{barcode}
            preprocess.out.reads
                .set{reads}
        }
        if(params.assemble){
            assemble(barcode, reads)
            if(params.polish){
                polish(assemble.out.barcode, assemble.out.reads, assemble.out.contigs)
                polish.out.polished_contigs
                    .set{processed_contigs}
                polish.out.barcode
                    .set{processed_barcode}
            } else {
                assemble.out.contigs
                    .set{processed_contigs}
                assemble.out.barcode
                    .set{processed_barcode}
            }
            classify_contigs(processed_barcode, processed_contigs)
        }
        classify_and_estimate_abundance(barcode, reads)
}