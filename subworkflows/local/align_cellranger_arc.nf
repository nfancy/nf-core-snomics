/*
 * Alignment with Cellranger
 */

include {CELLRANGER_MKGTF} from "../../modules/nf-core/cellranger_arc/arc_mkgtf/main.nf"
include {CELLRANGER_MKREF} from "../../modules/nf-core/cellranger_arc/arc_mkref/main.nf"
include {CELLRANGER_COUNT} from "../../modules/nf-core/cellranger_arc/arc_count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ARC_ALIGN {
    take:
        fasta
        gtf
        cellranger_index
        ch_fastq
        samplesheet

    main:
        ch_versions = Channel.empty()

        assert cellranger_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--fasta') and a gtf file ('--gtf')."

        if (!cellranger_index) {
            // Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_MKGTF( gtf )
            ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

            // Make reference genome
            CELLRANGER_MKREF( fasta, CELLRANGER_MKGTF.out.gtf, "cellranger_reference" )
            ch_versions = ch_versions.mix(CELLRANGER_MKREF.out.versions)
            cellranger_index = CELLRANGER_MKREF.out.reference
        }

        // Obtain read counts
        CELLRANGER_COUNT (
            samplesheet
            cellranger_index
        )
        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_out  = CELLRANGER_COUNT.out.outs
}
