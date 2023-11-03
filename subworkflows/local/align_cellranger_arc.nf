/*
 * Alignment with Cellranger-arc
 */

include { CELLRANGER_ARC_MKGTF } from "../../modules/local/cellranger_arc/arc_mkgtf/main.nf"
include { CELLRANGER_ARC_MKREF } from "../../modules/local/cellranger_arc/arc_mkref/main.nf"
include { CELLRANGER_ARC_COUNT } from "../../modules/local/cellranger_arc/arc_count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ARC_ALIGN {

    take:
        fasta
        gtf
        cellranger_arc_index
        ch_fastq
    
    main:
        ch_versions = Channel.empty()

        
        assert cellranger_arc_index || (fasta && gtf):
            "Must provide either a cellranger index or both a fasta file ('--fasta') and a gtf file ('--gtf')."

        if ( !cellranger_arc_index ) {
            //Filter GTF based on gene biotypes passed in params.modules
            CELLRANGER_ARC_MKGTF( gtf )
            ch_versions = ch_versions.mix( CELLRANGER_ARC_MKGTF.out.versions )

            // Make reference genome
            CELLRANGER_ARC_MKREF( fasta, CELLRANGER_ARC_MKGTF.out.gtf, "cellranger_arc_reference" )
            ch_versions          = ch_versions.mix( CELLRANGER_ARC_MKREF.out.versions )
            cellranger_arc_index = CELLRANGER_ARC_MKREF.out.reference
        }
        
        // Group accessibility and expression fastqs by ID
        // Currently blocking by groupTuple

        ch_fastq
            .map { [ it[0].id, [ it[0].library_type, it[1].flatten() ] ] } // Create list with library type key
            .groupTuple() // Group ID
            .map { [ it[0], it[1].sort{ a, b ->  a[0] <=> b[0] } ] } // Sort by library type key
            .map { [ it[0], it[1][0][1], it[1][1][1] ] } // Expand lists to fastqs
            .set { ch_fastq }

        // Obtain read counts
        CELLRANGER_ARC_COUNT (
            ch_fastq,
            cellranger_arc_index
        )
        
        ch_versions = ch_versions.mix( CELLRANGER_ARC_COUNT.out.versions )

    emit:
       ch_versions
       cellranger_out = CELLRANGER_ARC_COUNT.out.outs

}