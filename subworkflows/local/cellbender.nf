/*
 * Alignment with Cellranger
 */

include {CELLBENDER_REMOVE_BACKGROUND} from "../../modules/local/cellbender_remove_background.nf"
include {CONVERT_H5} from "../../modules/local/convert_h5_to_sparse_mat.nf"

// Define workflow to remove background using cellbender and convert h5 matrix to feature_bc_matrix
workflow CELLBENDER {
    take:
        input_h5
        ensembl_mapping

    main:
        ch_versions = Channel.empty()

        // Obtain read counts
        CELLBENDER_REMOVE_BACKGROUND (
            input_h5
        )
        ch_versions = ch_versions.mix(CELLBENDER_REMOVE_BACKGROUND.out.versions)
        ch_cellbender_h5 = CELLBENDER_REMOVE_BACKGROUND.out.cellbender_h5


        CONVERT_H5 ( 
            ch_cellbender_h5,
            ensembl_mapping
        )

        ch_versions = ch_versions.mix(CONVERT_H5.out.versions)


    emit:
        ch_versions
}
