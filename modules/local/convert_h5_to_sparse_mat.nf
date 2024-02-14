process CONVERT_H5 {
    tag "$meta.id"
    label 'process_medium'

    container "nfancy/snomics:4.2.3"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "this pipeline does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(cellbender_h5)
    path ensembl_mapping

    output:
    tuple val(meta), path("cellbender_feature_bc_matrix/*")  ,        emit: cellbender_feature_bc_matrix
    path "versions.yml"                     ,		 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/snomics/bin/
    """
    convert_h5_to_sparse_mat.r \\
        --input_file $cellbender_h5 \\
        --ensembl_mapping $ensembl_mapping

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | grep "R version" | cut -d' ' -f3 )
    END_VERSIONS
    """
}