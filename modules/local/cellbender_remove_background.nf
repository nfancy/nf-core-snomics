process CELLBENDER_REMOVE_BACKGROUND {
    tag "$meta.id"
    label (params.gpu ? 'with_gpus': 'process_high')

    container "us.gcr.io/broad-dsde-methods/cellbender:0.3.0"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLBENDER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(input_h5)

    output:
    tuple val(meta), path("${meta.id}*")  ,         emit: cellbender_out
    tuple val(meta), path("${meta.id}*_filtered.h5")  ,    emit: cellbender_h5
    path "versions.yml"                     ,		 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """

    cellbender \\
        remove-background \\
            --input $input_h5 \\
            --output '${meta.id}' \\
            --expected-cells 7000 \\
            --total-droplets-included 30000 \\
            --cuda \\
            --fpr 0.01 \\
            --epochs 150   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(echo \$( cellbender --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}"
    touch ${meta.id}/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(echo \$( cellbender --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
