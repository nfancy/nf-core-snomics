process SAMPLESHEET_CHECK {
    tag       "${samplesheet}|${aligner}"
    label     'process_single'
    container "nfancy/snomics:4.2.3"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "this pipeline does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path samplesheet
    val aligner

    output:
    path "checked_samplesheet.csv", emit: checked_samplesheet
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/snomics/bin/
    """
    check_samplesheet.r \\
        --input $samplesheet \\
        --aligner $aligner

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         R: \$( R --version | grep "R version" | cut -d' ' -f3 )
    END_VERSIONS
    """
}
