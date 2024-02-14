process CELLRANGER_N_CELL {
    tag       "$meta.id"
    label     'process_low'
    container "nfancy/snomics:4.2.3"

    input:
    tuple val(meta), path(metrics_summary)

    output:
    tuple val(meta), stdout,                emit: n_cells
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/snomics/bin/
    """
    get_cellranger_out_n_cell.r \\
        --metrics_summary $metrics_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         R: \$( R --version | grep "R version" | cut -d' ' -f3 )
    END_VERSIONS
    """
}
