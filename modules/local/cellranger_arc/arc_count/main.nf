process CELLRANGER_ARC_COUNT {
    tag       "$meta.id"
    label     'process_high'
    container "austins2/cellranger-arc:v2.0.0" 

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(acc_reads, stageAs: 'fastqs/accessibility/*'), path(gex_reads, stageAs: 'fastqs/expression/*')
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs/*"), emit: outs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def all_reads = (acc_reads + gex_reads).join(' ')
    """
    # rename to cellranger file names
    for fastq in ${all_reads}; do
        renamed=\$(echo \$fastq | sed -r "s@_([0-9]).merged.fastq.gz@_S1_L001_R\\1_001.fastq.gz@")
        mv \$fastq \$renamed
    done

    ### output cellranger-arc csv

    echo "fastqs, sample, library_type" > cellranger_arc_samplesheet.csv
    echo "\${PWD}/fastqs/expression,${meta.id},Gene Expression" >> cellranger_arc_samplesheet.csv
    echo "\${PWD}/fastqs/accessibility,${meta.id},Chromatin Accessibility" >> cellranger_arc_samplesheet.csv

    ### run cellranger-arc

    cellranger-arc \\
        count \\
        --id=$meta.id \\
        --reference=$reference \\
        --libraries=cellranger_arc_samplesheet.csv \\
        --localcores=$task.cpus \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/outs/
    touch ${meta.id}/outs/fake-gex_possorted_bam.bam.bai \
          ${meta.id}/outs/fake-atac_fragments.tsv.gz \
          ${meta.id}/outs/fake-per_barcode_metrics.csv \
          ${meta.id}/outs/fake-atac_fragments.tsv.gz.tbi \
          ${meta.id}/outs/fake-atac_cut_sites.bigwig" \
          ${meta.id}/outs/fake-web_summary.html" \
          ${meta.id}/outs/fake-atac_possorted_bam.bam" \
          ${meta.id}/outs/fake-atac_possorted_bam.bam.bai" \
          ${meta.id}/outs/fake-atac_peak_annotation.tsv" \
          ${meta.id}/outs/fake-filtered_feature_bc_matrix.h5 \
          ${meta.id}/outs/fake-cloupe.cloupe \
          ${meta.id}/outs/fake-raw_feature_bc_matrix.h5 \
          ${meta.id}/outs/fake-summary.csv \
          ${meta.id}/outs/fake-gex_molecule_info.h5 \
          ${meta.id}/outs/fake-gex_possorted_bam.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

}
