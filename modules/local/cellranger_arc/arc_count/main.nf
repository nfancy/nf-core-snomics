process CELLRANGER_ARC_COUNT {
    tag       "$id"
    label     'process_high'
    container "austins2/cellranger-arc:v2.0.0" 

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(id), path(acc_reads, stageAs: 'fastqs/accessibility/*'), path(gex_reads, stageAs: 'fastqs/expression/*')
    path  reference

    output:
    tuple val(id), path("${id}/outs/*"), emit: outs
    path "versions.yml"                , emit: versions

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
    echo "\${PWD}/fastqs/expression,${id},Gene Expression" >> cellranger_arc_samplesheet.csv
    echo "\${PWD}/fastqs/accessibility,${id},Chromatin Accessibility" >> cellranger_arc_samplesheet.csv

    ### run cellranger-arc

    cellranger-arc \\
        count \\
        --id=$id \\
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
    mkdir -p "${id}/outs/"
    touch sample-${id}/outs/fake_file.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

}
