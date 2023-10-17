process CELLRANGER_ARC_COUNT {
    tag       "$meta_gex.id"
    label     'process_high'
    container "austins2/cellranger-arc:v2.0.0" 

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta_gex), path("cellranger_arc_fastqs/gex/${meta_gex.id}_??_.fastq.gz")
    tuple val(meta_acc), path("cellranger_arc_fastqs/acc/${meta_acc.id}_??_.fastq.gz")
    path  reference

    output:
    tuple val(meta_gex), path("${meta_gex.id}/outs/*"), emit: outs
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    echo "fastqs, sample, library_type" > cellranger_arc_samplesheet.csv
    echo "\${PWD}/cellranger_arc_fastqs/gex,${meta_gex.id},Gene Expression" >> cellranger_arc_samplesheet.csv
    echo "\${PWD}/cellranger_arc_fastqs/acc,${meta_acc.id},Chromatin Accessibility" >> cellranger_arc_samplesheet.csv

    mapfile -d ' ' -t fastqs < <(find . -name '*.fastq.gz')
    
    for fastq in \$fastqs; do
        end=\$(echo \$(readlink \$fastq) | sed -r 's@^.*(_S[0-9]_L[0-9]{3}_R[1-2]_001.fastq.gz\$)@\\1@')
        renamed=\$(echo \$fastq | sed -r "s@_[0-9]{2}_.fastq.gz\\\$@\$end@")
        mv \$fastq \$renamed
    done

    cellranger-arc \\
        count \\
        --id=$meta_gex.id \\
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

}
