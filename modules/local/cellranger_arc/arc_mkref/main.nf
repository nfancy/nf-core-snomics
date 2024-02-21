process CELLRANGER_ARC_MKREF {
    tag       "$fasta"
    label     'process_high'
    container "cumulusprod/cellranger-arc:2.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path fasta
    path gtf
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo -e "{\n\
    \tgenome: ['$reference_name']\n\
    \tinput_fasta: ['$fasta']\n\
    \tinput_gtf: ['$gtf']\n\
    \tnon_nuclear_contigs: ['$params.non_nuclear_contigs']\n\
    \tinput_motifs: '$params.input_motifs'\n}" > arc_mkref.config

    grep -v "\"null\"" arc_mkref.config > tmpfile && mv tmpfile arc_mkref.config
    
    # need check for required arguments

    cellranger-arc \\
        mkref \\
        --config=arc_mkref.config \\
        --memgb=${task.memory.toGiga()}\\
        --nthreads=$task.cpus\\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
    
}
