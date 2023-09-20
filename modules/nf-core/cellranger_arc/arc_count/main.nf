process CELLRANGER_COUNT {
    tag "$meta.gem"
    label 'process_high'

    container "cumulusprod/cellranger-arc:2.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(reads)
    path  reference
    path  samplesheet

    output:
    tuple val(meta), path("sample-${meta.gem}/outs/*"), emit: outs
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_name = reference.name

    """
    mkdir -p cellranger_arc_fastqs
    mkdir -p cellranger_arc_fastqs/gex
    mkdir -p cellranger_arc_fastqs/acc

    # add check fastq naming structure

    mapfile -t info < <(awk -F ',' -v OFS='\n' 'NR!=1 {print $1, $2, $3, $4}' $samplesheet)

    for ((i=0;i<${#info[@]}/4;i++))
    do 
        renamed_fq1=$(echo "${info[$i*4+1]}" | sed "s/.*_S1/${info[$i*4]}_S1/")
        renamed_fq2=$(echo "${info[$i*4+2]}" | sed "s/.*_S1/${info[$i*4]}_S1/")

        if [[ "${info[(($i*4+3))]}" = "Gene Expression"* ]]
            then
                ln -sf ${info[$i*4+1]} cellranger_arc_fastqs/gex/$renamed_fq1
                ln -sf ${info[$i*4+2]} cellranger_arc_fastqs/gex/$renamed_fq2
            else
                ln -sf ${info[$i*4+1]} cellranger_arc_fastqs/acc/$renamed_fq1
                ln -sf ${info[$i*4+2]} cellranger_arc_fastqs/acc/$renamed_fq2
        fi
    done

    echo "fastqs, sample, library_type" > cellranger_arc_samplesheet.csv
    echo "${PWD}/cellranger_arc_fastqs/gex,,Gene Expression" >> cellranger_arc_samplesheet.csv
    echo "${PWD}/cellranger_arc_fastqs/acc,,Chromatin Accessibility" >> cellranger_arc_samplesheet.csv

    # run
    cellranger-arc \\
        count \\
        --id='sample-${meta.gem}' \\
        --reference=$reference_name \\
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
    mkdir -p "sample-${meta.gem}/outs/"
    touch sample-${meta.gem}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
