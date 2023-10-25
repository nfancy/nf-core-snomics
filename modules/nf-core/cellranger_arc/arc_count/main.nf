process CELLRANGER_ARC_COUNT {
    tag       "$meta_gex.id"
    label     'process_high'
    container "austins2/cellranger-arc:v2.0.0" 

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGER_ARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta_gex), path("gex/${meta_gex.id}__??__.fastq.gz")
    tuple val(meta_acc), path("acc/${meta_acc.id}__??__.fastq.gz")
    path  reference

    output:
    tuple val(meta_gex), path("${meta_gex.id}/outs/*"), emit: outs
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    fastqs=(gex/*.fastq.gz acc/*.fastq.gz)

    ### obtain full fastq paths

    for fastq in \${fastqs[@]}; do
        link=\$(readlink \$fastq) 
        echo \$link >> tmp_links
        links+=(\$link)
    done

    ### iterate over directories, finding the first directory not in common (indicating more than one flowcell)

    dir_depth=\$(head -n 1 tmp_links | tr '/' '\n' | wc -l)

    for col_num in \$(seq 1 \$dir_depth); do
        unique=(\$(cut -d '/' -f \$col_num tmp_links | sort | uniq))
        unique_count=\${#unique[@]}
        if [ \$unique_count -gt 1 ]; then
            break
        fi
    done
    rm tmp_links

    ### Match the fastq files belonging to each flowcell

    fc=1
    fastq_tracker=()
    for dir in \${unique[@]}; do
        match_idx=()
        match_fastqs=()
        for i in \${!links[@]}; do
            link=\${links[\$i]}
            if [[ "\$link" == *"\$dir"* ]]; then
                match_idx+=("\$i")
                bn=\$(basename \$link)
                match_fastqs+=(\$bn)
            fi
        done

        ### Check if names are repeated (flowcell check)

        common=(\$(echo \${fastq_tracker[@]} \${match_fastqs[@]} | tr ' ' '\n' | sort | uniq -d))
        if [ \${#common[@]} -gt 0 ]; then
            fc=\$(( \$fc+1 ))
        fi
        fastq_tracker+=(\${match_fastqs[@]})

        ### Move symlink fastqs to new location with corect name

        mkdir -p ./flowcell_\${fc}/acc ./flowcell_\${fc}/gex
        for idx in \${match_idx[@]}; do
            ending=\$(echo \${links[\$idx]} | sed -r 's@^.*(_S[0-9]_L[0-9]{3}_R[1-2]_001.fastq.gz)@\\1@')
            renamed=\$(echo \${fastqs[\$idx]} | sed -r "s@__[0-9]{2}__.fastq.gz@\$ending@")
            mv \${fastqs[\$idx]} flowcell_\${fc}/\${renamed}
        done
        rmdir ./flowcell_\${fc}/acc ./flowcell_\${fc}/gex --ignore-fail-on-non-empty
    done

    rmdir acc gex

    ### output cellranger-arc csv

    echo "fastqs, sample, library_type" > cellranger_arc_samplesheet.csv
    dirs=(flowcell*/*)
    for dir in \${dirs[@]}; do
        if [[ \$dir == *"gex" ]]; then
            lib_type='Gene Expression'
        else
            lib_type='Chromatin Accessibility'
        fi
        echo "\${PWD}/\${dir},${meta_gex.id},\${lib_type}" >> cellranger_arc_samplesheet.csv
    done

    ### run cellranger-arc

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

    stub:
    """
    mkdir -p "${meta_gex.id}/outs/"
    touch sample-${meta_gex.id}/outs/fake_file.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

}
