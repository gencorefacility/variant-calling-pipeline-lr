process MINIMAP2 {
    tag "$meta.id"
    label 'process_high'
    
    //conda "bioconda::minimap2=2.26 bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    tuple val(meta), path("*.stats"), emit: stats
    path "minimap2_versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    minimap2 \\
        -ax $params.minimap2_preset \\
        -t $task.cpus \\
        $args \\
        $reference \\
        $reads | \\
    samtools sort \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -
    
    samtools index -@ $task.cpus ${prefix}.bam

    samtools stats -@ $task.cpus ${prefix}.bam > ${prefix}.stats

    cat <<-END_VERSIONS > minimap2_versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}