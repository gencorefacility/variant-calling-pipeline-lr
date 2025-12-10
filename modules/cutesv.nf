process CUTESV {
    tag "$meta.id"
    label 'process_medium'
    
    //conda "bioconda::cutesv=2.0.3"
    container 'quay.io/biocontainers/cutesv:2.1.3--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "cutesv_versions.yml"    , emit: versions

    when:
    params.run_cutesv == true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    mkdir -p cutesv_tmp

    cuteSV \\
        $bam \\
        $reference \\
        ${prefix}.cutesv.vcf \\
        cutesv_tmp \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > cutesv_versions.yml
    "${task.process}":
        cutesv: \$(cuteSV --version 2>&1 | sed 's/cuteSV //')
    END_VERSIONS
    """
}