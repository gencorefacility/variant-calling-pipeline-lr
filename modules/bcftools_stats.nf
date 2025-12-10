process BCFTOOLS_STATS {
    tag "$meta.id"
    
    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
    tuple val(meta), path(vcf)
    val(caller)

    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    path "bcftools_stats_*.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${caller}"
    """
    bcftools stats ${vcf} > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > bcftools_stats_${caller}_versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}