process RTG_VCFSTATS {

    //conda "bioconda::rtg-tools=3.12.1"
    container 'quay.io/biocontainers/rtg-tools:3.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(vcf)
    val(caller)
    
    output:
    tuple val(meta), path("*.txt"), emit: rtg_stats
    path "rtg_vcfstats_*.yml", emit: versions
    
    script:
    """
    rtg vcfstats ${vcf} > ${meta.id}.${caller}.rtg_stats.txt

    cat <<-END_VERSIONS > rtg_vcfstats_${caller}_versions.yml
    "${task.process}":
        rtg-tools: \$(rtg --version | head -n1 | sed 's/^.*rtg-tools //; s/ .*\$//')
    END_VERSIONS
    """
}
