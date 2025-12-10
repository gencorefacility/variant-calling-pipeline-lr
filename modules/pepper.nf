process PEPPER {
    tag "$meta.id"
    label 'process_high'
    
    //conda "bioconda::pepper-deepvariant=0.8.0"
    container 'docker://kishwars/pepper_deepvariant:r0.8'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "pepper_versions.yml"           , emit: versions

    when:
    params.run_pepper == true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    run_pepper_margin_deepvariant call_variant \\
        -b $bam \\
        -f $reference \\
        -o output \\
        -t $task.cpus \\
        $params.pepper_preset \\
        $args

    mv output/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz ${prefix}.pepper.vcf.gz
    mv output/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz.tbi ${prefix}.pepper.vcf.gz.tbi

    cat <<-END_VERSIONS > pepper_versions.yml
    "${task.process}":
        pepper: \$(run_pepper_margin_deepvariant --version 2>&1 | grep "PEPPER" | sed 's/PEPPER: //')
    END_VERSIONS
    """
}