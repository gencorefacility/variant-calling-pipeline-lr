process CLAIR3 {
    tag "$meta.id"
    label 'process_high'
    
    //conda "bioconda::clair3=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.4--py39hf5e1c6e_0' :
        'hkubal/clair3:v1.0.4' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai
    path model

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "clair3_versions.yml"           , emit: versions

    when:
    params.run_clair3 == true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    run_clair3.sh \\
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=$task.cpus \\
        --platform=$params.clair3_platform \\
        --model_path=$model \\
        --output=clair3_output \\
        $args

    mv clair3_output/merge_output.vcf.gz ${prefix}.clair3.vcf.gz
    mv clair3_output/merge_output.vcf.gz.tbi ${prefix}.clair3.vcf.gz.tbi

    cat <<-END_VERSIONS > clair3_versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | grep "Clair3" | sed 's/Clair3 //')
    END_VERSIONS
    """
}