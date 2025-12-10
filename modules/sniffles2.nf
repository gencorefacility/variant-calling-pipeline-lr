// Use task.ext.args to pass additional arguments to Sniffles2
// Example in config: process { withName: 'SNIFFLES2' { ext.args = '--minsvlen 30 --minsupport 3' } }
// Note: tandem_repeats input is OPTIONAL - pass an empty list [] if not available
process SNIFFLES2 {
    tag "$meta.id"
    label 'process_medium'
    
    //conda "bioconda::sniffles=2.2"
    container 'quay.io/biocontainers/sniffles:2.7.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_fai
    path tandem_repeats  // OPTIONAL: pass [] or empty channel if not available

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.snf")   , emit: snf
    path "sniffles2_versions.yml"    , emit: versions

    when:
    params.run_sniffles2 == true

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Only add tandem-repeats flag if file is provided (not empty list or NO_FILE placeholder)
    def tandem_repeats_arg = (tandem_repeats && tandem_repeats.name != 'NO_FILE') ? "--tandem-repeats $tandem_repeats" : ""
    
    """
    sniffles \\
        --input $bam \\
        --vcf ${prefix}.sniffles2.vcf.gz \\
        --snf ${prefix}.snf \\
        --threads $task.cpus \\
        --reference $reference \\
        $tandem_repeats_arg \\
        $args

    cat <<-END_VERSIONS > sniffles2_versions.yml
    "${task.process}":
        sniffles2: \$(sniffles --version 2>&1 | sed 's/Sniffles //')
    END_VERSIONS
    """
}