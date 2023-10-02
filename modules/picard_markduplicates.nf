process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    container 'biocontainers/picard:3.0.0--hdfd78af_1'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"),   emit: bam
    tuple val(meta), path("*.metrics.txt"),          emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.bam \\
        --REFERENCE_SEQUENCE $ref \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt
    """
}
