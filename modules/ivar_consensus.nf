process IVAR_CONSENSUS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/ivar:1.4--h6b7c446_1'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*consensus*.fa"),    emit: consensus
    tuple val(meta), path("*.qual.txt"),        emit: qual
    tuple val(meta), path("*.mpileup"),         emit: mpileup, optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        mpileup \\
        --reference $ref \\
        $args2 \\
        $bam \\
        | ivar \\
            consensus \\
            $args \\
            -p $prefix \\
            -i $prefix

    # removing leading and trailing Ns
    awk '/^>/ {print; next} {sub(/^N+/, ""); printf "%s", \$0} END {print ""}' ${prefix}.fa > ${prefix}_temp.fa
    sed '/^>/!s/N\\+\$//' ${prefix}_temp.fa > ${prefix}.fa
    rm ${prefix}_temp.fa
    """
}
