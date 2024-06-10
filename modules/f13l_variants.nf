process F13L_VARIANTS {
    tag "${meta.id}}"
    label 'process_medium'
    container 'sereewit/ivar:1.4.2_49079e2'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path ref_index
    path gff
    val save_mpileup

    output:
    tuple val(meta), path("*.tsv"),     emit: tsv
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def f13l_region = 'ON563414.3:39094-40212'
    """
    samtools \\
        mpileup \\
        $args2 \\
        --reference $ref \\
        --region ${f13l_region} \\
        $bam \\
        $mpileup \\
        | ivar \\
            variants \\
            $args \\
            -g $gff \\
            -r $ref \\
            -p $prefix

    edit_f13l_variants.py \\
        ${prefix}.tsv \\
        $gff \\
        ${prefix}.reformatted
    """
}
