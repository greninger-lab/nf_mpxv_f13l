process BBMAP_ALIGN_CONSENSUS {
    tag "${meta.id}}"
    label 'process_medium_java'
    container 'quay.io/jefffurlong/revica:ubuntu-20.04'

    input:
    tuple val(meta), path(fastq)
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"),    emit: bam
    tuple val(meta), path(consensus),                                   emit: consensus
    tuple val(meta), path(fastq),                                       emit: reads
    tuple val(meta), path("*.log"),                                     emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"

    // Set the ref_seq variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containg a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( consensus.isDirectory() ) {
        if ( consensus ==~ /(.\/)?consensus\/?/ ) {
            ref_seq = ''
        } else {
            ref_seq = "path=${consensus}"
        }
    } else {
        ref_seq = "ref=${consensus}"
    }

    """
    bbmap.sh \\
        ${ref_seq} \\
        $input \\
        out=${prefix}.bam \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbmap.log

    # remove unmapped reads, sort and index bam files  
    samtools view -b -F 4 -@ ${task.cpus} ${prefix}.bam -o ${prefix}_mapped.bam           
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${prefix}_mapped.bam  
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam 
    """
}
