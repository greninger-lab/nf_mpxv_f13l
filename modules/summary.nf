process SUMMARY {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jefffurlong/revica:ubuntu-20.04'

    input:
    tuple val(meta), path(consensus), path(map_ref_bam), path(map_ref_bai), path(fastp_trim_log)

    output:
    path("*.tsv"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def f13l_length = '1119'

    """
    # raw reads and trimmed reads
    raw_reads=\$(grep -A1 "before filtering:" ${fastp_trim_log} | grep 'total reads:' | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}')
    trimmed_reads=\$(grep -A1 "after filtering:" ${fastp_trim_log} | grep 'total reads:' | cut -d: -f2 | tr -d " " | awk 'NF{sum+=\$1} END {print sum}')
    pct_reads_trimmed=\$(echo "\${trimmed_reads}/\${raw_reads}*100" | bc -l)
    pct_reads_trimmed_formatted=\$(printf "%.2f" "\${pct_reads_trimmed}")

    # mapped reads
    mapped_reads=\$(samtools view -F 4 -c ${map_ref_bam})
    pct_reads_mapped=\$(echo "\${mapped_reads}/\${raw_reads}*100" | bc -l)
    pct_reads_mapped_formatted=\$(printf "%.2f" "\${pct_reads_mapped}")

    # whole genome coverage
    pct_genome_covered=\$(samtools coverage ${map_ref_bam} | awk 'NR>1' | cut -f6)
    pct_genome_covered_formatted=\$(printf "%.2f" "\${pct_genome_covered}")
    mean_genome_coverage=\$(samtools coverage ${map_ref_bam} | awk 'NR>1' | cut -f7)
    mean_genome_coverage_formatted=\$(printf "%.2f" "\${mean_genome_coverage}")
    
    # F13L coverage
    pct_F13L_covered=\$(samtools coverage -r "ON563414.3:39094-40212" ${map_ref_bam} | awk 'NR>1' | cut -f6)
    pct_F13L_covered_formatted=\$(printf "%.2f" "\${pct_F13L_covered}")
    mean_F13L_coverage=\$(samtools coverage -r "ON563414.3:39094-40212" ${map_ref_bam} | awk 'NR>1' | cut -f7)
    mean_F13L_coverage_formatted=\$(printf "%.2f" "\${mean_F13L_coverage}")
    num_bases_F13L_50x=\$(samtools depth -r "ON563414.3:39094-40212" ${map_ref_bam} | awk '{if(\$3>50)print\$3}' | wc -l)
    pct_F13L_50x=\$(echo "\${num_bases_F13L_50x}/${f13l_length}*100" | bc -l)
    pct_F13L_50x_formatted=\$(printf "%.2f" "\${pct_F13L_50x}")
    num_bases_F13L_100x=\$(samtools depth -r "ON563414.3:39094-40212" ${map_ref_bam} | awk '{if(\$3>100)print\$3}' | wc -l)
    pct_F13L_100x=\$(echo "\${num_bases_F13L_100x}/${f13l_length}*100" | bc -l)
    pct_F13L_100x_formatted=\$(printf "%.2f" "\${pct_F13L_100x}")

    # consensus genome
    consensus_length=\$(awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${consensus} | awk 'FNR==2{print val,\$1}')
    num_ns_consensus=\$(grep -v "^>" ${consensus} | tr -c -d N | wc -c)
    pct_ns=\$(echo "\${num_ns_consensus}/\${consensus_length}*100" | bc -l | awk 'FNR==1{print val,\$1}')
    pct_ns_formatted=\$(printf "%.4f" "\${pct_ns}")
    num_as_consensus=\$(grep -v "^>" ${consensus} | tr -c -d A | wc -c)
    num_cs_consensus=\$(grep -v "^>" ${consensus} | tr -c -d C | wc -c)
    num_gs_consensus=\$(grep -v "^>" ${consensus} | tr -c -d G | wc -c)
    num_ts_consensus=\$(grep -v "^>" ${consensus} | tr -c -d T | wc -c)
    num_non_ns_ambiguous=\$(echo "\${consensus_length}-\${num_as_consensus}-\${num_cs_consensus}-\${num_gs_consensus}-\${num_ts_consensus}-\${num_ns_consensus}" | bc -l)
    

    echo "sample_name\traw_reads\ttrimmed_reads\tpct_reads_trimmed\tmapped_reads\tpct_reads_mapped\tpct_genome_covered\tmean_genome_coverage\tpct_F13L_covered\tmean_F13L_coverage\tpct_F13L_50x\tpct_F13L_100x\tconsensus_length\tnum_ns\tpct_ns\tnum_ambiguous" > ${prefix}_summary.tsv
    echo "${prefix}\t\$raw_reads\t\$trimmed_reads\t\${pct_reads_trimmed_formatted}\t\${mapped_reads}\t\${pct_reads_mapped_formatted}\t\${pct_genome_covered_formatted}\t\${mean_genome_coverage_formatted}\t\${pct_F13L_covered_formatted}\t\${mean_F13L_coverage_formatted}\t\${pct_F13L_50x_formatted}\t\${pct_F13L_100x_formatted}\t\${consensus_length}\t\${num_ns_consensus}\t\${pct_ns_formatted}\t\${num_non_ns_ambiguous}" >> ${prefix}_summary.tsv
    """
}
