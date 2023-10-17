version 1.0

# align clustered hifi reads to reference and generate the final bam using pbaa bampaint
task alignClusteredHifiReads {
    
    input {
        File genome_reference
        File pbmm2_index
        File clustered_hifi_reads
        File pbaa_read_info
        String file_label
    }  

    String log_level = "DEBUG"

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{pbmm2_index} genome_reference.mmi
        
        pbmm2 align \
            --preset hifi \
            --sort \
            ~{genome_reference} \
            ~{clustered_hifi_reads} \
            ~{file_label}_clustered_hifi_to_reference_alignment.bam \
            --log-level ~{log_level} \
            --log-file ~{file_label}_clustered_hifi_to_reference_alignment.log
        
        pbaa bampaint \
            ~{pbaa_read_info} \
            ~{file_label}_clustered_hifi_to_reference_alignment.bam \
            ~{file_label}_clustered_hifi_to_reference_alignment_painted.bam \
            --log-level ~{log_level} \
            --log-file ~{file_label}_clustered_hifi_to_reference_alignment_pbaa_bampaint.log
        
        seqkit stats -a -T ~{clustered_hifi_reads} > ~{file_label}_clustered_hifi_reads_fastq_stats.tab
        
        samtools flagstat \
            ~{file_label}_clustered_hifi_to_reference_alignment_painted.bam > ~{file_label}_clustered_hifi_to_reference_alignment_painted_flagstat.txt
        
        samtools idxstat \
            ~{file_label}_clustered_hifi_to_reference_alignment_painted.bam > ~{file_label}_clustered_hifi_to_reference_alignment_idx.txt

        printf "ref_contig\tref_contig_len\ttotal_mapped_reads\ttotal_unmapped_reads\n" | cat - ~{file_label}_clustered_hifi_to_reference_alignment_idx.txt > ~{file_label}_clustered_hifi_to_reference_alignment_painted_idxstat.txt
    >>>

    output {
        File clustered_hifi_to_reference_alignment_painted_bam = file_label + "_clustered_hifi_to_reference_alignment_painted.bam"
        File clustered_hifi_to_reference_alignment_painted_bam_index = file_label + "_clustered_hifi_to_reference_alignment_painted.bam.bai"
        File clustered_hifi_to_reference_alignment_log = file_label + "_clustered_hifi_to_reference_alignment.log"
        File clustered_hifi_to_reference_alignment_bampaint_log = file_label + "_clustered_hifi_to_reference_alignment_pbaa_bampaint.log"
        File clustered_hifi_to_reference_alignment_painted_flagstat = file_label + "_clustered_hifi_to_reference_alignment_painted_flagstat.txt"
        File clustered_hifi_to_reference_alignment_painted_idxstat = file_label + "_clustered_hifi_to_reference_alignment_painted_idxstat.txt"
        File clustered_hifi_reads_fastq_stats = file_label + "_clustered_hifi_reads_fastq_stats.tab"
    }
}