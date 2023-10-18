version 1.0

# cluster hifi reads using pbAA
task clusterReads {
    
    input {
        File guide_seq
        File hifi_reads_fastq_gz
        String file_label
        Int max_amplicon_size = 20000
        Float min_cluster_frequency = 0.125
        String docker
    }  

    String log_level = "DEBUG"

    command <<<
        set -euo pipefail

        hifireads_file_base=$(basename ~{hifi_reads_fastq_gz} .fastq.gz)
        
        ln -s ~{guide_seq} pbaa_guide.fasta
              
        samtools faidx pbaa_guide.fasta -o pbaa_guide.fasta.fai

        gunzip -c ~{hifi_reads_fastq_gz} > ${hifireads_file_base}.fastq

        mv ${hifireads_file_base}.fastq ~{file_label}.hifi_reads.fastq
        
        samtools faidx --fastq ~{file_label}.hifi_reads.fastq
        
        pbaa cluster --min-cluster-frequency ~{min_cluster_frequency} \
            --max-amplicon-size ~{max_amplicon_size} \
            pbaa_guide.fasta \
            ~{file_label}.hifi_reads.fastq \
            ~{file_label}_pbaa \
            --log-level ~{log_level} \
            --log-file ~{file_label}_pbaa_cluster.log
        
        seqkit stats -a -T ~{file_label}.hifi_reads.fastq > ~{file_label}_fastq_seq_stats.tab
        
        seqkit stats -a -T ~{file_label}_pbaa_passed_cluster_sequences.fasta > ~{file_label}_pbaa_passed_cluster_sequences_stats.tab
        
        seqkit stats -a -T ~{file_label}_pbaa_failed_cluster_sequences.fasta > ~{file_label}_pbaa_failed_cluster_sequences_stats.tab
    >>>

    output {
        File pbaa_passed_cluster_sequences = file_label + "_pbaa_passed_cluster_sequences.fasta"
        File pbaa_failed_cluster_sequences = file_label + "_pbaa_failed_cluster_sequences.fasta"
        File pbaa_read_info = file_label + "_pbaa_read_info.txt"
        File pbaa_run_log = file_label + "_pbaa_cluster.log"
        File fastq_seq_stats = file_label + "_fastq_seq_stats.tab"
        File pbaa_passed_cluster_sequences_stats = file_label + "_pbaa_passed_cluster_sequences_stats.tab"
        File pbaa_failed_cluster_sequences_stats = file_label + "_pbaa_failed_cluster_sequences_stats.tab"
        File hifi_reads_fastq = file_label + ".hifi_reads.fastq"
        File hifi_reads_fastq_index = file_label + ".hifi_reads.fastq.fai"
    }

    runtime {
        docker: "~{docker}"
    }
}