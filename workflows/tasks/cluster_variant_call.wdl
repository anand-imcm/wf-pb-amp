version 1.0

# cluster using pbAA and variant calling
task amplicon_analysis {
    input {
        File pbaa_guide_fasta
        File genome_reference
        File amplicons_fastq_gz
        String file_label
        String docker
        Int min_cluster_read_count = 2
    }  

    String log_level = "DEBUG"

    command <<<

        set -euo pipefail
        
        fastq_gz_inp=~{file_label}.fastq.gz # replace ".fastq.gz" with ".hifi_reads.fastq.gz"
        
        fastq_inp=~{file_label}.fastq # .hifi_reads.fastq

        samtools faidx ~{pbaa_guide_fasta}

        samtools faidx ~{genome_reference}

        gunzip -c ~{amplicons_fastq_gz} > $fastq_inp

        samtools faidx --fastq $fastq_inp

        pbaa cluster \
            --log-level ~{log_level} \
            --log-file ~{file_label}_pbaa.log \
            --trim-ends 5 \
            --min-cluster-read-count ~{min_cluster_read_count} \
            ~{pbaa_guide_fasta} ${fastq_inp} ~{file_label}_pbaa

        # convert pbaa outcome to VCF
        python3 /scripts/consensusVariants.py \
            --runName ~{file_label}_pbaa \
            --prefix ~{file_label}_pbaa \
            --read_info ~{file_label}_pbaa_read_info.txt \
            --hifiSupport ${fastq_inp} \
            ~{genome_reference} \
            ~{file_label}_pbaa_passed_cluster_sequences.fasta > ~{file_label}_consensusVariants.log
        
        python3 /scripts/pbaa2vcf.py \
            --passOnly -s barcode \
            -o ~{file_label}.vcf \
            ~{file_label}_pbaa_alleles.csv \
            ~{file_label}_pbaa_variants.csv \
            ~{genome_reference} > ~{file_label}_pbaa2vcf.log
        

        seqkit stats -a -T ${fastq_inp} > ~{file_label}_fastq_seq_stats.tab
        seqkit stats -a -T ~{file_label}_pbaa_passed_cluster_sequences.fasta > ~{file_label}_pbaa_passed_cluster_sequences_stats.tab
        seqkit stats -a -T ~{file_label}_pbaa_failed_cluster_sequences.fasta > ~{file_label}_pbaa_failed_cluster_sequences_stats.tab
    >>>

    output {
        File pbaa_passed_cluster_sequences = file_label + "_pbaa_passed_cluster_sequences.fasta"
        File pbaa_failed_cluster_sequences = file_label + "_pbaa_failed_cluster_sequences.fasta"
        File pbaa_read_info = file_label + "_pbaa_read_info.txt"
        File pbaa_run_log = file_label + "_pbaa.log"
        File pbaa_alleles = file_label + "_pbaa_alleles.csv"
        File pbaa_variants = file_label + "_pbaa_variants.csv"
        File pbaa_vcf = file_label + ".vcf"
        File fastq_seq_stats = file_label + "_fastq_seq_stats.tab"
        File pbaa_passed_cluster_sequences_stats = file_label + "_pbaa_passed_cluster_sequences_stats.tab"
        File pbaa_failed_cluster_sequences_stats = file_label + "_pbaa_failed_cluster_sequences_stats.tab"
    }

    parameter_meta {

    }

    runtime {
        docker: "~{docker}"
    }
}