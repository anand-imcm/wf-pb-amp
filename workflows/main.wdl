version 1.0
import "./tasks/cluster_variant_call.wdl" as cluster_to_vcf
import "./tasks/alignment_and_metrics.wdl" as reads_to_bam
import "./tasks/consensus_vcf_and_metrics.wdl" as call_variants

workflow main {

    String container_src = "pbaaenv2:latest"

    input {
        File reads_fastq_gz
        File guide_fasta
        File genome_ref
        File genome_index_pbmm
        String prefix
    }

    call cluster_to_vcf.amplicon_analysis {
        input : pbaa_guide_fasta = guide_fasta, genome_reference = genome_ref, amplicons_fastq_gz = reads_fastq_gz, file_label = prefix, docker = container_src
    }
    call reads_to_bam.alignment_metrics {
        input: genome_reference = genome_ref, genome_index = genome_index_pbmm, amplicons_fastq_gz = reads_fastq_gz, file_label = prefix, docker = container_src
    }
    call call_variants.consensus_variant_calling {
        input : genome_reference = genome_ref, genome_index = genome_index_pbmm, bam_depth = alignment_metrics.bam_depth_report, pbaa_vcf = amplicon_analysis.pbaa_vcf, file_label = prefix, docker = container_src
    }

    output {
        File pbaa_passed_cluster_sequences = amplicon_analysis.pbaa_passed_cluster_sequences
        File pbaa_failed_cluster_sequences = amplicon_analysis.pbaa_failed_cluster_sequences
        File pbaa_read_info = amplicon_analysis.pbaa_read_info
        File pbaa_run_log = amplicon_analysis.pbaa_run_log
        File pbaa_alleles = amplicon_analysis.pbaa_alleles
        File pbaa_variants = amplicon_analysis.pbaa_variants
        File pbaa_vcf = amplicon_analysis.pbaa_vcf
        File fastq_seq_stats = amplicon_analysis.fastq_seq_stats
        File pbaa_passed_cluster_sequences_stats = amplicon_analysis.pbaa_passed_cluster_sequences_stats
        File pbaa_failed_cluster_sequences_stats = amplicon_analysis.pbaa_failed_cluster_sequences_stats

        File fastqc_report = alignment_metrics.fastqc_report
        File aligned_sorted_bam = alignment_metrics.aligned_sorted_bam
        File aligned_amplicon_cluster_log = alignment_metrics.aligned_amplicon_cluster_log
        File bam_mpileup_report = alignment_metrics.bam_mpileup_report
        File bam_depth_report = alignment_metrics.bam_depth_report
        File flagstat_report = alignment_metrics.flagstat_report
        File idxstat_report = alignment_metrics.idxstat_report

        File vcfcons_fasta = consensus_variant_calling.vcfcons_fasta
        File vcfcons_frag_fasta = consensus_variant_calling.vcfcons_frag_fasta
        File vcfcons_info = consensus_variant_calling.vcfcons_info
        File vcfcons_vcf = consensus_variant_calling.vcfcons_vcf
        File vcfcons_csv = consensus_variant_calling.vcfcons_csv
        File vcfcons_log = consensus_variant_calling.vcfcons_log
        File consensus_reference_aligned_bam = consensus_variant_calling.consensus_reference_aligned_bam
        File consensus_reference_alignment_log = consensus_variant_calling.consensus_reference_alignment_log
    }
}






# to do
# bam - alignment metrics: samtools flagstat, idxstat
# vcf - bcftools stats

# get the summary metrics
# task summary {

# }