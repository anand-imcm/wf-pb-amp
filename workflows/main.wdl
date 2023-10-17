version 1.0

import "./tasks/cluster_sequences.wdl" as generate_consensus
import "./tasks/align_concensus.wdl" as align_consensus_reads
import "./tasks/call_variants.wdl" as variant_stats
import "./tasks/extract_reads.wdl" as extract
import "./tasks/align_hifi_cluster.wdl" as align_clustered_reads

workflow main {

    # String pipeline_version = "1.0.0"
    # String container_src = "ghcr.io/anand-imcm/wf-pb-amp:~{pipeline_version}"

    input {
        File reads_fastq_gz
        File guide_fasta
        File genome_ref
        File genome_index_pbmm
        File clinvar_vcf
        File features_gff
        File target_bed
        String prefix
    }

    parameter_meta {
        reads_fastq_gz : "Input PacBio HiFi reads in .fastq.gz format."
        guide_fasta : "Amplicon reference .fasta file. This reference file is used for pbAA clustering."
        genome_ref : "Human reference genome .fasta file."
        genome_index_pbmm : "Reference index generated through pbmm2 in .mmi format."
        prefix : "Sample name. This will be used as prefix for all the output files."
    }

    call generate_consensus.clusterReads {
        input: guide_seq = guide_fasta, hifi_reads_fastq_gz = reads_fastq_gz, file_label = prefix
    }

    call align_consensus_reads.alignConsensus {
        input: genome_reference = genome_ref, pbmm2_index = genome_index_pbmm, passed_cluster_seq = clusterReads.pbaa_passed_cluster_sequences, file_label = prefix
    }

    call variant_stats.variantCall {
        input: consensus_to_ref_aligned_bam = alignConsensus.consensus_to_reference_alignment_bam, genome_reference = genome_ref, clinvar = clinvar_vcf, gff = features_gff, bed = target_bed, file_label = prefix
    }

    call extract.extractClusteredHifiReads {
        input: hifi_fastq = clusterReads.hifi_reads_fastq, hifi_fastq_index = clusterReads.hifi_reads_fastq_index, pbaa_read_info = clusterReads.pbaa_read_info, file_label = prefix
    }

    call align_clustered_reads.alignClusteredHifiReads {
        input: genome_reference= genome_ref, pbmm2_index = genome_index_pbmm, clustered_hifi_reads = extractClusteredHifiReads.clustered_hifi_fastq, pbaa_read_info = clusterReads.pbaa_read_info, file_label = prefix
    }

    # output {
    #     File pbaa_passed_cluster_sequences = amplicon_analysis.pbaa_passed_cluster_sequences
    #     File pbaa_failed_cluster_sequences = amplicon_analysis.pbaa_failed_cluster_sequences
    #     File pbaa_read_info = amplicon_analysis.pbaa_read_info
    #     File pbaa_run_log = amplicon_analysis.pbaa_run_log
    #     File pbaa_alleles = amplicon_analysis.pbaa_alleles
    #     File pbaa_variants = amplicon_analysis.pbaa_variants
    #     File pbaa_vcf = amplicon_analysis.pbaa_vcf
    #     File fastq_seq_stats = amplicon_analysis.fastq_seq_stats
    #     File pbaa_passed_cluster_sequences_stats = amplicon_analysis.pbaa_passed_cluster_sequences_stats
    #     File pbaa_failed_cluster_sequences_stats = amplicon_analysis.pbaa_failed_cluster_sequences_stats

    #     File fastqc_report = alignment_metrics.fastqc_report
    #     File aligned_sorted_bam = alignment_metrics.aligned_sorted_bam
    #     File aligned_amplicon_cluster_log = alignment_metrics.aligned_amplicon_cluster_log
    #     File bam_mpileup_report = alignment_metrics.bam_mpileup_report
    #     File bam_depth_report = alignment_metrics.bam_depth_report
    #     File flagstat_report = alignment_metrics.flagstat_report
    #     File idxstat_report = alignment_metrics.idxstat_report

    #     File vcfcons_fasta = consensus_variant_calling.vcfcons_fasta
    #     File vcfcons_frag_fasta = consensus_variant_calling.vcfcons_frag_fasta
    #     File vcfcons_info = consensus_variant_calling.vcfcons_info
    #     File vcfcons_vcf = consensus_variant_calling.vcfcons_vcf
    #     File vcfcons_csv = consensus_variant_calling.vcfcons_csv
    #     File vcfcons_log = consensus_variant_calling.vcfcons_log
    #     File consensus_reference_aligned_bam = consensus_variant_calling.consensus_reference_aligned_bam
    #     File consensus_reference_alignment_log = consensus_variant_calling.consensus_reference_alignment_log
    # }

    meta {
        description: "A WDL-based workflow for Variant calling using PacBio HiFi CCS data."
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }

}