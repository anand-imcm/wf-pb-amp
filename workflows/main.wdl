version 1.0

import "./tasks/cluster_sequences.wdl" as generate_consensus
import "./tasks/align_concensus.wdl" as align_consensus_reads
import "./tasks/call_variants.wdl" as variant_stats
import "./tasks/extract_reads.wdl" as extract
import "./tasks/align_hifi_cluster.wdl" as align_clustered_reads
import "./tasks/cluster_barcode_qc.wdl" as cluster_qc

workflow main {

    String pipeline_version = "1.2.0"
    String container_src = "ghcr.io/anand-imcm/wf-pb-amp:latest:~{pipeline_version}"

    input {
        File reads_fastq_gz
        File guide_fasta
        File genome_ref
        File genome_index_pbmm
        File clinvar_vcf
        File features_gff
        File target_bed
        File lima_report
        String prefix
    }

    parameter_meta {
        reads_fastq_gz : "Input PacBio HiFi reads in .fastq.gz format."
        guide_fasta : "Amplicon reference .fasta file. This reference file is used for pbAA clustering."
        genome_ref : "Human reference genome .fasta file."
        genome_index_pbmm : "Reference index generated through pbmm2 in .mmi format."
        clinvar_vcf : "Clinvar vcf for annotation in .gz format."
        features_gff : "Reference geneome annotation in .gff3.gz format."
        target_bed : "Coordinates for the amplified regions (target) in .bed format."
        lima_report : "Lima report file obtained from the demultiplexing process in .lima.report format."     
        prefix : "Sample name. This will be used as prefix for all the output files."
    }

    call generate_consensus.clusterReads {
        input: guide_seq = guide_fasta, hifi_reads_fastq_gz = reads_fastq_gz, file_label = prefix, docker = container_src
    }

    call align_consensus_reads.alignConsensus {
        input: genome_reference = genome_ref, pbmm2_index = genome_index_pbmm, passed_cluster_seq = clusterReads.pbaa_passed_cluster_sequences, file_label = prefix, docker = container_src
    }

    call variant_stats.variantCall {
        input: consensus_to_ref_aligned_bam = alignConsensus.consensus_to_reference_alignment_bam, genome_reference = genome_ref, clinvar = clinvar_vcf, gff = features_gff, bed = target_bed, file_label = prefix, docker = container_src
    }

    call extract.extractClusteredHifiReads {
        input: hifi_fastq = clusterReads.hifi_reads_fastq, hifi_fastq_index = clusterReads.hifi_reads_fastq_index, pbaa_read_info = clusterReads.pbaa_read_info, file_label = prefix, docker = container_src
    }

    call align_clustered_reads.alignClusteredHifiReads {
        input: genome_reference = genome_ref, pbmm2_index = genome_index_pbmm, clustered_hifi_reads = extractClusteredHifiReads.clustered_hifi_fastq, pbaa_read_info = clusterReads.pbaa_read_info, file_label = prefix, docker = container_src
    }

    call cluster_qc.clusterMetrics {
        input: clustered_holes = extractClusteredHifiReads.clustered_holes, lima_report = lima_report, pbaa_read_info = clusterReads.pbaa_read_info, file_label = prefix, docker = container_src
    }

    output {
        File fastq_seq_stats = clusterReads.fastq_seq_stats
        File pbaa_failed_cluster_sequences = clusterReads.pbaa_failed_cluster_sequences
        File pbaa_failed_cluster_sequences_stats = clusterReads.pbaa_failed_cluster_sequences_stats
        File pbaa_passed_cluster_sequences = clusterReads.pbaa_passed_cluster_sequences
        File pbaa_passed_cluster_sequences_stats = clusterReads.pbaa_passed_cluster_sequences_stats
        File pbaa_read_info = clusterReads.pbaa_read_info
        File pbaa_run_log = clusterReads.pbaa_run_log

        File consensus_to_reference_alignment_bam = alignConsensus.consensus_to_reference_alignment_bam
        File consensus_to_reference_alignment_flagstat = alignConsensus.consensus_to_reference_alignment_flagstat
        File consensus_to_reference_alignment_idxstat = alignConsensus.consensus_to_reference_alignment_idxstat
        File consensus_to_reference_alignment_log = alignConsensus.consensus_to_reference_alignment_log

        Array[File] raw_vcf = variantCall.raw_vcf
        Array[File] annotated_vcf = variantCall.annotated_vcf
        File variant_summary = variantCall.variant_summary
        File variant_on_target_summary = variantCall.variant_on_target_summary

        File clustered_hifi_fastq = extractClusteredHifiReads.clustered_hifi_fastq

        File clustered_hifi_reads_fastq_stats = alignClusteredHifiReads.clustered_hifi_reads_fastq_stats
        File clustered_hifi_to_reference_alignment_painted_bam = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_painted_bam
        File clustered_hifi_to_reference_alignment_painted_bam_index = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_painted_bam_index
        File clustered_hifi_to_reference_alignment_painted_idxstat = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_painted_idxstat
        File clustered_hifi_to_reference_alignment_painted_flagstat = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_painted_flagstat
        File clustered_hifi_to_reference_alignment_bampaint_log = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_bampaint_log
        File clustered_hifi_to_reference_alignment_log = alignClusteredHifiReads.clustered_hifi_to_reference_alignment_log

        File clusterQC_report = clusterMetrics.clusterQC_report
    }

    meta {
        description: "A WDL-based workflow for Variant calling and annotation using PacBio HiFi reads."
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }

}