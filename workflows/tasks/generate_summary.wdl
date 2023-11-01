version 1.0

# extract hifi reads from cluster read info
task report {
    
    input {
        File fastq_seq_stats
        File pbaa_passed_cluster_sequences_stats
        File pbaa_failed_cluster_sequences_stats
        File clusterQC_report
        File variant_summary
        File variant_on_target_summary
        File raw_hifi_reads_fastq_stats
        File raw_hifi_to_reference_alignment_log
        File raw_hifi_to_reference_alignment_pass_variants_annotated_summary
        File raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary
        String file_label
        String docker
    }

    command <<<
        set -euo pipefail

        perl /scripts/pbaa_report.pl \
            --fastq ~{fastq_seq_stats} \
            --passSeq ~{pbaa_passed_cluster_sequences_stats} \
            --failSeq ~{pbaa_failed_cluster_sequences_stats} \
            --clusterQC ~{clusterQC_report} \
            --clusterVariants ~{variant_summary} \
            --clusterOnTargetVariants ~{variant_on_target_summary} \
            --prefix ~{file_label}
        
        perl /scripts/deepvariant_report.pl \
            --fastq ~{raw_hifi_reads_fastq_stats} \
            --pbmmlog ~{raw_hifi_to_reference_alignment_log} \
            --allVariants ~{raw_hifi_to_reference_alignment_pass_variants_annotated_summary} \
            --onTargetVariants ~{raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary} \
            --prefix ~{file_label}
    >>>

    output {
        File pbaa_sequence_summary = file_label + "_amplicon_workflow_sequence_summary.tsv"
        File pbaa_unique_variant_summary = file_label + "_amplicon_workflow_variants_summary.tsv"
        File germline_sequence_summary = file_label + "_germline_workflow_sequence_summary.tsv"
        File germline_unique_variant_summary = file_label + "_germline_workflow_variants_summary.tsv"
    }

    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}