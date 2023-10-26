perl /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/scripts/pbaa_report.pl \
    --fastq /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-clusterReads/execution/hifi_reads-4_fastq_seq_stats.tab \
    --passSeq /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-clusterReads/execution/hifi_reads-4_pbaa_passed_cluster_sequences_stats.tab \
    --failSeq /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-clusterReads/execution/hifi_reads-4_pbaa_failed_cluster_sequences_stats.tab \
    --clusterQC /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-clusterMetrics/execution/hifi_reads-4_clusterQC_report.tsv \
    --clusterVariants /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-variantCall/execution/hifi_reads-4_variant_summary.tsv \
    --clusterOnTargetVariants /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-variantCall/execution/hifi_reads-4_variant_on_target_summary.tsv \
    --prefix hifi_reads-4

perl /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/scripts/deepvariant_report.pl \
    --fastq /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-HifiReadsAlign/execution/hifi_reads-4_raw_hifi_reads_fastq_stats.tab \
    --idxstat /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-HifiReadsAlign/execution/hifi_reads-4_raw_hifi_to_reference_alignment_idxstat.txt \
    --allVariants /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-HifiReadsVarCallDV/execution/hifi_reads-4_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv \
    --onTargetVariants /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/cromwell-executions/main/6cdf7d5a-ca36-4344-af9a-be42c2811edf/call-HifiOnTargetVarsDV/execution/hifi_reads-4_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv \
    --prefix hifi_reads-4