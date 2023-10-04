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
}






# to do
# bam - alignment metrics: samtools flagstat, idxstat
# vcf - bcftools stats

# get the summary metrics
# task summary {

# }