# wf-pb-amp : A WDL-based workflow for Variant calling using PacBio HiFi reads.

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/anand-imcm/wf-pb-amp/publish.yml)
![GitHub release (with filter)](https://img.shields.io/github/v/release/anand-imcm/wf-pb-amp)

Import the workflow to your Terra workspace using the link below.

- [Dockstore](https://dockstore.org/workflows/github.com/anand-imcm/wf-pb-amp:main?tab=info)

Locate the 'Launch with' widget at the top right of the Dockstore workflow page, and select the 'Terra' platform option. 


## About

A WDL-based workflow that takes PacBio HiFi reads as input and runs PacBio Amplicon Analysis under the hood.

This pipeline is based on the amplicon analysis outlined in the Official PacBio GitHub [hifi-amplicon-workflow](https://github.com/PacificBiosciences/hifi-amplicon-workflow), which is originally a Snakemake-based workflow. Our implementation is based on WDL (Workflow Description Language) and includes some customizations to meet the specific requirements of our data.


## Workflow steps

- Cluster hifi reads with [`pbaa`](https://github.com/PacificBiosciences/pbAA) tool.
- Align cluster consensus to the reference.
- Call variants per cluster.
- Annotate variants per cluster and generate a summary.
- Extract the hifi reads based on pbaa clusters.
- Align the clustered hifi reads to the reference.
- Color-code BAM records, after aligning the clustered hifi reads.
- Generate the cluster qc report using the lima summary records.
- Align the raw hifi reads to reference.
- Call variants after aligning the raw reads using [`deepvariant`](https://github.com/google/deepvariant) tool.
- Annotate the variants reported by `deepvariant`
- Generate summary table using total and on-target variants.


## Inputs

- `reads_fastq_gz` : Input PacBio HiFi reads in .fastq.gz format. "File"
- `guide_fasta` : Amplicon reference .fasta file. This reference file is used for pbAA clustering. "File"
- `genome_ref` : Human reference genome .fasta format. "File"
- `genome_index_pbmm` : Reference index generated through pbmm2 in .mmi format. "File"
- `clinvar_vcf` : Clinvar vcf for annotation in .gz format. "File"
- `features_gff` : Reference geneome annotation in .gff3.gz format. "File"
- `target_bed` : Coordinates for the amplified regions (target) in .bed format. "File"
- `lima_report` : Lima report file obtained from the demultiplexing process in .lima.report format. "File"
- `prefix` : Sample name. This will be used as prefix for all the output files. "String"
- `max_amplicon_size` : `pbaa cluster` filter option. "Int (optional, default = 20000)"
- `min_cluster_frequency` : `pbaa cluster` filter option. "Float (optional, default = 0.125)"


## Output

- Input Fastq stats
  - `fastq_seq_stats`
- `pbaa` cluster results
  - `pbaa_failed_cluster_sequences`
  - `pbaa_failed_cluster_sequences_stats`
  - `pbaa_passed_cluster_sequences`
  - `pbaa_passed_cluster_sequences_stats`
  - `pbaa_read_info`
  - `pbaa_run_log`
- Variant call and annotation from the consensus sequences per cluster
  - `consensus_to_reference_alignment_bam`
  - `consensus_to_reference_alignment_flagstat`
  - `consensus_to_reference_alignment_idxstat`
  - `consensus_to_reference_alignment_log`
  - `raw_vcf`
  - `annotated_vcf`
- Variant summary from all clusters
  - `variant_summary`
  - `variant_on_target_summary`
- Alignment results using the hifi reads clustered by `pbaa`
  - `clustered_hifi_fastq`
  - `clustered_hifi_reads_fastq_stats`
  - `clustered_hifi_to_reference_alignment_painted_bam`
  - `clustered_hifi_to_reference_alignment_painted_bam_index`
  - `clustered_hifi_to_reference_alignment_painted_idxstat`
  - `clustered_hifi_to_reference_alignment_painted_flagstat`
  - `clustered_hifi_to_reference_alignment_bampaint_log`
  - `clustered_hifi_to_reference_alignment_log`
  - `clusterQC_report`
- Alignment results of hifi reads to a reference
  - `raw_hifi_to_reference_alignment_bam`
  - `raw_hifi_to_reference_alignment_index`
  - `raw_hifi_to_reference_alignment_log`
  - `raw_hifi_to_reference_alignment_flagstat`
  - `raw_hifi_to_reference_alignment_idxstat`
  - `raw_hifi_reads_fastq_stats`
- Variant call and annotation from the hifi reads
  - `raw_hifi_to_reference_alignment_all_variants_vcf`
  - `raw_hifi_to_reference_alignment_all_variants_annotated_vcf`
  - `raw_hifi_to_reference_alignment_all_variants_annotated_summary`
  - `raw_hifi_to_reference_alignment_pass_variants_annotated_vcf`
  - `raw_hifi_to_reference_alignment_pass_variants_annotated_summary`
- Variant summary from the hifi reads
  - `raw_hifi_to_reference_alignment_all_variants_stats`
  - `raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_vcf`
  - `raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary`
- Summary
  - `pbaa_sequence_summary`
  - `pbaa_unique_variant_summary`
  - `germline_sequence_summary`
  - `germline_unique_variant_summary`


## Components

- Python packages
  - pandas<=2.1.1
  - pysam<=0.21.0
- Tools
  - bcftools<=1.18
  - bedtools<=2.31.0
  - bzip2<=1.0.8
  - fastqc<=0.12.1
  - pbaa<=1.0.3
  - pbmm2<=1.13.0
  - samtools<=1.18
  - seqkit<=2.5.1
  - seqtk<=1.4
- Containers
  - ghcr.io/anand-imcm/wf-pb-amp
  - google/deepvariant