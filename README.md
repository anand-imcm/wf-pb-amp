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

- Cluster hifi reads with [pbaa](https://github.com/PacificBiosciences/pbAA) tool.
- Align cluster consensus to the reference.
- Call variants per cluster.
- Annotate variants per cluster and generate a summary.
- Extract the hifi reads based on pbaa clusters.
- Align the clustered hifi reads to the reference.
- Color-code BAM records, after aligning the clustered hifi reads.
- Generate the cluster qc report using the lima summary records.


## Inputs

- `reads_fastq_gz` : "Input PacBio HiFi reads in .fastq.gz format."
- `guide_fasta` : "Amplicon reference .fasta file. This reference file is used for pbAA clustering."
- `genome_ref` : "Human reference genome .fasta file."
- `genome_index_pbmm` : "Reference index generated through pbmm2 in .mmi format."
- `clinvar_vcf` : "Clinvar vcf for annotation in .gz format."
- `features_gff` : "Reference geneome annotation in .gff3.gz format."
- `target_bed` : "Coordinates for the amplified regions (target) in .bed format."
- `lima_report` : "Lima report file obtained from the demultiplexing process in .lima.report format."     
- `prefix` : "Sample name. This will be used as prefix for all the output files."
- `max_amplicon_size`": "Int (optional, default = 20000)",
- `min_cluster_frequency`": "Float (optional, default = 0.125)"


## Output

- `fastq_seq_stats`
- `pbaa_failed_cluster_sequences`
- `pbaa_failed_cluster_sequences_stats`
- `pbaa_passed_cluster_sequences`
- `pbaa_passed_cluster_sequences_stats`
- `pbaa_read_info`
- `pbaa_run_log`
- `consensus_to_reference_alignment_bam`
- `consensus_to_reference_alignment_flagstat`
- `consensus_to_reference_alignment_idxstat`
- `consensus_to_reference_alignment_log`
- `raw_vcf`
- `annotated_vcf`
- `variant_summary`
- `variant_on_target_summary`
- `clustered_hifi_fastq`
- `clustered_hifi_reads_fastq_stats`
- `clustered_hifi_to_reference_alignment_painted_bam`
- `clustered_hifi_to_reference_alignment_painted_bam_index`
- `clustered_hifi_to_reference_alignment_painted_idxstat`
- `clustered_hifi_to_reference_alignment_painted_flagstat`
- `clustered_hifi_to_reference_alignment_bampaint_log`
- `clustered_hifi_to_reference_alignment_log`
- `clusterQC_report`


## Components

- Python packages
  - pandas<=2.1.1
  - pysam<=0.21.0
- Tools
  - bcftools<=1.18
  - bedtools<=2.31.0
  - bzip2<=1.0.8
  - pbaa<=1.0.3
  - pbmm2<=1.13.0
  - samtools<=1.18
  - seqkit<=2.5.1
  - seqtk<=1.4
