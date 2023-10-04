# wf-pb-amp : A WDL-based workflow for Variant calling using PacBio HiFi CCS data.

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/anand-imcm/wf-pb-amp/publish.yml)
![GitHub release (with filter)](https://img.shields.io/github/v/release/anand-imcm/wf-pb-amp)

Import the workflow to your Terra workspace using the link below.

- [Dockstore](https://dockstore.org/workflows/github.com/anand-imcm/wf-pb-amp:main?tab=info)

Locate the 'Launch with' widget at the top right of the Dockstore workflow page, and select the 'Terra' platform option. 


## About

A WDL-based workflow that takes PacBio HiFi reads as input and runs PacBio Amplicon Analysis under the hood.

The steps and commands for running the amplicon analysis implemented in this workflow are outlined in the Official PacBio Github [Wiki page](https://github.com/PacificBiosciences/CoSA/wiki/Variant-calling-using-PacBio-HiFi-CCS-data#4c-variant-calling-using-pbaa)


## WDL tasks

- `cluster_to_vcf` [__Task 1__]  Variant calling using [pbaa](https://github.com/PacificBiosciences/pbAA) tool. `pbaa` accepts an input fastq and a guide reference, and outputs a cluster sequence. It runs additional scripts to convert the cluster sequence to VCF format which is used by `VCFCons` tool later in the workflow.
- `reads_to_bam`  [__Task 2__] Aligning HiFi reads to the reference and generate the alignment and coverage summary.
- `call_variants` [__Task 3__] The final task to generate the VCF file after removing low qual variants.


## Inputs

- `reads_fastq_gz` : Input PacBio HiFi reads in .fastq.gz format. ["File (required)"]
- `guide_fasta` : Amplicon reference .fasta file. This reference file is used for pbAA clustering. ["File (required)"]
- `genome_ref` : Reference genome .fasta file. ["File (required)"]
- `genome_index_pbmm` : Reference genome index file generated through pbmm2 in .mmi format. ["File (required)"]
- `prefix` : Sample name. This will be used as prefix for all the output files ["String (required)"]
- `consensus_variant_calling.min_alt_freq` : ["Float (optional, default = 0.5)"]
- `consensus_variant_calling.sort_thread` : ["Int (optional, default = 4)"]
- `consensus_variant_calling.min_coverage` : ["Int (optional, default = 4)"]
- `consensus_variant_calling.alignment_thread` : ["Int (optional, default = 4)"]
- `alignment_metrics.sort_thread` : ["Int (optional, default = 4)"]
- `alignment_metrics.alignment_thread` : ["Int (optional, default = 4)"]
- `amplicon_analysis.min_cluster_read_count` : ["Int (optional, default = 2)"]


## Output

- `pbaa_passed_cluster_sequences`
- `pbaa_failed_cluster_sequences`
- `pbaa_read_info`
- `pbaa_run_log`
- `pbaa_alleles`
- `pbaa_variants`
- `pbaa_vcf`
- `fastq_seq_stats`
- `pbaa_passed_cluster_sequences_stats`
- `pbaa_failed_cluster_sequences_stats`
- `fastqc_report`
- `aligned_sorted_bam`
- `aligned_amplicon_cluster_log`
- `bam_mpileup_report`
- `bam_depth_report`
- `flagstat_report`
- `idxstat_report`
- `vcfcons_fasta`
- `vcfcons_frag_fasta`
- `vcfcons_info`
- `vcfcons_vcf`
- `vcfcons_csv`
- `vcfcons_log`
- `consensus_reference_aligned_bam`
- `consensus_reference_alignment_log`


## Components
- Scripts from `cosa-py3` package extracted from `smrttools-release_12.0.0.177059` bundle
  - consensusVariants.py
  - parser.py
  - pbaa2vcf.py
  - VCFCons.py
- Python packages
  - bio<=1.5.9
  - mappy<=2.26
  - pandas<=2.1.1
  - pysam<=0.21.0
  - PyVCF<=0.6.8
  - scipy<=1.11.3
- Tools
  - samtools<=1.17
  - pbaa<=1.0.3
  - pbmm2<=1.13.0
  - seqkit<=2.5.1
  - bcftools<=1.17
  - setuptools<58
  - gcc<=13.2.0
  - fastqc<=0.12.1
