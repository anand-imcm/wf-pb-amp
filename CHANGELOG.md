# Changelog

All notable changes to this project will be documented in this file.

## [1.3.5] - 2023-11-01

### Fixed
- Alignment metrics from pbmm2 log file

## [1.3.4] - 2023-10-26

### Added
- New task to generate summary:
  - sequence level summary for amplicon and germline methods
  - variant summary for amplicon and germline variant calls

## [1.3.3] - 2023-10-24

### Fixed
- pbaa bampaint issue -> `ERROR -|- BamPainterRunner -|- 0x7f8bce5ab740|| -|- No reads tagged` . It results in an exit code of 1, leading to the workflow being marked as failed.
- pbaa bampaint output is made optional since the input bam files may have no reads.

## [1.3.2] - 2023-10-24

### Fixed
- Bcftools annotation issue by adding a flag to process the GTs as is, create haplotypes regardless of their phase. 

## [1.3.1] - 2023-10-23

### Added
- New modules for variant calling using [`deepvariant`](https://github.com/google/deepvariant).
- Added a summary of variant calls using total and ontarget variant.

## [1.2.0] - 2023-10-18

### Added

- Breaking changes in this release
- New workflow based on PacBio's [hifi-amplicon-workflow](https://github.com/PacificBiosciences/hifi-amplicon-workflow)
- New modules
  - clusterReads
  - alignConsensus
  - variantCall
  - extractClusteredHifiReads
  - alignClusteredHifiReads
  - clusterMetrics
- The required inputs for the workflow has also been changed. Refer the documentation for usage instructions.


## [1.0.0] - 2023-10-04

### Added

- Initial release