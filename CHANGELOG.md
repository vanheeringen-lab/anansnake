# Changelog

Here, the changes to `anansnake` will be summarized.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

## Unreleased

## [0.1.0] - 2023-11-23

### Changed

- all output folders can now be configured
- running gimme maelstrom is now optional

### Fixed

- ANANSE errors ("No regions overlap")

## [0.0.1] - 2023-01-12

### Added

- run differential gene expression analysis using `DESeq2` from `Seq2Science` RNA-seq output
- run peak enrichment analysis using `gimmemotifs maelstrom` from `Seq2Science` ATAC-seq output
- obtain nonhuman transcription factor binding motifs using `gimmemotifs motif2factors`
- run network differentiation analysis using `ANANSE` from `Seq2Science` ATAC-seq & RNA-seq output
