# ATAC-seq Analysis Pipeline for Arabidopsis thaliana

This repository contains a modular pipeline for processing ATAC-seq data from Arabidopsis thaliana recombinant inbred lines (RILs). The pipeline includes quality control, adapter trimming, read alignment, deduplication, BAM-to-BED conversion, quality summary, and peak calling using MACS2.

## Pipeline Overview

| Step | Description | Script |
|------|-------------|--------|
| Step 0 | Run FastQC on raw FASTQ files | `step0_run_fastqc.sh` |
| Step 1 | Trim adapters using Trimmomatic | `step1_trim_adapters.sh` |
| Step 2 | Align reads to reference genome using BWA and sort/index with samtools | `step2_bismark_alignment.sh` |
| Step 3.1 | Remove duplicates with Picard | `step3-1_remove_duplicates.sh` |
| Step 3.2 | Convert BAM to BED and perform post-dedup FastQC | `step3-2_convert_and_fastqc.sh` |
| Step 4 | Call peaks using MACS2 | `step4_peak_calling.sh` |

## Requirements

- BWA
- Trimmomatic v0.39
- Picard
- Samtools
- Bedtools
- FastQC
- MACS2
- Arabidopsis thaliana genome (FASTA file)

## File Structure

RIL_meth/
├── ATAC/
│ ├── raw_fastq/ # Input FASTQ files
│ ├── pipeline_output/
│ │ ├── fastqc/ # FastQC results (Step 0)
│ │ ├── step1_trimmed/ # Trimmed FASTQs (Step 1)
│ │ ├── step2_alignment/ # BAM alignments (Step 2)
│ │ ├── step3_dedup/ # Deduplicated BAMs (Step 3.1)
│ │ │ └── bed_fastqc/ # BED and QC (Step 3.2)
│ │ └── step4_macs2/ # MACS2 peaks (Step 4)



## Notes

- Input FASTQ files should be named as `*_1.fastq.gz` and `*_2.fastq.gz` for paired-end reads.
- Reference genome should be in FASTA format and indexed for BWA.
- You may need to adjust file paths and software locations within each script before execution.

## Author

Zhilin Zhang  
Technical University of Munich  
Contact: zhilin.zhang@tum.de


