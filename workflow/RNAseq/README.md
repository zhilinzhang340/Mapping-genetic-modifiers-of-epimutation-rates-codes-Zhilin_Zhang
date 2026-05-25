# RNA-seq Workflow


This directory documents the RNA-seq preprocessing pipeline for generating gene-level count tables from raw paired-end FASTQ files.

This README is intended as a concise workflow-level description. Some scripts contain absolute paths from the original analysis environment and should be updated before rerunning on another system.

---



## Pipeline Overview

| Step | Description | Tool |
|:-----|:------------|:-----|
| Step 0 | FastQC on raw FASTQ files | FastQC |
| Step 1 | Adapter trimming | Trimmomatic (v0.38) |
| Step 2 | Read alignment to TAIR10 | STAR |
| Step 3 | BAM sorting and indexing | Samtools |
| Step 4 | Read counting per gene | HTSeq-count |


## Files

| File | Description |
|:-----|:------------|
| `run_rnaseq_pipeline.sh` | Bash pipeline: merge FASTQs, QC, trim, align, count (Steps 0-4) |


---

## Input data

The preprocessing script expects paired-end FASTQ files in a user-defined input directory:

```text
raw_data/total_fastq/
```

Expected FASTQ naming pattern:

```text
*_1.fq.gz
*_2.fq.gz
```

The script automatically identifies sample prefixes, merges multiple sequencing chunks belonging to the same sample when needed, or creates symbolic links when only one chunk is present.


---

## Output structure

The preprocessing script writes results under an analysis output directory such as:

```text
analysis_results/
├── 00_merged_fastq/
├── 00_fastqc_raw/
├── 01_trimmed/
├── 02_aligned/
└── 03_counts/
```

| Directory | Description |
|---|---|
| `00_merged_fastq/` | Merged or linked paired-end FASTQ files |
| `00_fastqc_raw/` | FastQC reports for raw FASTQ files |
| `01_trimmed/` | Trimmomatic paired and unpaired trimmed FASTQ files |
| `02_aligned/` | STAR-aligned BAM files, sorted BAM files, and BAM indexes |
| `03_counts/` | HTSeq-count gene-level count files |

---


## Preprocessing workflow

Run the RNA-seq preprocessing pipeline after updating all input paths, software paths, reference genome paths, and thread settings in `run_rnaseq_pipeline.sh`.

Example:

```bash
bash run_rnaseq_pipeline.sh
```

Main settings used in the script:

| Parameter | Value or description |
|---|---|
| Read type | Paired-end RNA-seq |
| Adapter trimming | Trimmomatic |
| Alignment | STAR |
| Reference genome | Arabidopsis TAIR10 STAR index |
| Gene annotation | Arabidopsis TAIR10 GTF |
| BAM processing | samtools sort and index |
| Read counting | HTSeq-count |
| HTSeq strandedness | `--stranded=reverse` |
| HTSeq feature type | `--type=exon` |
| HTSeq mode | `--mode=union` |



---

## Requirements

- STAR (with TAIR10 genome index)
- Trimmomatic (v0.38)
- Samtools
- HTSeq
- FastQC
- Arabidopsis TAIR10 GTF annotation file


## Notes

- Paired-end 150 bp reads, sequenced on Illumina NovaSeq 6000 (Novogene).
- Samples with multiple FASTQ slices are automatically merged by the pipeline script.
- HTSeq-count uses `--stranded=reverse --mode=union --type=exon`.
- File paths in scripts need adjustment for local use.
- Raw sequencing files are not stored directly in this GitHub repository.
- Sample-level information is provided in the top-level metadata file:
  `RIL_WGBS_ATAC-seq-RNA-seq_sample_metadata.csv`.



---
## Citation

Please cite the bioRxiv preprint if you use this repository:

Zhang Z., Wanney W., Xu Y., Zicola J., Hancock A. M., Schmitz R. J., and Johannes F.  
**Mapping genetic modifiers of epimutation rates reveals a punctuated-equilibrium model of CG methylome evolution.**  
bioRxiv. doi: `10.1101/2025.06.14.659605`

---

## Author

Zhilin Zhang, Technical University of Munich
Email: zhilin.zhang@tum.de; zhangzhilin94@gmail.com
