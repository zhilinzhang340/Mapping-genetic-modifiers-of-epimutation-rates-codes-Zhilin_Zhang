# Mapping Genetic Modifiers of Epimutation Rates

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20384172.svg)](https://doi.org/10.5281/zenodo.20384172)

Code and data for:

> Zhang Z, Wanney W, Xu Y, Zicola J, Hancock AM, Schmitz RJ\*, Johannes F\*. **Mapping genetic modifiers of epimutation rates reveals a punctuated-equilibrium model of CG methylome evolution.** *bioRxiv* (2025). https://doi.org/10.1101/2025.06.14.659605




## Overview

This repository is organized for the **June 2025 preprint version** of the manuscript. It documents how the main preprint figures were generated and provides workflow-level documentation for WGBS preprocessing, epimutation-rate estimation, QTL mapping, ATAC-seq analysis, and RNA-seq analysis.


## Purpose of this repository

This repository is intended to make the analyses underlying the preprint transparent and reproducible. It provides:

1. A figure-by-figure map linking each script to the corresponding preprint figure panel.
2. The processed data tables used by each script.
3. The expected output figure files.
4. Workflow documentation for WGBS preprocessing, methylation-state calling, AlphaBeta rate estimation, QTL mapping, ATAC-seq analysis, and RNA-seq analysis.

## Data and Code Availability

| Resource | Description | Link |
|:---------|:------------|:-----|
| Project code | Analysis scripts and workflow notebooks | This repository |
| WGBS data | 371 MAML samples | GEO [GSE296957](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296957) |
| ATAC-seq data | Selected RILs (RIL50, RIL149) | GEO [GSE297077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297077) |
| MethylStar | WGBS preprocessing and Methimpute calling | [GitHub](https://github.com/jlab-code/MethylStar) |
| AlphaBeta | Pedigree-based epimutation rate estimation | [GitHub](https://github.com/jlab-code/AlphaBeta) |

## Sample Metadata

`RIL_WGBS_ATAC-seq-RNA-seq_sample_metadata.csv` provides a complete sample index for all sequencing data (WGBS, ATAC-seq, RNA-seq), including sample names, RIL IDs, genotype backgrounds, generation numbers, lineage information, processed file names, and raw FASTQ file names.

---

## Repository Structure

```
.
├── README.md
├── RIL_WGBS_ATAC-seq-RNA-seq_sample_metadata.csv
├── scripts/                  # R scripts used to generate the main preprint figures 
├── data/                     # Processed input tables and selected intermediate outputs used by the figure scripts
├── figures/                  # PDF outputs generated from the figure scripts 
└── workflow/                 # Preprocessing pipelines
    ├── WGBS/                 # WGBS preprocessing, methylation-state calling, AlphaBeta rate estimation, and QTL mapping workflow documentation 
    ├── ATAC/                 # ATAC-seq preprocessing and peak-calling workflow documentation 
    └── RNAseq/               # RNA-seq preprocessing workflow documentation 
```

---

---

## Script-to-Figure Mapping

All scripts are in `scripts/` and read input from `data/`. Panels without scripts (1b, 1f, 4a, 4e) were prepared manually.

### Figure 1: Genotype-dependent methylation variation and MAML population

| Figure | Script | Input data | Description |
|:-------|:-------|:-----------|:------------|
| 1a | `plot_Figure1A.R` | `850_samples_CG_methy_level_used_2026-01.txt`, `RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt`, `new_observed_methy_level_using_proportion_MA-accession_2023-11.txt` | mCG distribution across 850 accessions |
| 1c | `plot_Figure1C.R` | `Col-0_Site_global_ABneutral_CG_estimates.Rdata`, `Kn-0_Site_global_ABneutral_CG_estimates.Rdata`, `Mt-0_Site_global_ABneutral_CG_estimates.Rdata`, `Tsu-0_Site_global_ABneutral_CG_estimates.Rdata` | mCG divergence in 4 MA accessions |
| 1d | `plot_Figure1D.R` | `gain_loss-rate_MA-accessions_used_2026-01.txt` | Gain/loss rates in 4 MA accessions |
| 1e | `plot_Figure1E.R` | `new_observed_methy_level_using_proportion_all_samples_4_MA-accession_2025-05.txt`, `ABneutral_Boot_CG_estimates_Col0_global_2000.Rdata`, `ABneutral_Boot_CG_estimates_Kn0_global_2000.Rdata`, `ABneutral_Boot_CG_estimates_Mt0_global_2000.Rdata`, `ABneutral_Boot_CG_estimates_Tsu0_global_2000.Rdata` | Observed vs predicted mCG levels |
| 1g | `plot_Figure1G.R` | `known_makers_all_Beta_rate_add_genotype_matrix_without_ID_only_68RIL_lines_all_annotations_2023-11.csv`, `newrun_2023-11_68_RIL_all_annotations_alpha_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2023-11.csv` | Genetic map comparison |

### Figure 2: Mapping genetic mediators of methylation levels

| Figure | Script | Input data | Description |
|:-------|:-------|:-----------|:------------|
| 2a | `plot_Figure2A.R` | `RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt` | mCG levels across 68 RILs |
| 2b | `plot_Figure2B.R` | `RIL_371_all_sample_CG_new_methy_level_without_C_M_global_2026-01.txt` | mCG stability: G2 vs G9, L1 vs L2 |
| 2c | `plot_Figure2C.R` | `newrun_2024-11_68_RIL_global_gbM_gene_UM_TE_intergenic_teM_new_methy_lvl_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2024-11.csv`, `newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt` | mCG-level QTL (global + annotations) |
| 2d | `plot_Figure2D.R` | `newrun_2024-11_68_RIL_global_CHG_CWA_CHH_nonCWA_new_methy_lvl_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2024-11.csv`, `newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt` | Non-CG methylation QTL |
| 2e | `plot_Figure2E.R` | `68RIL_global_per5kb_qtl_peak_marker_with_pos.txt` | Cis-trans 5-kb window QTL map |
| 2f | `plot_Figure2F.R` | `68RIL_global_per5kb_qtl_peak_marker_chr1_with_trans_pos.txt`, `68RIL_global_per5kb_qtl_peak_marker_chr2_with_trans_pos.txt`, `gbM_gene_anotation_extract_Arabidopsis.bed`, `teM_gene_anotation_extract_Arabidopsis.bed`, `UM_gene_anotation_extract_Arabidopsis.bed`, `gene_anotation_extract_Arabidopsis.bed`, `other_gene_execpt_gbM_teM_UM_TAIR10_Arabidopsis_whole_gene.bed`, `TAIR10_TE.bed`, `TAIR10_whole_chr_position.txt`, `intergenic.gff3` | Feature enrichment in trans-bands |
| 2g | `plot_Figure2G.R` | `global_methy_lvl_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt` | Haplotype effects on global mCG |
| 2h | `plot_Figure2H.R` | `gbM_methy_lvl_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt` | Haplotype effects on gbM mCG |
| 2i | `plot_Figure2I.R` | `global_methy_lvl_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_only_G2_only_G9_2026-01.txt` | G2 vs G9 by Chr1 QTL haplotype |

### Figure 3: Mapping genetic modifiers of CG epimutation rates

| Figure | Script | Input data | Description |
|:-------|:-------|:-----------|:------------|
| 3a | `plot_Figure3A.R` | `Alpha_Beta_rate__70_RIL_lines.txt` | Gain/loss rates across 68 RILs |
| 3b | `plot_Figure3B.R` | `RIL15_Site_global_ABneutral_CG_estimates.Rdata`, `RIL50_Site_global_ABneutral_CG_estimates.Rdata`, `RIL58_Site_global_ABneutral_CG_estimates.Rdata`, `RIL149_Site_global_ABneutral_CG_estimates.Rdata`, `RIL150_Site_global_ABneutral_CG_estimates.Rdata`, `RIL155_Site_global_ABneutral_CG_estimates.Rdata`, `RIL162_Site_global_ABneutral_CG_estimates.Rdata`, `RIL182_Site_global_ABneutral_CG_estimates.Rdata` | Divergence in high- vs low-rate RILs |
| 3c | `plot_Figure3C.R` | `newrun_2023-11_68_RIL_all_annotations_add_intergeneic_alpha_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2025-04.csv`, `newrun_2023-11_68_RIL_all_annotations_add_intergeneic_beta_rate_without_ID_SNP_markers_after_using_SNP_sliding_windows_cleaning_2025-04.csv`, `newrun_2023-11_68RIL_with_position_sliding_windows_200SNP_5SNP_step_for_SNP_matrix_with_counts_0_2_NA_as_equal_matrix.txt` | Alpha/beta rate QTL (global + annotations) |
| 3d | `plot_Figure3D.R` | `global_gain_loss_rate_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt` | Haplotype effects on global rates |
| 3e | `plot_Figure3E.R` | `gbM_gain_loss_rate_Lerseg_withoutchr1_vs_Cviseg_withoutchr1_RIL_merge_peakLer_peakCvi_2026-01.txt` | Haplotype effects on gbM rates |

### Figure 4: VIM2/VIM4 deletion perturbs CG methylation maintenance

| Figure | Script | Input data | Description |
|:-------|:-------|:-----------|:------------|
| 4b | `plot_Figure4B.R` | `AT1TE80755_methy_lvl_vs_expression_VIM4_2026-01.txt` | VIM2/4 expression vs TE methylation |
| 4c | `plot_Figure4C.R` | `expression_VIM4_vs_global_methy_lvl_2026-01.txt` | VIM2/4 expression vs global mCG |
| 4d | `plot_Figure4D.R` | `expression_VIM4_vs_global_rate_2026-01.txt` | VIM2/4 expression vs epimutation rate |
| 4f | `plot_Figure4F.R` | Hardcoded in script (qPCR values) | qPCR of VIM2/4 |
| 4g | `plot_Figure4G.R` | `preTableSNew_expression_level_across_67RILs.txt` | Co-expression heatmaps |
| 4h | `plot_Figure4H.R` | `preTableSNew_expression_level_across_67RILs.txt` | DEG barplot by peak genotype |

### Figure 5: Epigenomic, chromatin, and phenotypic consequences

| Figure | Script | Input data | Description |
|:-------|:-------|:-----------|:------------|
| 5a | `plot_Figure5A.R` | `TableSNew_Cvi_accssions_global_mCG_level_deletion.txt` | mCG in Cvi accessions +/- deletion |
| 5b | `plot_Figure5B.R` | `TableSNew_Cvi_accssions_global_mCG_level_deletion.txt` | mCG variance in Cvi accessions |
| 5c | `plot_Figure5C.R` | `TableSNew_SNPs_within_deletion_global_mCG_level_across_799natural_accessions.txt` | mCG by SNP genotype in deletion |
| 5d | `plot_Figure5D.R` | `natural_accessions_global_Methylation_level_variance_across_SNPs_within_deletion_by_major_minor.txt` | mCG variance by SNP genotype |
| 5e | `plot_Figure5E.R` | `RIL149_ATAC_ABneutral_CG_estimates.Rdata`, `RIL50_ATAC_ABneutral_CG_estimates.Rdata` | ATAC-seq divergence: RIL50 vs RIL149 |
| 5f | `plot_Figure5F.R` | `pheno_peak_status_68RIL_312phenos_prepare_lm_1_2623_0as-1_2as1.csv` | Phenotypic variance across 312 traits |
| 5g | `plot_Figure5G.R` | `phenotype_QTL_68_RIL_312_traits.txt` | pQTL associations at Chr1 peak |

Helper scripts:

| Script | Used by |
|:-------|:--------|
| `plotFits-overlay_interceptZero_formal_2021-09_using_genotypes_than_stress_adjust_2023-11.R` | Fig. 1c |
| `plotfits_divergence_adjust_fig3b_two_color_2025-05.R` | Fig. 3b |
| `plotFits-overlay_interceptZero_formal_2021-09_using_genotypes_than_stress.R` | Fig. 5e |

---

## Workflow Documentation

Detailed pipeline documentation is provided in `workflow/`.

| Directory | Purpose | Key tools |
|:----------|:--------|:----------|
| `workflow/WGBS/` | WGBS preprocessing, Methimpute calling, AlphaBeta rate estimation, QTL mapping | MethylStar, AlphaBeta, R/qtl |
| `workflow/ATAC/` | ATAC-seq alignment, peak calling, accessibility divergence analysis | BWA, MACS2, AlphaBeta |
| `workflow/RNAseq/` | RNA-seq alignment, expression quantification |  HTSeq |

## Workflow-specific README files

Detailed README files for each sequencing workflow are provided in the corresponding subdirectories:

- `workflow/WGBS/README.md`: WGBS preprocessing, methylation-state calling, AlphaBeta epimutation-rate estimation, and QTL mapping.
- `workflow/ATAC/README.md`: ATAC-seq quality control, trimming, alignment, duplicate removal, peak calling, and accessibility-divergence analysis.
- `workflow/RNAseq/README.md`: RNA-seq preprocessing and expression quantification.

The top-level README provides an overview of the repository and maps the main manuscript figures to scripts and processed input tables. The workflow-specific README files provide more detailed instructions for each data type.

### WGBS workflow

The WGBS workflow is documented in:

```text
workflow/WGBS/
```

This workflow covers:

### Step 1: WGBS preprocessing (MethylStar)

Raw WGBS reads were processed using [MethylStar](https://github.com/jlab-code/MethylStar): trimming, Bismark alignment (TAIR10), methylation extraction, and Methimpute calling. Only calls with `posteriorMax >= 0.99` were retained. See `workflow/WGBS/`.

### Step 2: Epimutation rate estimation (AlphaBeta)

Per-RIL CG gain (alpha) and loss (beta) rates were estimated using [AlphaBeta](https://github.com/jlab-code/AlphaBeta) with the neutral model, F-test validation, and 1000 bootstrap replicates. See `workflow/WGBS/`.

### Step 3: QTL mapping (R/qtl)

QTL scans were performed using R/qtl (Haley-Knott regression). The genetic map was re-estimated with `error.prob = 1e-7`. Genome-wide significance thresholds were determined from 1000 permutations. LOD confidence intervals were obtained with `lodint(drop = 2)`.

| QTL scan | Phenotype | Figure |
|:---------|:----------|:-------|
| Methylation-level QTL | Steady-state mCG levels (global + annotations) | Fig. 2c |
| Non-CG methylation QTL | CHG, CHH, CWA, non-CWA levels | Fig. 2d |
| Alpha-rate QTL | Methylation gain rates (global + annotations) | Fig. 3c |
| Beta-rate QTL | Methylation loss rates (global + annotations) | Fig. 3c |

A detailed R Markdown workflow notebook is provided in `workflow/WGBS/`.


### ATAC-seq workflow

The ATAC-seq workflow is documented in:

```text
workflow/ATAC/
```

This workflow covers read quality control, adapter trimming, alignment, duplicate removal, BAM-to-BED conversion, peak calling, and downstream accessibility-divergence analysis.

Main purpose:

- Process ATAC-seq reads.
- Align reads using BWA.
- Remove duplicates.
- Convert BAM to BED.
- Call peaks using MACS2.
- Generate accessibility-state matrices.
- Estimate accessibility divergence across generations.

Main downstream figure panel:

| Figure | Related output |
|---|---|
| Fig. 5e | ATAC-seq accessibility divergence in selected high- and low-rate RILs |



### RNA-seq workflow

The RNA-seq workflow is documented in:

```text
workflow/RNAseq/
```

This workflow covers RNA-seq preprocessing and expression quantifications.


Main purpose:

- Process RNA-seq data.
- Quantify expression.

Main downstream figure panels:

| Figure | Related output |
|---|---|
| Fig. 4b-d | Expression of VIM4 versus methylation and epimutation-rate phenotypes |
| Fig. 4g | Co-expression heatmaps among VIM genes and methylation regulators |
| Fig. 4h | Expression differences between peak genotypes |

---



## Key software and tools

| Tool or package | Purpose |
|---|---|
| MethylStar | WGBS preprocessing and Methimpute-based methylation-state calling |
| Bismark | Bisulfite read alignment and methylation extraction |
| FastQC | FASTQ quality control |
| Trimmomatic | Adapter and quality trimming |
| Samtools | BAM processing |
| Bedtools | Genomic interval operations |
| BWA-MEM | ATAC-seq read alignment |
| Methimpute | Cytosine methylation-state classification |
| AlphaBeta | Pedigree-based estimation of methylation gain and loss rates |
| R/qtl | QTL mapping using Haley-Knott regression |
| BWA | ATAC-seq read alignment |
| MACS2 | ATAC-seq peak calling |
| ggplot2, data.table, dplyr | Data processing and visualization |

---


## Reproducing figure outputs

Most figure panels can be regenerated by running the corresponding R scripts from the repository root after placing the processed input files in `data/`.

Example:

```bash
Rscript scripts/plot_Figure2C.R
Rscript scripts/plot_Figure3C.R
```

Output PDFs are written to the `figures/` directory.

Some panels are schematic or composite figure panels, including experimental-design diagrams, genome-browser snapshots, and PCR-gel panels. These do not have a single standalone plotting script but are based on the processed outputs and source data documented above.


## Notes
- This repository is organized for the June 2025 preprint version.
- File paths in scripts may need adjustment for local use.
- Large raw sequencing files are not stored directly in this repository and should be obtained from GEO or the linked data archive.
- All analyses restricted to Chr1-5.
- In legacy file names, `UM` = lowly methylated genes (`LM` in manuscript).

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


