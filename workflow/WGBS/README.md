# WGBS Workflow: Methylation-State Calling, Epimutation-Rate Estimation, and QTL Mapping

This directory documents the WGBS analysis pipeline from methylome preprocessing to QTL mapping of methylation levels and epimutation rates.

This directory documents the WGBS-based analysis workflow used for the bioRxiv preprint:

**Mapping genetic modifiers of epimutation rates reveals a punctuated-equilibrium model of CG methylome evolution**  
Zhang et al., bioRxiv 2025  
DOI: `10.1101/2025.06.14.659605`


## Pipeline Overview

| Step | Description | Tool |
|:-----|:------------|:-----|
| Step 1 | WGBS alignment, methylation extraction, Methimpute calling | MethylStar / Bismark |
| Step 2 | Pedigree-based epimutation rate estimation | AlphaBeta |
| Step 3 | QTL mapping of mCG levels and epimutation rates | R/qtl |

---
## Workflow details

### Step 1: WGBS preprocessing and methylation-state calling

Raw WGBS reads were processed using MethylStar, which includes Bismark alignment and Methimpute methylation-state calling.

Relevant external software:

- MethylStar: `https://github.com/jlab-code/MethylStar`
- Bismark
- Bowtie2
- Methimpute

Methimpute assigns each cytosine to one of three methylation states:

- `U`: unmethylated
- `I`: intermediate
- `M`: methylated

Only high-confidence calls were retained for downstream analyses.

---

### Step 2: AlphaBeta rate estimation

CG methylation gain and loss rates were estimated with AlphaBeta using the pedigree structure of each mutation accumulation mapping line.

Main scripts:

```text
run_AB_global_rate_check_samples_RIL4.R
function_run_alphabeta_nstarts_2000_CG_progenitor_intermediate_without_selection_using_newp0uu_rmCM_2023-11.R
newbuildPedigree_stand_alone_rmCM_2023_11.R
convertDMATRIX.R
time.R
```

Main analysis settings:

| Parameter | Value |
|---|---|
| Cytosine context | `CG` |
| Methimpute posterior filter | `posteriorMax >= 0.99` |
| AlphaBeta model | `ABneutral` |
| Starting points | `Nstarts = 2000` |
| Bootstrap replicates | `Nboot = 1000` |
| Model comparison | `ABneutral` versus `ABnull` using F-test |

The example script `run_AB_global_rate_check_samples_RIL4.R` shows how the workflow was run for one RIL. The same logic was applied across RILs and annotation classes.

---

### Step 3: QTL mapping

QTL scans were performed using the R package `qtl`.

Main script:

```text
plot_Figure3C.R
```

Main analysis settings:

| Parameter | Value |
|---|---|
| Genotype coding | `0 = Ler`, `2 = Cvi` |
| QTL method | Haley-Knott regression (`method = "hk"`) |
| Genetic-map estimation | `est.map(error.prob = 1e-7)` |
| Genotype probabilities | `calc.genoprob(step = 1)` |
| Permutations | `n.perm = 1000` |
| Confidence interval | `lodint(drop = 2)` |

QTL mapping was performed for:

- steady-state mCG levels;
- non-CG methylation levels;
- AlphaBeta-estimated CG gain rates;
- AlphaBeta-estimated CG loss rates.

The main downstream figures are:

| Figure | Analysis |
|---|---|
| Fig. 2c | mCG-level QTL scan |
| Fig. 2d | non-CG methylation QTL scan |
| Fig. 3c | alpha-rate and beta-rate QTL scans |

---
## Files

| File | Description |
|:-----|:------------|
| `WGBS_rate_QTL_workflow_pipeline_final.Rmd` | R Markdown notebook documenting the full workflow (Steps 1-3) |
| `run_AB_global_rate_check_samples_RIL4.R` | Example script: run AlphaBeta for RIL4 global CG rates |
| `function_run_alphabeta_nstarts_2000_CG_progenitor_intermediate_without_selection_using_newp0uu_rmCM_2023-11.R` | AlphaBeta wrapper function (ABneutral + ABnull + bootstrap) |
| `newbuildPedigree_stand_alone_rmCM_2023_11.R` | Modified `buildPedigree` function (removes chloroplast/mitochondria) |
| `convertDMATRIX.R` | Helper: convert divergence matrix for AlphaBeta |
| `time.R` | Helper: timing utility |
| `plot_Figure3C.R` | QTL LOD curve plotting (Fig. 3c, alpha and beta rates) |


---

## Input data

The workflow expects the following input types.

| Input | Description |
|---|---|
| Methimpute methylation-state files | Per-sample cytosine-level methylation-state files generated from WGBS data |
| Node list | Table listing sampled individuals in each pedigree, generation information, and corresponding methylation-state files |
| Edge list | Table describing parent-offspring relationships in the MA pedigree |
| QTL cross files | R/qtl-formatted CSV files containing genotype markers and methylation or rate phenotypes |
| Annotation files | Gene, gbM, teM, TE, intergenic, and other genomic-feature annotations |

High-confidence Methimpute calls with `posteriorMax >= 0.99` were used for AlphaBeta rate estimation.

---


## Expected output

| Output | Description |
|---|---|
| `ABneutral_CG_estimates_*.Rdata` | AlphaBeta neutral-model output |
| `ABnull_CG_estimates_*.Rdata` | AlphaBeta null-model output |
| `ABneutral_Boot_CG_estimates_*.Rdata` | Bootstrap output for AlphaBeta estimates |
| `ABneutral_estimatats_*.txt` | Estimated alpha and beta parameters |
| `out_ABneutral_ABnull_Ftest_p_value.txt` | F-test comparison between neutral and null models |
| QTL scan tables | R/qtl scan results for methylation levels, alpha rates, and beta rates |
| LOD curve PDFs | QTL plots used for preprint Fig. 2c and Fig. 3c |

---

## Requirements

- R (4.x) with packages: AlphaBeta, qtl, data.table, ggplot2, stringr, scales
- [MethylStar](https://github.com/jlab-code/MethylStar)
- [AlphaBeta](https://github.com/jlab-code/AlphaBeta)

## Usage

The R Markdown notebook can be knitted as-is (code blocks default to `eval = FALSE`). To rerun the analysis, set `run_code: true` in the YAML header and update paths in `params`.

To run AlphaBeta for a single RIL:

```r
source("function_run_alphabeta_nstarts_2000_CG_progenitor_intermediate_without_selection_using_newp0uu_rmCM_2023-11.R")
run.alphabeta.new(
  nodelist   = "nodelist_update_RIL4.txt",
  edelist    = "RILRaw_edgelist.txt",
  name       = "RIL4_global_rate",
  input.dir  = "path/to/node_edg/",
  output.dir = "path/to/output/"
)
```

---

## Notes

- Some scripts contain absolute paths from the original analysis environment. These paths should be updated before running the scripts on another machine.
- The GitHub repository contains workflow documentation and scripts. Large raw sequencing files should be downloaded from GEO or the linked data archive.
- Downstream analyses were restricted to the five nuclear chromosomes unless otherwise stated.

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

