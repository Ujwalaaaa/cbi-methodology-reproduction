# Reproduction of the CBI Methodology on a Synthetic Ovarian-Cancer Benchmark

## Overview

This repository contains the implementation of the Chemotherapy Benefit
Index (CBI) methodology described in Ma et al. (2025), applied to a
synthetic ovarian-cancer-like benchmark with known ground truth.

This is a **methodology reproduction**, not a numerical reproduction.
The paper's reported numbers depend on the original TCGA cohort and its
under-documented preprocessing, which was not reproducible.

## Pipeline

Five-stage pipeline, faithful to the paper:

1. Correlation network (WGCNA, weight > 0.85)
2. Module assignment (5 biological themes)
3. Bayesian network + consensus graph (bnlearn, R=100, hc, bic-g)
4. Node selection (total degree >= 3)
5. Multivariable Cox regression (P < 0.05) + CBI construction
6. Random Forest classification (paper hyperparameters)

## Repository Structure

    code/
      01_wgcna_network_qc.R          WGCNA network and power scan
      02_blind_cbi_paper_pipeline.R  BN + consensus graph + node selection
      03_stage3a_cox_cbi.R           Cox regression + CBI + tertiles
      04_stage3b_rf_seed_robust.py   Random Forest across 10 seeds
      05_auc_evaluation.py           Internal + external testing AUC

    metadata/
      predefined_module_map.csv      49 genes assigned to 5 modules

## Key Results

| Metric | Value |
|---|---|
| Candidate genes entering Cox | 30 |
| Cox-selected genes | 7 |
| Survival genes recovered | 5 / 10 |
| Decoys eliminated | 8 / 8 |
| Internal testing AUC (mean, 10 seeds) | 0.9564 |
| External testing AUC cohortA (N=100) | 0.9518 |
| External testing AUC cohortB (N=150) | 0.9511 |
| External testing AUC cohortC (N=200) | 0.9745 |

## How to Run

### Requirements

- R 4.6.1 with packages: WGCNA, bnlearn, survival
- Python 3.10+ with packages: numpy, pandas, scikit-learn

### Install dependencies

    Rscript -e 'install.packages(c("WGCNA","bnlearn","survival"))'
    pip install -r requirements.txt

### Execution order

    Rscript code/01_wgcna_network_qc.R
    Rscript code/02_blind_cbi_paper_pipeline.R
    Rscript code/03_stage3a_cox_cbi.R
    python  code/04_stage3b_rf_seed_robust.py
    python  code/05_auc_evaluation.py

### Configure paths

The scripts contain a ROOT variable that points to the benchmark
directory. Change it to match your local installation.

## Data

The pipeline uses a synthetic ovarian-cancer-like benchmark, not the
original TCGA cohort.

## Documented Deviations

1. WGCNA power = 0.30 (paper does not specify its power)
2. Consensus graph threshold = 0.5 (paper does not specify)
3. CBI positive shift (required because standardised Risk Score crosses zero)
4. GO enrichment substituted with predefined module map (synthetic genes
   lack real Gene Ontology annotations)

## Reference

Ma S, Zhou L, et al. (2025). Identification of a novel chemotherapy benefit
index for patients with advanced ovarian cancer based on Bayesian network
analysis. PLoS One 20(5): e0322130.
