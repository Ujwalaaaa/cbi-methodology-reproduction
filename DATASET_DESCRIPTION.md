# Dataset Description

## Type

Synthetic benchmark. This is **not** the original TCGA ovarian cancer
cohort. It is a controlled synthetic dataset generated to mirror the
paper's structure.

## Structure

| Property | Value |
|---|---|
| Patients | 146 |
| Genes | 3,000 |
| Favorable / poor split | 76 / 70 |
| Structured genes | 49 |
| Survival genes | 10 |
| Decoy genes | 8 |
| Network-other genes | 31 |
| Hidden confounders | 4 (H1-H4) |
| Biological modules | 5 |
| External cohorts | A (100), B (150), C (200) |

## Files

| File | Content |
|---|---|
| blinded/synthetic_expression.csv | 146 x 3001 |
| blinded/synthetic_clinical.csv | 146 x 5 |
| metadata/predefined_module_map.csv | 49 genes to 5 modules |
| ground_truth/ground_truth.json | planted answer (audit only) |
| external_validation/cohortA-C | external validation cohorts |

## Why Synthetic

The original TCGA cohort could not be reproduced because the paper does
not document its preprocessing pipeline or its sample-selection criteria
at the level required. The synthetic benchmark provides a controlled
setting with known ground truth, allowing the pipeline implementation to
be validated.

## Calibration

The benchmark was calibrated before freezing so the planted signal is
recoverable:

| Parameter | Value |
|---|---|
| Coefficient scale kappa | 0.95 |
| Weibull shape | 1.4 |
| Median baseline | 20 months |
| Measurement noise (latent SD) | 0.10 - 0.22 |
| Named-edge correlations | 0.80 - 0.92 |
| Decoy correlations | 0.94 - 0.97 |

## Reproducibility

Generator: benchmark v1.0.0
Master seed: 20261005
Byte-identical regeneration verified via SHA256 checksums.
