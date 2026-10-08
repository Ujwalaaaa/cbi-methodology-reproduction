"""
stage3b_rf_seed_robust.py

Random Forest AUC across 10 random seeds, using the same
gene set and tertile split as the primary run.

Purpose: report mean / median / min / max / std of held-out AUC
         to show the result is not a lucky single split.

Paper settings (Section 2.3 of the CBI paper):
    class_weight = balanced
    criterion    = entropy
    n_estimators = 230
    max_features = 10
    min_samples_leaf = 3
    max_depth = 5
    train/test = 70/30

Inputs:
    results/stage3_cox_cbi_rf/tertiles.csv
    results/stage3_cox_cbi_rf/cox_selected_genes.csv
    blinded/synthetic_expression.csv

Outputs:
    results/stage3_cox_cbi_rf/rf_auc_seed_robust.csv
    results/stage3_cox_cbi_rf/rf_auc_seed_summary.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


# ------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------
ROOT  = Path("C:/Users/MARAO/Downloads/benchmark/benchmark")
OUT   = ROOT / "results" / "stage3_cox_cbi_rf"
BLIND = ROOT / "blinded"


# ------------------------------------------------------------
# 1. LOAD
# ------------------------------------------------------------
expr = pd.read_csv(BLIND / "synthetic_expression.csv", index_col="patient_id")
tert = pd.read_csv(OUT / "tertiles.csv").set_index("patient_id")
sel  = pd.read_csv(OUT / "cox_selected_genes.csv")

print(f"Expression matrix: {expr.shape}")
print(f"Tertiles:          {tert.shape}")
print(f"Cox-selected genes: {sel['gene'].tolist()}")


# ------------------------------------------------------------
# 2. ALIGN
# ------------------------------------------------------------
common = expr.index.intersection(tert.index)
expr = expr.loc[common]
tert = tert.loc[common]


# ------------------------------------------------------------
# 3. CBI-HIGH vs CBI-LOW
# ------------------------------------------------------------
mask = tert["tertile"].isin(["low", "high"])
X = expr.loc[mask, sel["gene"].tolist()].values
y = (tert.loc[mask, "tertile"] == "high").astype(int).values

print(f"\nN high: {int(y.sum())}, N low: {int((1 - y).sum())}")


# ------------------------------------------------------------
# 4. 10-SEED RF
# ------------------------------------------------------------
N_SEEDS = 10
aucs = []

for seed in range(N_SEEDS):
    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y,
        test_size=0.30,
        random_state=seed,
        stratify=y,
    )

    clf = RandomForestClassifier(
        class_weight="balanced",
        criterion="entropy",
        n_estimators=230,
        max_features=10,
        min_samples_leaf=3,
        max_depth=5,
        random_state=seed,
        n_jobs=-1,
    )
    clf.fit(X_tr, y_tr)

    proba = clf.predict_proba(X_te)[:, 1]
    auc_seed = roc_auc_score(y_te, proba)
    aucs.append(auc_seed)

    print(f"  seed {seed:2d}: AUC = {auc_seed:.4f}")


# ------------------------------------------------------------
# 5. SUMMARY
# ------------------------------------------------------------
aucs_arr = np.array(aucs)

print("\n================ SEED ROBUSTNESS ================")
print(f"Mean AUC:   {np.mean(aucs_arr):.4f}")
print(f"Median AUC: {np.median(aucs_arr):.4f}")
print(f"Min AUC:    {np.min(aucs_arr):.4f}")
print(f"Max AUC:    {np.max(aucs_arr):.4f}")
print(f"Std AUC:    {np.std(aucs_arr, ddof=1):.4f}")
print(f"Seeds > 0.90: {int(np.sum(aucs_arr > 0.90))}/{N_SEEDS}")
print(f"All: {[round(a, 4) for a in aucs]}")


# ------------------------------------------------------------
# 6. WRITE OUTPUTS
# ------------------------------------------------------------
pd.DataFrame({
    "seed": list(range(N_SEEDS)),
    "auc":  aucs,
}).to_csv(OUT / "rf_auc_seed_robust.csv", index=False)

summary = pd.DataFrame({
    "metric": [
        "mean_auc",
        "median_auc",
        "min_auc",
        "max_auc",
        "std_auc",
        "n_seeds",
        "n_seeds_above_0.90",
    ],
    "value": [
        round(float(np.mean(aucs_arr)), 4),
        round(float(np.median(aucs_arr)), 4),
        round(float(np.min(aucs_arr)), 4),
        round(float(np.max(aucs_arr)), 4),
        round(float(np.std(aucs_arr, ddof=1)), 4),
        N_SEEDS,
        int(np.sum(aucs_arr > 0.90)),
    ],
})
summary.to_csv(OUT / "rf_auc_seed_summary.csv", index=False)

print(f"\nWrote: {OUT / 'rf_auc_seed_robust.csv'}")
print(f"Wrote: {OUT / 'rf_auc_seed_summary.csv'}")