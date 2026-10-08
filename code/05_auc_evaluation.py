"""
auc_test_only.py

Reports testing AUC only (mean and min) across
internal and three external cohorts.

No training AUC. No max. No std.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


ROOT = Path("C:/Users/MARAO/Downloads/benchmark/benchmark")
OUT  = ROOT / "results" / "auc_evaluation"
OUT.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Load frozen model
# ------------------------------------------------------------
cox_sel = pd.read_csv(ROOT / "results" / "stage3_cox_cbi_rf" / "cox_selected_genes.csv")
genes = cox_sel["gene"].tolist()
betas = cox_sel["beta"].values

print(f"Frozen model: {len(genes)} genes")
print(f"  {genes}")
print()


# ------------------------------------------------------------
# Compute CBI tertiles
# ------------------------------------------------------------
def compute_cbi_tertiles(expr_df, genes, betas):
    X = expr_df[genes].copy()
    X = (X - X.mean()) / (X.std(ddof=1) + 1e-12)
    rs = X.values @ betas
    rs_shift = rs - rs.min() + 0.05 * (rs.max() - rs.min())
    cbi = 1.0 / rs_shift
    tert = pd.qcut(cbi, q=3, labels=["low", "medium", "high"])
    return X, tert


# ------------------------------------------------------------
# Testing AUC for one seed
# ------------------------------------------------------------
def test_auc(X, tert, seed=42, test_size=0.30):
    mask = np.asarray(tert.isin(["low", "high"]), dtype=bool)
    X_rf = X.values[mask]
    y_rf = (np.asarray(tert)[mask] == "high").astype(int)

    if y_rf.sum() < 5 or (1 - y_rf).sum() < 5:
        return None

    X_tr, X_te, y_tr, y_te = train_test_split(
        X_rf, y_rf,
        test_size=test_size,
        random_state=seed,
        stratify=y_rf,
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

    return roc_auc_score(y_te, clf.predict_proba(X_te)[:, 1])


# ------------------------------------------------------------
# Evaluate one cohort across 10 seeds
# ------------------------------------------------------------
def evaluate_cohort(expr_df, label):
    X, tert = compute_cbi_tertiles(expr_df, genes, betas)

    aucs = []
    for seed in range(10):
        a = test_auc(X, tert, seed=seed)
        if a is not None:
            aucs.append(a)

    return {
        "cohort":  label,
        "n":       int(len(expr_df)),
        "auc_mean": round(float(np.mean(aucs)), 4),
        "auc_min":  round(float(np.min(aucs)),  4),
    }


# ------------------------------------------------------------
# Run all cohorts
# ------------------------------------------------------------
rows = []

expr_int = pd.read_csv(
    ROOT / "blinded" / "synthetic_expression.csv",
    index_col="patient_id",
)
rows.append(evaluate_cohort(expr_int, "internal"))

for cohort in ["cohortA", "cohortB", "cohortC"]:
    expr_ext = pd.read_csv(
        ROOT / "external_validation" / cohort / "synthetic_expression.csv",
        index_col="patient_id",
    )
    rows.append(evaluate_cohort(expr_ext, cohort))


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "auc_test_only.csv", index=False)

print("=" * 60)
print("TESTING AUC ONLY")
print("=" * 60)
print(summary.to_string(index=False))
print()
print(f"Saved: {OUT / 'auc_test_only.csv'}")