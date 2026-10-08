# ============================================================
# STAGE 3A: COX -> CBI -> TERTILES  (R, as in the paper)
# ============================================================
rm(list = ls()); gc()
set.seed(20261005)

suppressPackageStartupMessages({ library(survival) })

ROOT  <- "C:/Users/MARAO/Downloads/benchmark/benchmark"
BLIND <- file.path(ROOT, "blinded")
CAND_FILE <- file.path(ROOT, "results", "pooled_49_module_audit",
                       "module_union_candidates_cpdag_TRUE.csv")
OUT <- file.path(ROOT, "results", "stage3_cox_cbi_rf")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

expr_raw <- read.csv(file.path(BLIND, "synthetic_expression.csv"),
                     check.names = FALSE)
clin <- read.csv(file.path(BLIND, "synthetic_clinical.csv"),
                 check.names = FALSE)
cand <- read.csv(CAND_FILE, check.names = FALSE)

expr_raw <- expr_raw[match(clin$patient_id, expr_raw$patient_id), ]
stopifnot(all(expr_raw$patient_id == clin$patient_id))

gene_names <- setdiff(names(expr_raw), "patient_id")
expr <- as.matrix(expr_raw[, gene_names, drop = FALSE])
rownames(expr) <- expr_raw$patient_id

# Paper-faithful 76/70 split (Version A)
fav_idx  <- clin$prognosis_controlled_A == "favorable"
poor_idx <- clin$prognosis_controlled_A == "poor"
cat("Favorable:", sum(fav_idx), " Poor:", sum(poor_idx), "\n")

available <- intersect(cand$gene, colnames(expr))
cat("Candidate genes available:", length(available), "\n")
X <- scale(expr[, available, drop = FALSE])

cox_df <- data.frame(time = clin$os_months, event = clin$os_event,
                     X, check.names = FALSE)
cox_fit <- coxph(Surv(time, event) ~ ., data = cox_df)
s <- summary(cox_fit)

coef_table <- data.frame(
  gene    = rownames(s$coefficients),
  coef    = s$coefficients[, "coef"],
  hr      = s$coefficients[, "exp(coef)"],
  p_value = s$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)
write.csv(coef_table, file.path(OUT, "cox_results.csv"), row.names = FALSE)

sel <- coef_table$gene[coef_table$p_value < 0.05]
cat("\nCox-selected (p<0.05):", length(sel), "\n"); print(sel)
if (length(sel) == 0) stop("No genes selected by Cox.")

beta <- coef_table$coef[match(sel, coef_table$gene)]
names(beta) <- sel
X_sel <- X[, sel, drop = FALSE]

rs <- as.numeric(X_sel %*% beta)
rs_shift <- rs - min(rs) + 0.05 * (max(rs) - min(rs))
cbi <- 1 / rs_shift

cbi_df <- data.frame(patient_id = rownames(X_sel),
                     risk_score = rs, rs_shifted = rs_shift,
                     cbi = cbi)
write.csv(cbi_df, file.path(OUT, "cbi.csv"), row.names = FALSE)

tert <- cut(cbi, breaks = quantile(cbi, probs = seq(0, 1, 1/3)),
            labels = c("low","medium","high"), include.lowest = TRUE)
tert_df <- data.frame(patient_id = rownames(X_sel),
                      cbi = cbi, tertile = as.character(tert))
write.csv(tert_df, file.path(OUT, "tertiles.csv"), row.names = FALSE)
cat("\nTertile counts:\n"); print(table(tert))

# Save selected genes and beta for the Python RF step
write.csv(data.frame(gene = sel, beta = as.numeric(beta)),
          file.path(OUT, "cox_selected_genes.csv"), row.names = FALSE)
cat("\nStage 3A complete.\n")