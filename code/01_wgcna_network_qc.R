# ============================================================
# 01_wgcna_network_qc.R
#
# Paper-faithful upstream network reconstruction
# Paper:
# "Identification of a novel chemotherapy benefit index..."
#
# Purpose:
#   1. Load ONLY blinded expression + clinical data
#   2. Reproduce 76/70 prognosis split
#   3. Run WGCNA separately in both cohorts
#   4. Quantify strong gene-pairs at weight > 0.85
#   5. Compare against paper-reported:
#        Better prognosis = 4,415 pairs
#        Poor prognosis   = 4,050 pairs
#
# IMPORTANT:
#   We do NOT tune a WGCNA parameter to force 4,415/4,050.
#   The paper does not report enough information to justify that.
#
#   Instead, we perform an explicit parameter sensitivity QC.
#   The power scan is limited to 0.1–1.0 because integer powers
#   above 1 produce very few strong pairs on this dataset and
#   are not informative for the comparison.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(WGCNA)
})

allowWGCNAThreads()

# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

ROOT <- "C:/Users/MARAO/Downloads/benchmark/benchmark"

EXPR_FILE <- file.path(ROOT, "blinded", "synthetic_expression.csv")
CLIN_FILE <- file.path(ROOT, "blinded", "synthetic_clinical.csv")

OUT <- file.path(ROOT, "results", "wgcna_network_qc")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Paper targets
# ------------------------------------------------------------

PAPER_FAVORABLE_N <- 76
PAPER_POOR_N      <- 70

PAPER_FAVORABLE_EDGES <- 4415
PAPER_POOR_EDGES      <- 4050

WEIGHT_THRESHOLD <- 0.85

# ------------------------------------------------------------
# 2. Load blinded data
# ------------------------------------------------------------

expr_raw <- read.csv(EXPR_FILE, check.names = FALSE)
clinical <- read.csv(CLIN_FILE, check.names = FALSE)

cat("\n================ DATA QC ================\n")
cat("Expression dimensions:", nrow(expr_raw), "x", ncol(expr_raw), "\n")
cat("Clinical dimensions:",   nrow(clinical), "x", ncol(clinical), "\n")

stopifnot(
  "patient_id" %in% names(expr_raw),
  "patient_id" %in% names(clinical)
)

# ------------------------------------------------------------
# 3. Match patients
# ------------------------------------------------------------

common_ids <- intersect(expr_raw$patient_id, clinical$patient_id)

stopifnot(length(common_ids) == nrow(clinical))
stopifnot(length(common_ids) == nrow(expr_raw))

expr_raw <- expr_raw[match(common_ids, expr_raw$patient_id), ]
clinical <- clinical[match(common_ids, clinical$patient_id), ]

stopifnot(identical(
  as.character(expr_raw$patient_id),
  as.character(clinical$patient_id)
))

# ------------------------------------------------------------
# 4. Expression matrix
# ------------------------------------------------------------

gene_names <- setdiff(colnames(expr_raw), "patient_id")
expr <- as.matrix(expr_raw[, gene_names, drop = FALSE])
storage.mode(expr) <- "numeric"

cat("Patients:", nrow(expr), "\n")
cat("Genes:",    ncol(expr), "\n")
cat("Missing expression values:", sum(is.na(expr)), "\n")

stopifnot(!anyNA(expr))

# ------------------------------------------------------------
# 5. Reproduce prognosis split
# ------------------------------------------------------------

stopifnot("prognosis_controlled_A" %in% names(clinical))

prognosis <- clinical$prognosis_controlled_A

cat("\n================ PROGNOSIS QC ================\n")
print(table(prognosis, useNA = "ifany"))

fav_idx  <- which(prognosis == "favorable")
poor_idx <- which(prognosis == "poor")

cat("Favorable:", length(fav_idx), "\n")
cat("Poor:",      length(poor_idx), "\n")

stopifnot(length(fav_idx)  == PAPER_FAVORABLE_N)
stopifnot(length(poor_idx) == PAPER_POOR_N)

# ------------------------------------------------------------
# 6. WGCNA input matrices
# ------------------------------------------------------------

datExpr_fav  <- expr[fav_idx,  , drop = FALSE]
datExpr_poor <- expr[poor_idx, , drop = FALSE]

gqc_fav  <- goodSamplesGenes(datExpr_fav,  verbose = 0)
gqc_poor <- goodSamplesGenes(datExpr_poor, verbose = 0)

cat("\nFavorable good samples:", gqc_fav$allOK,  "\n")
cat("Poor good samples:",       gqc_poor$allOK, "\n")

stopifnot(gqc_fav$allOK)
stopifnot(gqc_poor$allOK)

# ------------------------------------------------------------
# 7. Helper: count upper-triangle weights > threshold
# ------------------------------------------------------------

count_strong_pairs <- function(A, threshold = 0.85) {
  A <- as.matrix(A)
  diag(A) <- NA
  vals <- A[upper.tri(A, diag = FALSE)]
  vals <- vals[is.finite(vals)]
  strong <- vals > threshold
  list(
    n_pairs    = length(vals),
    n_strong   = sum(strong),
    proportion = mean(strong),
    max_weight = max(vals),
    median_weight = median(vals),
    q95_weight = as.numeric(quantile(vals, 0.95))
  )
}

# ------------------------------------------------------------
# 8. Baseline raw correlation diagnostic
# ------------------------------------------------------------

cat("\n================ RAW CORRELATION DIAGNOSTIC ================\n")

R_fav  <- cor(datExpr_fav,  method = "pearson", use = "pairwise.complete.obs")
R_poor <- cor(datExpr_poor, method = "pearson", use = "pairwise.complete.obs")

raw_fav  <- count_strong_pairs(abs(R_fav),  WEIGHT_THRESHOLD)
raw_poor <- count_strong_pairs(abs(R_poor), WEIGHT_THRESHOLD)

cat("Raw |cor| > 0.85 favorable:", raw_fav$n_strong,  "\n")
cat("Raw |cor| > 0.85 poor:",       raw_poor$n_strong, "\n")

# ------------------------------------------------------------
# 9. WGCNA adjacency sensitivity analysis
# ------------------------------------------------------------
#
# The paper says "weight > 0.85" but does NOT state:
#   - soft-thresholding power
#   - unsigned / signed network
#
# We scan fractional powers 0.1 to 1.0 (plus 1.0 exactly).
# Integer powers above 1 produce very few strong pairs on
# this dataset and are not informative for the comparison.
# ------------------------------------------------------------

powers <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

network_types <- c("unsigned", "signed")

qc_rows <- list()
row_id  <- 1

for (network_type in network_types) {
  for (power in powers) {
    
    cat("\nRunning:", network_type, "power =", power, "\n")
    
    # Favorable
    adj_fav <- adjacency(datExpr_fav, type = network_type, power = power)
    fav_qc  <- count_strong_pairs(adj_fav, WEIGHT_THRESHOLD)
    rm(adj_fav); gc()
    
    # Poor
    adj_poor <- adjacency(datExpr_poor, type = network_type, power = power)
    poor_qc  <- count_strong_pairs(adj_poor, WEIGHT_THRESHOLD)
    rm(adj_poor); gc()
    
    qc_rows[[row_id]] <- data.frame(
      network_type = network_type,
      power        = power,
      
      favorable_pairs         = fav_qc$n_pairs,
      favorable_strong_pairs  = fav_qc$n_strong,
      favorable_fraction      = fav_qc$proportion,
      
      poor_pairs              = poor_qc$n_pairs,
      poor_strong_pairs       = poor_qc$n_strong,
      poor_fraction           = poor_qc$proportion,
      
      favorable_target        = PAPER_FAVORABLE_EDGES,
      poor_target             = PAPER_POOR_EDGES,
      
      favorable_difference    = fav_qc$n_strong  - PAPER_FAVORABLE_EDGES,
      poor_difference         = poor_qc$n_strong - PAPER_POOR_EDGES,
      
      stringsAsFactors = FALSE
    )
    
    row_id <- row_id + 1
  }
}

wgcna_sensitivity <- do.call(rbind, qc_rows)

write.csv(
  wgcna_sensitivity,
  file.path(OUT, "wgcna_adjacency_weight_qc.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Find closest configurations
# ------------------------------------------------------------

wgcna_sensitivity$fav_abs_error  <- abs(wgcna_sensitivity$favorable_difference)
wgcna_sensitivity$poor_abs_error <- abs(wgcna_sensitivity$poor_difference)
wgcna_sensitivity$combined_abs_error <-
  wgcna_sensitivity$fav_abs_error + wgcna_sensitivity$poor_abs_error

closest_combined <- wgcna_sensitivity[
  order(wgcna_sensitivity$combined_abs_error),
]

write.csv(
  closest_combined,
  file.path(OUT, "wgcna_closest_configurations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Print key results
# ------------------------------------------------------------

cat("\n\n============================================================\n")
cat("WGCNA NETWORK QC SUMMARY\n")
cat("============================================================\n\n")

cat("Paper target — favorable:", PAPER_FAVORABLE_EDGES, "\n")
cat("Paper target — poor:",      PAPER_POOR_EDGES, "\n\n")

cat("Closest configurations:\n\n")
print(head(closest_combined, 15))

# ------------------------------------------------------------
# 12. Explicit exact-match check
# ------------------------------------------------------------

exact_match <- subset(
  wgcna_sensitivity,
  favorable_strong_pairs == PAPER_FAVORABLE_EDGES &
    poor_strong_pairs      == PAPER_POOR_EDGES
)

cat("\n================ EXACT MATCH CHECK ================\n")

if (nrow(exact_match) == 0) {
  cat(
    "NO WGCNA adjacency configuration among\n",
    "network_type = unsigned/signed and power = 0.1:1.0\n",
    "reproduced both 4,415 and 4,050 exactly.\n\n"
  )
} else {
  cat("Exact configuration(s) found:\n")
  print(exact_match)
}

# ------------------------------------------------------------
# 13. Save compact paper QC report
# ------------------------------------------------------------

report <- data.frame(
  metric = c(
    "patients_total",
    "favorable_patients",
    "poor_patients",
    "paper_favorable_strong_pairs",
    "paper_poor_strong_pairs",
    "raw_cor_favorable_strong_pairs",
    "raw_cor_poor_strong_pairs",
    "exact_wgcna_match_found"
  ),
  value = c(
    nrow(expr),
    length(fav_idx),
    length(poor_idx),
    PAPER_FAVORABLE_EDGES,
    PAPER_POOR_EDGES,
    raw_fav$n_strong,
    raw_poor$n_strong,
    nrow(exact_match) > 0
  )
)

write.csv(
  report,
  file.path(OUT, "wgcna_network_qc_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Save session information
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(OUT, "sessionInfo.txt")
)

cat("\nQC files written to:\n")
cat(OUT, "\n")

cat("\nDONE.\n")