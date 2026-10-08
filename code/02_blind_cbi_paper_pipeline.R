# ============================================================
# CBI PIPELINE — REPRODUCES THE 30-GENE CANDIDATE SET
# ============================================================
# Uses:
#   49 benchmark genes
#   All 146 patients pooled
#   bnlearn bootstrap (R = 100, score = bic-g, CPDAG = TRUE)
#   averaged.network(threshold = 0.5) — CONSENSUS NETWORK
#   total_degree >= 3
# ============================================================

library(bnlearn)

set.seed(20261005)

ROOT <- "C:/Users/MARAO/Downloads/benchmark/benchmark"

# --- load ---
expr <- read.csv(file.path(ROOT, "blinded", "synthetic_expression.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
clin <- read.csv(file.path(ROOT, "blinded", "synthetic_clinical.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
module_map <- read.csv(file.path(ROOT, "metadata", "predefined_module_map.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)

# --- 49 genes ---
all_genes <- intersect(module_map$gene, colnames(expr))
stopifnot(length(all_genes) == 49)

# --- pooled data (all 146 patients) ---
dat <- expr[, all_genes, drop = FALSE]
dat <- as.data.frame(lapply(dat, as.numeric), check.names = FALSE)

cat("Pooled data:", nrow(dat), "patients x", ncol(dat), "genes\n\n")

# --- bootstrap ---
set.seed(20261005)

cat("Running bootstrap structure learning (R = 100)...\n")
bs <- boot.strength(
  data           = dat,
  R              = 100,
  algorithm      = "hc",
  algorithm.args = list(score = "bic-g"),
  cpdag          = TRUE
)

cat("Bootstrap arcs evaluated:", nrow(bs), "\n\n")

# --- averaged network at threshold 0.5 ---
cat("Building averaged network at threshold = 0.5...\n")
avg_net <- averaged.network(bs, threshold = 0.5)

arc_mat <- arcs(avg_net)
arcs_df <- data.frame(
  from = arc_mat[, 1],
  to   = arc_mat[, 2],
  stringsAsFactors = FALSE
)

cat("Averaged network edges:", nrow(arcs_df), "\n\n")

# --- degree ---
nodes <- sort(unique(c(arcs_df$from, arcs_df$to)))

in_deg  <- table(factor(arcs_df$to,   levels = nodes))
out_deg <- table(factor(arcs_df$from, levels = nodes))

deg_df <- data.frame(
  gene         = nodes,
  in_degree    = as.integer(in_deg),
  out_degree   = as.integer(out_deg),
  total_degree = as.integer(in_deg) + as.integer(out_deg),
  stringsAsFactors = FALSE
)

deg_df <- deg_df[order(-deg_df$total_degree), ]

cat("Top 30 genes by degree:\n")
print(head(deg_df, 30))
cat("\n")

# --- candidates ---
candidates <- deg_df$gene[deg_df$total_degree >= 3]

cat("Genes with total degree >= 3:", length(candidates), "\n")
print(candidates)

# --- save ---
write.csv(data.frame(gene = candidates),
          file.path(ROOT, "results", "pooled_49_module_audit",
                    "module_union_candidates_cpdag_TRUE.csv"),
          row.names = FALSE)

cat("\nSaved candidate genes.\n")