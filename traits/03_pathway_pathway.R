#!/usr/bin/env Rscript
# Step 3: Pathway x pathway correlations across 1322 genomes.
# For each of the 528 unordered pairs of pathway completeness vectors:
#   - Pearson r (naive)
#   - Spearman rho (naive)
#   - Phylogenetic regression (phylolm, lambda model) -> slope p-value,
#     converted to a phylogeny-corrected correlation r via t-statistic.
# BH-FDR across all pairs.
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2/01_pathway_pathway"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

comp <- read.csv(file.path(PROC, "pathway_completeness.csv"), row.names = 1)
tree <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")
paths <- colnames(comp)

pair_rows <- list()
idx <- 0
for (i in seq_along(paths)) {
  for (j in (i + 1):length(paths)) {
    if (j > length(paths)) break
    p1 <- paths[i]; p2 <- paths[j]
    x <- comp[[p1]]; y <- comp[[p2]]
    pe <- suppressWarnings(cor.test(x, y, method = "pearson"))
    sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    d <- data.frame(y = y, x = x, row.names = rownames(comp))
    fit <- tryCatch(phylolm(y ~ x, data = d, phy = tree, model = "lambda"),
                    error = function(e) NULL)
    if (is.null(fit)) {
      slope <- NA; p_pgls <- NA; r_pgls <- NA; lambda <- NA
    } else {
      s <- summary(fit)
      slope <- s$coefficients["x", "Estimate"]
      p_pgls <- s$coefficients["x", "p.value"]
      tstat <- s$coefficients["x", "t.value"]
      dfres <- fit$n - 2
      r_pgls <- sign(tstat) * sqrt(tstat^2 / (tstat^2 + dfres))
      lambda <- s$optpar
    }
    idx <- idx + 1
    pair_rows[[idx]] <- data.frame(pathway_1 = p1, pathway_2 = p2,
      pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
      spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
      pgls_slope = slope, pgls_lambda = lambda, r_phylo = r_pgls, p_phylo = p_pgls)
  }
  if (i %% 10 == 0) cat("row", i, "of", length(paths), "- pairs done:", idx, "\n")
}
res <- do.call(rbind, pair_rows)
res$q_pearson <- p.adjust(res$pearson_p, method = "BH")
res$q_phylo   <- p.adjust(res$p_phylo, method = "BH")
write.csv(res, file.path(RES, "pathway_pathway_stats.csv"), row.names = FALSE)

cat("pairs tested:", nrow(res), "\n")
cat("phylolm fit failures:", sum(is.na(res$p_phylo)), "\n")
cat("significant q<0.05 (phylo):", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
cat("significant q<0.05 (naive pearson):", sum(res$q_pearson < 0.05, na.rm = TRUE), "\n")

# 33x33 symmetric matrices for plotting
n <- length(paths)
R <- diag(n); dimnames(R) <- list(paths, paths)
Rq <- matrix(NA_real_, n, n, dimnames = list(paths, paths))
for (k in seq_len(nrow(res))) {
  a <- res$pathway_1[k]; b <- res$pathway_2[k]
  R[a, b] <- R[b, a] <- res$r_phylo[k]
  Rq[a, b] <- Rq[b, a] <- res$q_phylo[k]
}
diag(Rq) <- 0
write.csv(R,  file.path(RES, "r_phylo_matrix.csv"))
write.csv(Rq, file.path(RES, "q_phylo_matrix.csv"))
cat("Step 3 done.\n")
