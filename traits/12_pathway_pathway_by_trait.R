#!/usr/bin/env Rscript
# Step 12: Pathway x pathway correlations computed separately within each
# selected trait's complete-case subset.
# For each trait: subset genomes to complete cases, prune tree, compute all
# 528 pathway-completeness pairs (Pearson, Spearman, phylolm lambda-model
# regression -> phylogeny-corrected r and p), BH-FDR within trait.
# Output: <trait>/pathway_pathway_<Trait>.csv (long format, one row per pair)
#         <trait>/pathway_pathway_r_matrix_<Trait>.csv (33x33, for plotting)
#         <trait>/pathway_pathway_q_matrix_<Trait>.csv
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2"

comp   <- read.csv(file.path(PROC, "pathway_completeness.csv"), row.names = 1)
traits <- read.csv(file.path(PROC, "traits_selected.csv"), row.names = 1,
                   na.strings = c("", "NA"))
tree   <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")
paths  <- colnames(comp)

trait_info <- data.frame(
  trait = c("Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"),
  stringsAsFactors = FALSE
)

for (trait in trait_info$trait) {
  keep <- rownames(comp)[!is.na(traits[rownames(comp), trait])]
  comp_sub <- comp[keep, ]
  tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep])
  cat("=== ", trait, ": ", length(keep), "genomes ===\n", sep = "")

  rows <- list(); idx <- 0
  for (i in seq_along(paths)) {
    for (j in (i + 1):length(paths)) {
      if (j > length(paths)) break
      p1 <- paths[i]; p2 <- paths[j]
      x <- comp_sub[[p1]]; y <- comp_sub[[p2]]
      pe <- suppressWarnings(cor.test(x, y, method = "pearson"))
      sp <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      d <- data.frame(y = y, x = x, row.names = rownames(comp_sub))
      fit <- tryCatch(phylolm(y ~ x, data = d, phy = tr, model = "lambda"),
                      error = function(e) NULL)
      if (is.null(fit)) {
        slope <- NA; p_pgls <- NA; r_pgls <- NA; lambda <- NA
      } else {
        s <- summary(fit)
        slope <- s$coefficients["x", "Estimate"]
        p_pgls <- s$coefficients["x", "p.value"]
        tstat <- s$coefficients["x", "t.value"]
        r_pgls <- sign(tstat) * sqrt(tstat^2 / (tstat^2 + fit$n - 2))
        lambda <- s$optpar
      }
      idx <- idx + 1
      rows[[idx]] <- data.frame(trait, pathway_1 = p1, pathway_2 = p2, n = length(x),
        pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
        spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
        pgls_slope = slope, pgls_lambda = lambda, r_phylo = r_pgls, p_phylo = p_pgls)
    }
  }
  res <- do.call(rbind, rows)
  res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
  res$q_pearson <- p.adjust(res$pearson_p, method = "BH")
  write.csv(res, file.path(RES, trait, paste0("pathway_pathway_", trait, ".csv")),
            row.names = FALSE)

  n <- length(paths)
  R <- diag(n); dimnames(R) <- list(paths, paths)
  Q <- matrix(NA_real_, n, n, dimnames = list(paths, paths))
  for (k in seq_len(nrow(res))) {
    a <- res$pathway_1[k]; b <- res$pathway_2[k]
    R[a, b] <- R[b, a] <- res$r_phylo[k]
    Q[a, b] <- Q[b, a] <- res$q_phylo[k]
  }
  diag(Q) <- 0
  write.csv(R, file.path(RES, trait, paste0("pathway_pathway_r_matrix_", trait, ".csv")))
  write.csv(Q, file.path(RES, trait, paste0("pathway_pathway_q_matrix_", trait, ".csv")))
  cat("  pairs:", nrow(res), "; failures:", sum(is.na(res$p_phylo)),
      "; sig q<0.05 (phylo):", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
}
cat("Step 12 done.\n")
