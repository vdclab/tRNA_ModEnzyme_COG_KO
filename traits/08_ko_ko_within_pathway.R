#!/usr/bin/env Rscript
# Step 8: Pairwise correlations between KOs within the same pathway.
# For every unordered pair of member KOs in each pathway with >= 2 KOs:
#   - phi coefficient (binary correlation)
#   - Spearman rho
#   - odds ratio (2x2 table, 0.5 continuity correction) + Fisher exact p
#   - phylogenetic logistic regression (phylolm, KO_2 ~ KO_1) -> p
# BH-FDR across all pairs. Output: 02_within_pathway/ko_ko_within_pathway_stats.csv
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2/02_within_pathway"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

ko   <- read.csv(file.path(PROC, "ko_binary.csv"), row.names = 1)
mem  <- read.csv(file.path(PROC, "pathway_members.csv"), stringsAsFactors = FALSE)
tree <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")

# pathway_members.csv: pathway column + KOs column (python-style list string);
# extract KO accessions (K + 5 digits) robustly
mem$KO_list <- regmatches(mem$KOs, gregexpr("K[0-9]{5}", mem$KOs))

rows <- list()
for (i in seq_len(nrow(mem))) {
  pathway <- mem$pathway[i]
  kos <- trimws(unlist(mem$KO_list[i]))
  kos <- kos[kos %in% colnames(ko)]
  if (length(kos) < 2) next
  prs <- combn(kos, 2, simplify = FALSE)
  for (pr in prs) {
    k1 <- pr[1]; k2 <- pr[2]
    d <- data.frame(x = ko[[k1]], y = ko[[k2]], row.names = rownames(ko))
    d <- d[complete.cases(d), ]
    n <- nrow(d)
    tab <- table(d$x, d$y)
    phi <- NA; or <- NA; p_fisher <- NA
    if (all(dim(tab) == 2)) {
      chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
      phi <- unname(sqrt(chi / n) * sign(tab[2, 2] * tab[1, 1] -
                                           tab[2, 1] * tab[1, 2]))
      ft <- fisher.test(tab)
      p_fisher <- ft$p.value
      or <- unname(ft$estimate)
    }
    sp <- suppressWarnings(cor.test(d$x, d$y, method = "spearman", exact = FALSE))
    # phylogenetic logistic regression (direction arbitrary for a pair)
    tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(d)])
    fit <- tryCatch(phyloglm(y ~ x, data = d, phy = tr, method = "logistic_MPLE",
                             boot = 0, btol = 50), error = function(e) NULL)
    if (!is.null(fit) && fit$convergence == 0) {
      method <- "phyloglm_logistic_MPLE"
      est <- unname(coef(fit)["x"]); p_phy <- unname(summary(fit)$coefficients["x", "p.value"])
    } else if (!is.null(fit)) {
      method <- "phyloglm_not_converged"
      est <- unname(coef(fit)["x"]); p_phy <- unname(summary(fit)$coefficients["x", "p.value"])
    } else {
      method <- "failed"
      est <- NA; p_phy <- NA
    }
    rows[[length(rows) + 1]] <- data.frame(pathway, KO_1 = k1, KO_2 = k2, n,
      method, phi, spearman_rho = unname(sp$estimate), odds_ratio = or,
      fisher_p = p_fisher, estimate = est, p_phylo = p_phy)
  }
}
res <- do.call(rbind, rows)
res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
res$q_fisher <- p.adjust(res$fisher_p, method = "BH")
write.csv(res, file.path(RES, "ko_ko_within_pathway_stats.csv"), row.names = FALSE)
cat("pairs tested:", nrow(res), "\n")
cat("phyloglm failures:", sum(res$method == "failed"), "\n")
cat("significant q<0.05 (phylo):", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
cat("significant q<0.05 (fisher):", sum(res$q_fisher < 0.05, na.rm = TRUE), "\n")
cat("Step 8 done.\n")
