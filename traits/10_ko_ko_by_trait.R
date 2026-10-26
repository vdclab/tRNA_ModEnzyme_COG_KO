#!/usr/bin/env Rscript
# Step 10: Pairwise KO-KO within-pathway correlations, computed separately
# within each selected trait's complete-case subset.
# For each trait: subset genomes to complete cases, prune tree, compute all
# unordered member-KO pairs per pathway (phi, Spearman, odds ratio + Fisher
# exact p, phylogenetic logistic p via phylolm), BH-FDR across that trait's
# pairs. Output: <trait>/ko_ko_within_pathway_<trait>.csv
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2"

ko     <- read.csv(file.path(PROC, "ko_binary.csv"), row.names = 1)
traits <- read.csv(file.path(PROC, "traits_selected.csv"), row.names = 1,
                   na.strings = c("", "NA"))
mem    <- read.csv(file.path(PROC, "pathway_members.csv"), stringsAsFactors = FALSE)
tree   <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")
mem$KO_list <- regmatches(mem$KOs, gregexpr("K[0-9]{5}", mem$KOs))

trait_info <- data.frame(
  trait = c("Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"),
  type  = c("binary", "binary", "Categorical", "Categorical"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(trait_info))) {
  trait <- trait_info$trait[i]
  keep <- rownames(ko)[!is.na(traits[rownames(ko), trait])]
  ko_sub <- ko[keep, ]
  tr_full <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% keep])
  cat("=== ", trait, ": ", length(keep), "genomes ===\n", sep = "")

  rows <- list()
  for (j in seq_len(nrow(mem))) {
    pathway <- mem$pathway[j]
    kos <- unlist(mem$KO_list[j]); kos <- kos[kos %in% colnames(ko_sub)]
    if (length(kos) < 2) next
    for (pr in combn(kos, 2, simplify = FALSE)) {
      k1 <- pr[1]; k2 <- pr[2]
      d <- data.frame(x = ko_sub[[k1]], y = ko_sub[[k2]], row.names = rownames(ko_sub))
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
      tr <- drop.tip(tr_full, tr_full$tip.label[!tr_full$tip.label %in% rownames(d)])
      fit <- tryCatch(phyloglm(y ~ x, data = d, phy = tr, method = "logistic_MPLE",
                               boot = 0, btol = 50), error = function(e) NULL)
      if (!is.null(fit) && fit$convergence == 0) {
        method <- "phyloglm_logistic_MPLE"
        est <- unname(coef(fit)["x"]); p_phy <- unname(summary(fit)$coefficients["x", "p.value"])
      } else if (!is.null(fit)) {
        method <- "phyloglm_not_converged"
        est <- unname(coef(fit)["x"]); p_phy <- unname(summary(fit)$coefficients["x", "p.value"])
      } else {
        method <- "failed"; est <- NA; p_phy <- NA
      }
      rows[[length(rows) + 1]] <- data.frame(trait, pathway, KO_1 = k1, KO_2 = k2, n,
        method, phi, spearman_rho = unname(sp$estimate), odds_ratio = or,
        fisher_p = p_fisher, estimate = est, p_phylo = p_phy)
    }
  }
  res <- do.call(rbind, rows)
  res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
  res$q_fisher <- p.adjust(res$fisher_p, method = "BH")
  write.csv(res, file.path(RES, trait,
                           paste0("ko_ko_within_pathway_", trait, ".csv")),
            row.names = FALSE)
  cat("  pairs:", nrow(res), "; failures:", sum(res$method == "failed"),
      "; sig q<0.05 (phylo):", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
}
cat("Step 10 done.\n")
