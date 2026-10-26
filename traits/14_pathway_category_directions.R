#!/usr/bin/env Rscript
# Step 14: Category-level direction of pathway-trait associations.
# For each selected trait and each pathway, test each trait category against
# the rest: phylogenetic logistic regression (phylolm::phyloglm,
# in_category ~ completeness) -> log OR per unit completeness, plus
# descriptive completeness stats per category. BH-FDR within trait across
# all pathway x category tests. Omnibus trait-level p (from step 2) merged
# for context.
# Output: <trait>/pathway_category_associations_<Trait>.csv
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2"

comp   <- read.csv(file.path(PROC, "pathway_completeness.csv"), row.names = 1)
traits <- read.csv(file.path(PROC, "traits_selected.csv"), row.names = 1,
                   na.strings = c("", "NA"))
tree   <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")

trait_info <- data.frame(
  trait  = c("Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"),
  type   = c("binary", "binary", "Categorical", "Categorical"),
  stringsAsFactors = FALSE
)

fit_cat_vs_rest <- function(d, tr) {
  # d: y = 1 if in category else 0; x = completeness
  fit <- tryCatch(phyloglm(y ~ x, data = d, phy = tr, method = "logistic_MPLE",
                           boot = 0, btol = 50), error = function(e) NULL)
  if (!is.null(fit) && fit$convergence == 0) {
    list(method = "phyloglm_logistic_MPLE", est = unname(coef(fit)["x"]),
         p = unname(summary(fit)$coefficients["x", "p.value"]))
  } else if (!is.null(fit)) {
    list(method = "phyloglm_not_converged", est = unname(coef(fit)["x"]),
         p = unname(summary(fit)$coefficients["x", "p.value"]))
  } else {
    fit2 <- tryCatch({
      setTimeLimit(elapsed = 60, transient = TRUE); on.exit(setTimeLimit(), add = TRUE)
      binaryPGLMM(y ~ x, data = d, phy = tr)
    }, error = function(e) NULL)
    if (!is.null(fit2)) list(method = "binaryPGLMM_fallback",
                             est = unname(fit2$B["x", 1]), p = unname(fit2$p.value[2]))
    else {
      s <- summary(glm(y ~ x, data = d, family = binomial))$coefficients
      list(method = "naive_glm_fallback", est = unname(s["x", "Estimate"]),
           p = unname(s["x", "Pr(>|z|)"]))
    }
  }
}

for (i in seq_len(nrow(trait_info))) {
  trait <- trait_info$trait[i]; type <- trait_info$type[i]
  # omnibus results from step 2 for context
  omni <- read.csv(file.path(RES, trait, "trait_pathway_stats.csv"))
  omni_map <- setNames(omni$p_phylo, omni$pathway)

  # category labels per genome (complete cases only)
  raw <- traits[[trait]]
  names(raw) <- rownames(traits)
  raw <- raw[!is.na(raw)]
  cats <- sort(unique(raw))

  rows <- list()
  for (pathway in colnames(comp)) {
    x_all <- comp[names(raw), pathway]
    for (cat in cats) {
      y <- as.numeric(raw == cat)
      d <- data.frame(y = y, x = x_all, row.names = names(raw))
      d <- d[complete.cases(d), ]
      n <- nrow(d); n_cat <- sum(d$y)
      in_cat <- d$y == 1
      mean_all <- mean(d$x); mean_cat <- mean(d$x[in_cat])
      mean_rest <- mean(d$x[!in_cat])
      if (n_cat < 10 || n - n_cat < 10 || length(unique(d$x)) < 2) {
        rows[[length(rows) + 1]] <- data.frame(trait, pathway, category = cat,
          n = n, n_category = n_cat, mean_completeness_in_category = mean_cat,
          mean_completeness_rest = mean_rest, diff_vs_rest = mean_cat - mean_rest,
          method = "skipped_small_category", estimate_logOR = NA,
          odds_ratio = NA, p_phylo = NA, omnibus_p_phylo = unname(omni_map[pathway]))
        next
      }
      tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(d)])
      out <- fit_cat_vs_rest(d, tr)
      rows[[length(rows) + 1]] <- data.frame(trait, pathway, category = cat,
        n = n, n_category = n_cat, mean_completeness_in_category = mean_cat,
        mean_completeness_rest = mean_rest, diff_vs_rest = mean_cat - mean_rest,
        method = out$method, estimate_logOR = out$est, odds_ratio = exp(out$est),
        p_phylo = out$p, omnibus_p_phylo = unname(omni_map[pathway]))
    }
  }
  res <- do.call(rbind, rows)
  res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
  outdir <- file.path(RES, trait)
  write.csv(res, file.path(outdir, paste0("pathway_category_associations_",
                                          trait, ".csv")), row.names = FALSE)
  cat("=== ", trait, ": ", nrow(res), " pathway x category tests; sig q<0.05: ",
      sum(res$q_phylo < 0.05, na.rm = TRUE), "\n", sep = "")
}
cat("Step 14 done.\n")
