#!/usr/bin/env Rscript
# Step 2: Trait x pathway associations with phylogenetic models.
# Binary traits  -> phylolm::phyloglm (fallback ape::binaryPGLMM, then naive glm)
# Categorical    -> geiger::phylo.anova (fallback naive Kruskal-Wallis flag)
# Naive stats computed alongside; BH-FDR within each trait family.

.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(ape); library(phylolm); library(phytools)
})

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

comp   <- read.csv(file.path(PROC, "pathway_completeness.csv"), row.names = 1)
traits <- read.csv(file.path(PROC, "traits_selected.csv"), row.names = 1,
                   na.strings = c("", "NA"))
tree   <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")

trait_info <- data.frame(
  trait = c("Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"),
  type  = c("binary", "binary", "Categorical", "Categorical"),
  stringsAsFactors = FALSE
)

run_one_binary <- function(trait, pathway) {
  ycol <- if (paste0(trait, "_coded") %in% colnames(traits)) paste0(trait, "_coded") else trait
  d <- data.frame(
    y = as.numeric(traits[[ycol]]),
    x = comp[[pathway]],
    row.names = rownames(comp)
  )
  d <- d[complete.cases(d), ]
  n <- nrow(d)
  if (length(unique(d$y)) < 2 || length(unique(d$x)) < 2)
    return(data.frame(trait, pathway, n, method = "skipped_degenerate",
                      estimate = NA, odds_ratio = NA, p_phylo = NA,
                      pearson_r = NA, pearson_p = NA,
                      spearman_rho = NA, spearman_p = NA))
  tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(d)])

  # naive stats
  pearson  <- suppressWarnings(cor.test(d$x, d$y, method = "pearson"))
  spearman <- suppressWarnings(cor.test(d$x, d$y, method = "spearman",
                                        exact = FALSE))

  # phylogenetic logistic regression
  fit <- tryCatch(
    phyloglm(y ~ x, data = d, phy = tr, method = "logistic_MPLE", boot = 0,
             btol = 50),
    error = function(e) NULL
  )
  if (!is.null(fit) && fit$convergence == 0) {
    method <- "phyloglm_logistic_MPLE"
    est <- coef(fit)["x"]; p <- summary(fit)$coefficients["x", "p.value"]
  } else {
    fit2 <- tryCatch(binaryPGLMM(y ~ x, data = d, phy = tr), error = function(e) NULL)
    if (!is.null(fit2)) {
      method <- "binaryPGLMM_fallback"
      est <- fit2$B["x", 1]; p <- fit2$p.value[2]
    } else {
      fit3 <- glm(y ~ x, data = d, family = binomial)
      method <- "naive_glm_fallback"
      s <- summary(fit3)$coefficients
      est <- s["x", "Estimate"]; p <- s["x", "Pr(>|z|)"]
    }
  }
  data.frame(trait, pathway, n, method, estimate = est,
             odds_ratio = exp(est), p_phylo = p,
             pearson_r = unname(pearson$estimate),
             pearson_p = pearson$p.value,
             spearman_rho = unname(spearman$estimate),
             spearman_p = spearman$p.value)
}

run_one_categorical <- function(trait, pathway) {
  d <- data.frame(
    g = factor(traits[[trait]]),
    x = comp[[pathway]],
    row.names = rownames(comp)
  )
  d <- d[complete.cases(d), ]
  n <- nrow(d)
  ng <- nlevels(droplevels(d$g))
  if (ng < 2 || length(unique(d$x)) < 2)
    return(data.frame(trait, pathway, n, n_groups = ng, method = "skipped_degenerate",
                      F_phylo = NA, p_phylo = NA,
                      kruskal_chi2 = NA, kruskal_p = NA, epsilon_sq = NA))
  tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(d)])
  d$g <- droplevels(d$g)

  kw <- kruskal.test(x ~ g, data = d)
  # epsilon-squared effect size
  rank_x <- rank(d$x); gmean <- mean(rank_x)
  eps <- sum(tapply(rank_x, d$g, function(v) length(v) * (mean(v) - gmean)^2)) /
    ((n - 1) * var(rank_x))

  # Simulation-based phylogenetic ANOVA (Adams 2012 style):
  # null distribution of ANOVA F from Brownian simulations on the tree,
  # with BM rate estimated from phylogenetically independent contrasts.
  # pic() requires a rooted tree -> midpoint-root a copy for the simulation
  # step only (no outgroup available; models themselves use the unrooted tree).
  F_obs <- unname(anova(lm(x ~ g, data = d))$`F value`[1])
  tr_r <- tr  # tree pre-rooted at Aquificota
  pic_x <- pic(d$x, tr_r, scaled = TRUE)
  sigma2 <- sum(pic_x^2) / (n - 1)
  sigma2 <- max(sigma2, .Machine$double.eps)
  set.seed(42)
  F_sim <- replicate(1000, {
    xs <- fastBM(tr_r, sig2 = sqrt(sigma2))[rownames(d)]
    dd <- data.frame(xs = as.numeric(xs), g = d$g)
    anova(lm(xs ~ g, data = dd))$`F value`[1]
  })
  p_sim <- mean(F_sim >= F_obs)
  data.frame(trait, pathway, n, n_groups = ng, method = "phylo_anova_BM_sim",
             F_phylo = F_obs, p_phylo = p_sim,
             kruskal_chi2 = unname(kw$statistic), kruskal_p = kw$p.value,
             epsilon_sq = eps)
}

for (i in seq_len(nrow(trait_info))) {
  trait <- trait_info$trait[i]
  type  <- trait_info$type[i]
  cat("=== trait:", trait, "(", type, ") ===\n")
  res <- do.call(rbind, lapply(colnames(comp), function(p) {
    if (type == "binary") run_one_binary(trait, p)
    else run_one_categorical(trait, p)
  }))
  if (type == "binary") {
    res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
    res$q_spearman <- p.adjust(res$spearman_p, method = "BH")
  } else {
    res$q_phylo <- p.adjust(res$p_phylo, method = "BH")
    res$q_kruskal <- p.adjust(res$kruskal_p, method = "BH")
  }
  outdir <- file.path(RES, trait)
  dir.create(outdir, showWarnings = FALSE)
  write.csv(res, file.path(outdir, "trait_pathway_stats.csv"), row.names = FALSE)
  cat("  fits:", sum(!is.na(res$p_phylo)), "/", nrow(res),
      "; significant q<0.05:", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
}
cat("Step 2 done.\n")
