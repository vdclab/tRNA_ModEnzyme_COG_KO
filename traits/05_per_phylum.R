#!/usr/bin/env Rscript
# Step 5: Per-phylum sensitivity analysis.
# For each globally significant trait-pathway pair (q_phylo < 0.05 in step 2),
# rerun the same model within each phylum with >= 10 genomes AND >= 10
# complete cases for that trait. Binary traits: phyloglm. Categorical traits:
# BM-simulation phylogenetic ANOVA (500 sims).
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({ library(ape); library(phylolm); library(phytools) })

PROC <- "/mnt/results/ko_trait_analysis_v2/00_processed"
RES  <- "/mnt/results/ko_trait_analysis_v2"
OUT  <- file.path(RES, "03_per_phylum")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

comp   <- read.csv(file.path(PROC, "pathway_completeness.csv"), row.names = 1)
traits <- read.csv(file.path(PROC, "traits_selected.csv"), row.names = 1,
                   na.strings = c("", "NA"))
tree   <- read.tree("/workspace/ko_trait_analysis_v2/bac120_aquificota_rooted.tree")

trait_info <- data.frame(
  trait = c("Gram_Stain", "Motility", "Oxygen_Req.", "Growth_Temp_Category"),
  type  = c("binary", "binary", "Categorical", "Categorical"),
  stringsAsFactors = FALSE
)

run_binary <- function(d, tr) {
  fit <- tryCatch(phyloglm(y ~ x, data = d, phy = tr, method = "logistic_MPLE",
                           boot = 0, btol = 50), error = function(e) NULL)
  if (!is.null(fit) && fit$convergence == 0) {
    c(method = "phyloglm", est = unname(coef(fit)["x"]),
      p = unname(summary(fit)$coefficients["x", "p.value"]))
  } else {
    fit2 <- tryCatch(binaryPGLMM(y ~ x, data = d, phy = tr), error = function(e) NULL)
    if (!is.null(fit2)) c(method = "binaryPGLMM", est = unname(fit2$B["x", 1]),
                          p = unname(fit2$p.value[2]))
    else c(method = "failed", est = NA, p = NA)
  }
}

run_categorical <- function(d, tr) {
  F_obs <- unname(anova(lm(x ~ g, data = d))$`F value`[1])
  tr_r <- tr
  pic_x <- pic(d$x, tr_r, scaled = TRUE)
  sigma2 <- max(sum(pic_x^2) / (nrow(d) - 1), .Machine$double.eps)
  set.seed(42)
  F_sim <- replicate(500, {
    xs <- fastBM(tr_r, sig2 = sqrt(sigma2))[rownames(d)]
    dd <- data.frame(xs = as.numeric(xs), g = d$g)
    anova(lm(xs ~ g, data = dd))$`F value`[1]
  })
  c(method = "phylo_anova_BM_sim", est = unname(F_obs),
    p = unname(mean(F_sim >= F_obs)))
}

all_rows <- list()
for (i in seq_len(nrow(trait_info))) {
  trait <- trait_info$trait[i]; type <- trait_info$type[i]
  stats <- read.csv(file.path(RES, trait, "trait_pathway_stats.csv"))
  sig <- stats[!is.na(stats$q_phylo) & stats$q_phylo < 0.05, ]
  if (nrow(sig) == 0) next
  ycol <- if (type == "binary" && paste0(trait, "_coded") %in% colnames(traits))
    paste0(trait, "_coded") else trait

  for (ph in unique(traits$Phylum)) {
    genomes <- rownames(traits)[!is.na(traits$Phylum) & traits$Phylum == ph]
    if (length(genomes) < 10) next
    for (k in seq_len(nrow(sig))) {
      pathway <- sig$pathway[k]
      yraw <- traits[genomes, ycol]
      yv <- if (type == "binary") as.numeric(yraw) else as.numeric(factor(yraw))
      d <- data.frame(y = yv,
                      x = comp[genomes, pathway],
                      g = factor(traits[genomes, trait]),
                      row.names = genomes)
      d <- d[complete.cases(d), c("y", "x", "g")]
      if (nrow(d) < 10) next
      if (type == "Categorical") d$g <- droplevels(d$g)
      if (length(unique(d$y)) < 2 || length(unique(d$x)) < 2 ||
          (type == "Categorical" && nlevels(d$g) < 2)) next
      tr <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rownames(d)])
      out <- if (type == "binary") run_binary(d, tr) else run_categorical(d, tr)
      all_rows[[length(all_rows) + 1]] <- data.frame(
        trait, pathway, phylum = ph, n = nrow(d), method = out["method"],
        estimate = unname(as.numeric(out["est"])),
        p_phylo = unname(as.numeric(out["p"])))
    }
  }
}
res <- do.call(rbind, all_rows)
if (nrow(res)) {
  for (tr_ in unique(res$trait)) {  # FDR within trait family
    m <- res$trait == tr_
    res$q_phylo[m] <- p.adjust(res$p_phylo[m], method = "BH")
  }
  write.csv(res, file.path(OUT, "per_phylum_stats.csv"), row.names = FALSE)
  cat("per-phylum tests:", nrow(res), "\n")
  cat("significant q<0.05:", sum(res$q_phylo < 0.05, na.rm = TRUE), "\n")
} else cat("no significant global pairs to rerun\n")
cat("Step 5 done.\n")
