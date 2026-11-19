## ============================================================================
## 02 Binary trait ~ modification correlation -- trait as RESPONSE, mod as predictor
## Rooting: midpoint (per pruned trait tree)
## Primary test: phyloglm (logistic MPLE, btol=10); effect size: phylolm (lambda) -> r_phylo
## Outputs: /mnt/results/trait_mod_correlation_v2_trait_response/
## ============================================================================
suppressPackageStartupMessages({
  library(phylolm)
  library(ape)
  library(phytools)
  library(ggplot2)
})

indir  <- "/mnt/user-uploads/input_v4"
outdir <- "/mnt/results/trait_mod_correlation_v2_trait_response"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot_msg <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

## ---------------------------------------------------------------- read inputs
mod <- read.delim(file.path(indir, "mod_AP.txt"), check.names = FALSE,
                  stringsAsFactors = FALSE)
stopifnot_msg(names(mod)[1] == "Assembly", "mod_AP.txt: first column is not 'Assembly'")
rownames(mod) <- mod$Assembly
mod$Assembly <- NULL
mod <- as.matrix(mod)

traits <- read.delim(file.path(indir, "genome_traits.txt"), check.names = FALSE,
                     stringsAsFactors = FALSE, na.strings = c("", "NA"))
rownames(traits) <- traits$Assembly
traits$Assembly <- NULL

tree <- read.tree(file.path(indir, "bac120.fasttree.unrooted.tree"))

sel <- read.delim(file.path(indir, "select_trait.txt"), header = FALSE,
                  stringsAsFactors = FALSE, col.names = c("trait", "type"))
binary_traits <- sel$trait[sel$type == "binary"]

cat("== QC ==\n")

## QC: modification matrix: exactly 36 columns, only 0/1, no blanks
stopifnot_msg(ncol(mod) == 36,
              sprintf("QC FAIL: expected 36 modification columns, found %d", ncol(mod)))
stopifnot_msg(!any(is.na(mod)), "QC FAIL: blank/missing values in mod_AP.txt")
stopifnot_msg(all(mod == 0 | mod == 1), "QC FAIL: mod_AP.txt contains values other than 0/1")
cat(sprintf("mod_AP.txt: %d genomes x %d modifications, all 0/1, no blanks -- OK\n",
            nrow(mod), ncol(mod)))

## QC: trait labels and exact expected counts
expected <- list(
  Gram_Stain = list(labels = c("negative", "positive"),
                    annotated = 649, missing = 673,
                    sizes = c(negative = 443, positive = 206)),
  Motility   = list(labels = c("nonmotile", "motile"),
                    annotated = 503, missing = 819,
                    sizes = c(nonmotile = 234, motile = 269))
)
for (tr in names(expected)) {
  stopifnot_msg(tr %in% colnames(traits), sprintf("QC FAIL: trait '%s' not in genome_traits.txt", tr))
  v <- traits[[tr]]
  lab <- unique(na.omit(v))
  stopifnot_msg(setequal(lab, expected[[tr]]$labels),
                sprintf("QC FAIL: trait '%s' has unexpected labels: %s",
                        tr, paste(setdiff(lab, expected[[tr]]$labels), collapse = ", ")))
  n_ann <- sum(!is.na(v)); n_mis <- sum(is.na(v))
  sz <- table(factor(v, levels = expected[[tr]]$labels))
  ok <- n_ann == expected[[tr]]$annotated && n_mis == expected[[tr]]$missing &&
    all(sz == expected[[tr]]$sizes)
  stopifnot_msg(ok, sprintf(
    "QC FAIL: trait '%s' counts mismatch (annotated=%d missing=%d; sizes: %s)",
    tr, n_ann, n_mis, paste(names(sz), sz, sep = "=", collapse = "; ")))
  cat(sprintf("trait '%s': annotated=%d missing=%d (%s) -- OK\n", tr, n_ann, n_mis,
              paste(names(sz), sz, sep = "=", collapse = "; ")))
}

## QC: genome IDs match across tree, traits, modification matrix; no duplicates
stopifnot_msg(!any(duplicated(rownames(mod))), "QC FAIL: duplicate genome IDs in mod_AP.txt")
stopifnot_msg(!any(duplicated(rownames(traits))), "QC FAIL: duplicate genome IDs in genome_traits.txt")
stopifnot_msg(!any(duplicated(tree$tip.label)), "QC FAIL: duplicate tips in tree")
stopifnot_msg(setequal(rownames(mod), rownames(traits)) &&
                setequal(rownames(mod), tree$tip.label),
              "QC FAIL: genome IDs differ across tree, genome_traits.txt and mod_AP.txt")
cat(sprintf("genome IDs identical across tree, traits and modification matrix (n=%d), no duplicates -- OK\n",
            length(tree$tip.label)))

## explicit response coding (label -> 0/1)
lev <- list(Gram_Stain = c(negative = 0, positive = 1),
            Motility   = c(nonmotile = 0, motile = 1))

## ---------------------------------------------------------------- helpers
fit_glm_once <- function(df, tre, start.beta = NULL, start.alpha = NULL) {
  warns <- character(0)
  obj <- withCallingHandlers(
    tryCatch(
      phyloglm(trait ~ mod, data = df, phy = tre, method = "logistic_MPLE",
               btol = 10, start.beta = start.beta, start.alpha = start.alpha),
      error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, sanitize(conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  list(fit = obj, warns = warns)
}

fit_lm_once <- function(df, tre) {
  warns <- character(0)
  obj <- withCallingHandlers(
    tryCatch(
      phylolm(trait ~ mod, data = df, phy = tre, model = "lambda",
              lower.bound = 1e-7, upper.bound = 1),
      error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, sanitize(conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  list(fit = obj, warns = warns)
}

sanitize <- function(x) gsub("[\r\n\t]+", " ", x)

## fit-flag conditions on a converged phyloglm fit (slope = modification predictor)
glm_flags <- function(f, warns) {
  se <- tryCatch(sqrt(diag(vcov(f))), error = function(e) rep(NA_real_, 2))
  cf <- coef(f)
  c(nonconvergence = !identical(as.integer(f$convergence), 0L),
    alpha_edge     = !identical(as.integer(f$alphaWarn), 0L),
    warning        = length(warns) > 0,
    big_coef       = !is.na(cf["mod"]) && abs(cf["mod"]) > 5,
    missing_coef_se = any(is.na(cf)) || any(is.na(se)))
}
flag_names <- c(nonconvergence = "non-convergence", alpha_edge = "alpha at search-range edge",
                warning = "phyloglm warning", big_coef = "|slope coefficient| > 5",
                missing_coef_se = "missing coefficient or SE")

## ---------------------------------------------------------------- per-trait analysis
flag_rows <- list()
res_list  <- list()

for (tr in binary_traits) {
  cat(sprintf("\n== trait: %s ==\n", tr))
  tdir <- file.path(outdir, tr)
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

  ## annotated genomes only; prune, then midpoint-root once
  ids <- rownames(traits)[!is.na(traits[[tr]])]
  trv <- unname(lev[[tr]][traits[ids, tr]])
  names(trv) <- ids
  stopifnot_msg(!any(is.na(trv)), sprintf("trait '%s': unexpected label after coding", tr))
  df <- data.frame(trait = trv, row.names = ids)     # rows named by genome ID

  tre <- keep.tip(tree, ids)
  tre <- midpoint.root(tre)
  treefile <- file.path(tdir, paste0(tr, "_pruned_midpoint_rooted.newick"))
  write.tree(tre, treefile)
  n <- length(ids)
  Tmean <- mean(node.depth.edgelength(tre)[1:n])
  la_lo <- -log(Tmean) - 4; la_hi <- -log(Tmean) + 4 # phyloglm log-alpha search range
  cat(sprintf("n=%d (class0=%d class1=%d); tree pruned+midpoint-rooted; T=%.4g; log-alpha range [%.3f, %.3f]\n",
              n, sum(trv == 0), sum(trv == 1), Tmean, la_lo, la_hi))

  mods <- colnames(mod)
  rows <- vector("list", length(mods))
  flag_start <- length(flag_rows) + 1L

  for (j in seq_along(mods)) {
    mname <- mods[j]
    y <- mod[ids, mname]
    df$mod <- y                                      # predictor: modification 0/1

    ## testability gates (within this trait's annotated genomes)
    n0 <- sum(trv == 0); n1 <- sum(trv == 1)
    pres <- mean(y)
    tab <- table(factor(trv, 0:1), factor(y, 0:1))
    untested_reasons <- character(0); det <- character(0)
    if (n0 < 10) { untested_reasons <- c(untested_reasons, "fewer than 10 genomes in a trait class")
      det <- c(det, sprintf("class0 n=%d", n0)) }
    if (n1 < 10) { untested_reasons <- c(untested_reasons, "fewer than 10 genomes in a trait class")
      det <- c(det, sprintf("class1 n=%d", n1)) }
    if (pres < 0.05) { untested_reasons <- c(untested_reasons, "presence < 5%")
      det <- c(det, sprintf("presence=%.3f", pres)) }
    if (pres > 0.95) { untested_reasons <- c(untested_reasons, "presence > 95%")
      det <- c(det, sprintf("presence=%.3f", pres)) }
    if (any(tab == 0)) { untested_reasons <- c(untested_reasons, "zero cell in 2x2 table")
      idx0 <- which(tab == 0)
      det <- c(det, sprintf("zero cell at trait=%s mod=%s",
                            dimnames(tab)[[1]][(idx0 - 1) %% 2 + 1],
                            dimnames(tab)[[2]][(idx0 - 1) %/% 2 + 1])) }
    tested <- length(untested_reasons) == 0

    ## ---- primary test: phyloglm logistic MPLE (response trait, predictor mod)
    or <- ci_lo <- ci_hi <- p_glm <- alpha <- phi <- NA_real_
    glm_failed <- FALSE; refitted <- FALSE; flag_persisted <- FALSE
    glm_flag_vec <- rep(FALSE, 5); names(glm_flag_vec) <- names(flag_names)
    fail_reason <- NA_character_

    if (tested) {
      r0 <- fit_glm_once(df, tre)
      if (inherits(r0$fit, "error")) {
        glm_flag_vec["warning"] <- TRUE
        err0 <- conditionMessage(r0$fit)
      } else {
        err0 <- NA_character_
        glm_flag_vec <- glm_flags(r0$fit, r0$warns)
      }

      if (any(glm_flag_vec) || !is.na(err0)) {
        refitted <- TRUE
        set.seed(123)   # only random step; immediately before the draws, per flagged fit
        starts <- lapply(1:9, function(i) {
          sb <- runif(2, -2, 2)               # coefficients drawn first
          sa <- exp(runif(1, la_lo, la_hi))   # then alpha, uniform on log scale
          list(start.beta = sb, start.alpha = sa)
        })
        cand <- list(list(fit = r0$fit, warns = r0$warns, start = "default"))
        for (s in starts) {
          rs <- fit_glm_once(df, tre, start.beta = s$start.beta, start.alpha = s$start.alpha)
          cand <- c(cand, list(list(fit = rs$fit, warns = rs$warns, start = "random")))
        }
        ok <- Filter(function(z) !inherits(z$fit, "error") &&
                       identical(as.integer(z$fit$convergence), 0L), cand)
        if (length(ok) == 0) {
          glm_failed <- TRUE
          msgs <- vapply(cand, function(z) if (inherits(z$fit, "error"))
            conditionMessage(z$fit) else sprintf("convergence code %s", z$fit$convergence), "")
          fail_reason <- paste(unique(msgs), collapse = " | ")
        } else {
          pen <- vapply(ok, function(z) as.numeric(z$fit$penlogLik), 0)
          keep <- ok[[which.max(pen)]]
          kf <- keep$fit
          glm_flag_vec <- glm_flags(kf, keep$warns)
          flag_persisted <- any(glm_flag_vec)
          cf <- coef(kf); se <- sqrt(diag(vcov(kf)))
          b <- cf["mod"]; s <- se["mod"]
          or <- exp(b); ci_lo <- exp(b - 1.96 * s); ci_hi <- exp(b + 1.96 * s)
          p_glm <- summary(kf)$coefficients["mod", "p.value"]
          alpha <- kf$alpha
        }
      } else {
        kf <- r0$fit
        cf <- coef(kf); se <- sqrt(diag(vcov(kf)))
        b <- cf["mod"]; s <- se["mod"]
        or <- exp(b); ci_lo <- exp(b - 1.96 * s); ci_hi <- exp(b + 1.96 * s)
        p_glm <- summary(kf)$coefficients["mod", "p.value"]
        alpha <- kf$alpha
      }
      if (glm_failed) {
        det <- c(det, paste0("logistic fit failed: ", fail_reason))
        glm_flag_vec["missing_coef_se"] <- TRUE
      }
    }

    ## ---- effect size: phylolm (lambda), r_phylo (response trait, predictor mod)
    r_phylo <- p_r <- lambda <- NA_real_
    lm_flags <- character(0)
    if (tested) {
      r1 <- fit_lm_once(df, tre)
      if (inherits(r1$fit, "error")) {
        lm_flags <- c(lm_flags, "phylolm error")
        det <- c(det, paste0("phylolm error: ", conditionMessage(r1$fit)))
      } else {
        if (length(r1$warns) > 0) {
          lm_flags <- c(lm_flags, "phylolm warning")
          det <- c(det, paste0("phylolm warning: ", paste(unique(r1$warns), collapse = " | ")))
        }
        lambda <- r1$fit$optpar
        if (abs(lambda) < 1e-6 || abs(lambda - 1) < 1e-6) {
          lm_flags <- c(lm_flags, "lambda within 1e-6 of 0 or 1")
          det <- c(det, sprintf("lambda=%.8f", lambda))
        }
        co <- summary(r1$fit)$coefficients
        tv <- co["mod", "t.value"]; p_r <- co["mod", "p.value"]
        r_phylo <- tv / sqrt(tv^2 + n - 2)
      }
    }

    ## ---- naive (non-phylogenetic) Pearson correlation
    naive_r <- naive_p <- NA_real_
    if (tested) {
      ct <- suppressWarnings(cor.test(trv, y))
      naive_r <- unname(ct$estimate); naive_p <- ct$p.value
    }

    rows[[j]] <- data.frame(
      modification = mname, n = n, n_class0 = n0, n_class1 = n1,
      n_present_class0 = sum(y[trv == 0]), n_present_class1 = sum(y[trv == 1]),
      OR = or, CI_low = ci_lo, CI_high = ci_hi, p_logistic = p_glm,
      alpha = alpha, phi = phi,
      r_phylo = r_phylo, p_rphylo = p_r, lambda = lambda,
      naive_r = naive_r, naive_p = naive_p,
      tested = tested, refitted = refitted, flag_persisted = flag_persisted,
      stringsAsFactors = FALSE)

    flag_rows[[length(flag_rows) + 1]] <- data.frame(
      trait = tr, modification = mname,
      tested = ifelse(tested, "yes", "no"),
      untested_reasons = if (length(untested_reasons)) paste(unique(untested_reasons), collapse = "; ") else "",
      result_flags = "", detail = paste(det, collapse = "; "),
      stringsAsFactors = FALSE)
    attr(flag_rows[[length(flag_rows)]], "glm_flags") <- if (tested) glm_flag_vec else NULL
    attr(flag_rows[[length(flag_rows)]], "lm_flags") <- lm_flags
    attr(flag_rows[[length(flag_rows)]], "glm_failed") <- glm_failed
  }

  res <- do.call(rbind, rows)
  tst <- res$tested

  ## BH per trait per layer, over the family actually tested
  res$q_logistic <- NA_real_
  fam_g <- tst & !is.na(res$p_logistic)
  res$q_logistic[fam_g] <- p.adjust(res$p_logistic[fam_g], "BH")
  m_logistic <- sum(fam_g)

  res$q_rphylo <- NA_real_
  fam_r <- tst & !is.na(res$p_rphylo)
  res$q_rphylo[fam_r] <- p.adjust(res$p_rphylo[fam_r], "BH")
  m_rphylo <- sum(fam_r)
  res$m_logistic <- m_logistic; res$m_rphylo <- m_rphylo

  ## disagreement flag: XOR of significance, or both significant with opposite signs
  sig_l <- !is.na(res$q_logistic) & res$q_logistic < 0.05
  sig_r <- !is.na(res$q_rphylo) & res$q_rphylo < 0.05
  opp <- sig_l & sig_r & !is.na(res$OR) & !is.na(res$r_phylo) &
    (sign(res$OR - 1) != sign(res$r_phylo))
  opp[is.na(opp)] <- FALSE
  disagree <- xor(sig_l, sig_r) | opp
  disagree[is.na(disagree)] <- FALSE

  ## assemble result flags into the flags file (this trait's rows only)
  for (kk in seq_along(mods)) {
    k <- flag_start + kk - 1L
    fr <- flag_rows[[k]]
    if (fr$tested[1] != "yes") next
    rf <- character(0); d <- character(0)
    gf <- attr(fr, "glm_flags"); lf <- attr(fr, "lm_flags")
    if (!is.null(gf)) for (fn in names(gf)) if (gf[fn]) rf <- c(rf, flag_names[fn])
    if (isTRUE(attr(fr, "glm_failed"))) rf <- c(rf, "logistic fit failed")
    for (l in lf) rf <- c(rf, l)
    if (disagree[kk]) {
      rf <- c(rf, "logistic/r_phylo disagreement")
      d <- c(d, sprintf("logistic q=%.3g (OR=%.3g), r_phylo q=%.3g (r=%.3g)",
                        res$q_logistic[kk], res$OR[kk], res$q_rphylo[kk], res$r_phylo[kk]))
    }
    if (isTRUE(res$flag_persisted[kk])) rf <- c(rf, "fit flag persisted after refit")
    fr$result_flags <- paste(unique(rf), collapse = "; ")
    fr$detail <- paste(Filter(function(z) nchar(z) > 0, c(fr$detail, d)), collapse = "; ")
    flag_rows[[k]] <- fr
  }

  res$any_flag <- vapply(seq_along(mods), function(kk) {
    fr <- flag_rows[[flag_start + kk - 1L]]
    if (fr$tested[1] == "yes") nchar(fr$result_flags) > 0 else NA
  }, TRUE)
  res$any_flag <- ifelse(tst, ifelse(res$any_flag %in% TRUE, "yes", "no"), "no")

  res <- res[, c("modification", "n", "n_class0", "n_class1", "n_present_class0",
                 "n_present_class1", "OR", "CI_low", "CI_high", "p_logistic", "q_logistic",
                 "alpha", "phi", "r_phylo", "p_rphylo", "q_rphylo", "lambda",
                 "naive_r", "naive_p", "m_logistic", "m_rphylo",
                 "tested", "refitted", "flag_persisted", "any_flag")]
  write.table(res, file.path(tdir, paste0(tr, "_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  res_list[[tr]] <- res

  cat(sprintf("tested: %d/36 | BH family m: logistic=%d, r_phylo=%d | significant (q<0.05): logistic=%d\n",
              sum(tst), m_logistic, m_rphylo, sum(res$q_logistic < 0.05, na.rm = TRUE)))

  ## ---- plot: r_phylo per tested modification, marked by logistic q < 0.05
  pd <- res[tst & !is.na(res$r_phylo), ]
  pd$sig <- ifelse(!is.na(pd$q_logistic) & pd$q_logistic < 0.05, "q < 0.05", "q >= 0.05")
  pd$modification <- factor(pd$modification, levels = pd$modification[order(pd$r_phylo)])
  p <- ggplot(pd, aes(x = r_phylo, y = modification, fill = sig)) +
    geom_col(width = 0.7, colour = "grey25", linewidth = 0.2) +
    scale_fill_manual(values = c("q < 0.05" = "#0279EE", "q >= 0.05" = "#ECE9E2"), name = "logistic") +
    labs(x = "r_phylo", y = NULL, title = tr) +
    theme_bw(base_size = 11) + theme(text = element_text(family = "Liberation Sans"),
                                     legend.position = "top")
  ggsave(file.path(tdir, paste0(tr, "_r_phylo_plot.png")), p, width = 6.5, height = 5.5, dpi = 300)
  ggsave(file.path(tdir, paste0(tr, "_r_phylo_plot.svg")), p, width = 6.5, height = 5.5)
}

## ---------------------------------------------------------------- flags file
flags <- do.call(rbind, flag_rows)
flags <- flags[, c("trait", "modification", "tested", "untested_reasons", "result_flags", "detail")]
write.table(flags, file.path(outdir, "flags.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")

cat("\n== package versions ==\n")
cat(R.version.string, "\n")
for (pk in c("phylolm", "ape", "phytools", "ggplot2"))
  cat(sprintf("%s: %s\n", pk, as.character(packageVersion(pk))))
cat("\nDone. Outputs in", outdir, "\n")
