## Modification-modification phylogenetic correlation -- Aquificota outgroup rooting
## Design: 630 unordered pairs among 36 modifications, all 1322 genomes, one midpoint-rooted tree.
## Presence filter: FRACTION RULE (< 5% or > 95% of 1,322) per user decision 2026-11-14.
suppressMessages({library(phylolm); library(ape); library(phytools)})

indir  <- "/mnt/user-uploads/input_v4"
outdir <- "/mnt/results/mod_mod_correlation_aquificota_root"
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

##### 1. QC #####
mod <- read.delim(file.path(indir, "mod_AP.txt"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(mod) == 37)                       # Assembly ID + 36 modification columns
rownames(mod) <- mod$Assembly
m <- as.matrix(mod[, -1, drop = FALSE])
storage.mode(m) <- "numeric"
stopifnot(nrow(m) == 1322)                       # 1.1
stopifnot(all(m %in% c(0, 1)), !any(is.na(m)))   # 1.1
tree <- read.tree(file.path(indir, "bac120.fasttree.unrooted.tree"))
stopifnot(setequal(rownames(m), tree$tip.label), anyDuplicated(rownames(m)) == 0)  # 1.2

##### 1.3 one tree, rooted once at Aquificota (outgroup rooting, per user request) #####
traits <- read.delim(file.path(indir, "genome_traits.txt"), check.names = FALSE,
                     stringsAsFactors = FALSE, na.strings = c("", "NA"))
rownames(traits) <- traits$Assembly
root_phylum <- trimws(readLines("/mnt/user-uploads/root.txt")[1])
stopifnot("Phylum" %in% colnames(traits))
aqu <- rownames(traits)[!is.na(traits$Phylum) & traits$Phylum == root_phylum]
stopifnot(length(aqu) > 0, all(aqu %in% tree$tip.label), is.monophyletic(tree, aqu))
cat(sprintf("outgroup '%s': %d genomes, monophyletic: TRUE\n", root_phylum, length(aqu)))
phy <- root(tree, outgroup = aqu, resolve.root = TRUE)
stopifnot(setequal(phy$tip.label, rownames(m)))
write.tree(phy, file.path(outdir, "tree_aquificota_rooted.nwk"))
n <- nrow(m)
Ttree <- mean(node.depth.edgelength(phy)[1:n])
alpha_lo <- -log(Ttree) - 4; alpha_hi <- -log(Ttree) + 4

##### 1.4 filters #####
pres <- colSums(m) / n
pass <- colnames(m)[pres >= 0.05 & pres <= 0.95]          # fraction rule
excluded <- setdiff(colnames(m), pass)
mods_sorted <- sort(pass, method = "radix")               # C-locale byte order (1.4/2.1)
all_pairs <- combn(mods_sorted, 2)                        # pr[1,k]=predictor (sorts first), pr[2,k]=response
npairs_total <- ncol(all_pairs)

zero_cell <- function(x, y) {
  tb <- table(factor(x, c(0, 1)), factor(y, c(0, 1)))
  if (any(tb == 0)) {
    idx <- which(tb == 0, arr.ind = TRUE)[1, ]
    paste0("zero cell: predictor=", rownames(tb)[idx[1]], " response=", colnames(tb)[idx[2]],
           " (count=0; predictor=", colnames(tb)[idx[2]], "=0/1 x response=", rownames(tb)[idx[1]], "=0/1)")
  } else NA_character_
}

##### helpers #####
sanitize <- function(x) gsub("[\r\n\t]+", " ", x)
capture_fit <- function(expr) {
  warns <- character(0)
  val <- tryCatch(
    withCallingHandlers(expr,
      warning = function(w) { warns <<- c(warns, sanitize(conditionMessage(w))); invokeRestart("muffleWarning") }),
    error = function(e) { attr(warns, "error") <- sanitize(conditionMessage(e)); NULL })
  if (!is.null(val) && inherits(val, "try-error")) val <- NULL  # try() inside: convert to error path
  attr(val, "warnings") <- warns
  val
}
fit_flags <- function(f, warns) {
  # flag conditions per 2.3, evaluated on a fitted phyloglm object
  fl <- character(0)
  if (!is.null(f)) {
    cf <- stats::coefficients(f); se <- summary(f)$coefficients[, "StdErr"]
    if (is.null(f$convergence) || f$convergence != 0) fl <- c(fl, "non-convergence")
    if (is.null(f$alphaWarn) || f$alphaWarn != 0)     fl <- c(fl, "alphaWarn!=0")
    if (any(!is.finite(cf)) || any(!is.finite(se)))   fl <- c(fl, "missing coefficient or SE")
    if (all(is.finite(cf)) && any(abs(cf) > 5))       fl <- c(fl, "|coef|>5")
  }
  if (length(warns)) fl <- c(fl, "phyloglm warning")
  unique(fl)
}

##### pair loop with resume #####
rawfile <- file.path(outdir, "pairs_raw.csv")
done <- character(0)
if (file.exists(rawfile)) {
  prev <- read.delim(rawfile, stringsAsFactors = FALSE)
  done <- paste(prev$predictor, prev$response, sep = "|")
  cat("resuming:", length(done), "pairs already computed\n")
}
raw <- file(rawfile, open = "a")
if (!file.exists(rawfile) || file.info(rawfile)$size == 0) {
  writeLines(paste(c("predictor","response","n","n00","n10","n01","n11",
    "OR","OR_lower","OR_upper","p","alpha","refitted","fit_flags_kept","refit_detail",
    "r_phylo","r_phylo_p","lambda","lambda_warn","phi","phylom_error","phylom_warn","status","status_detail"),
    collapse = "\t"), raw)
}

t0 <- Sys.time()
for (k in seq_len(npairs_total)) {
  pi_ <- all_pairs[1, k]; pj <- all_pairs[2, k]
  key <- paste(pi_, pj, sep = "|")
  if (key %in% done) next
  x <- m[, pi_]; y <- m[, pj]
  d <- data.frame(y = y, x = x); rownames(d) <- rownames(m)
  zc <- zero_cell(x, y)
  if (!is.na(zc)) {  # untested: zero cell
    writeLines(paste(c(pi_, pj, n, rep(NA, 4), rep(NA, 5), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "untested_zero_cell", zc),
                     collapse = "\t"), raw); flush(raw); next
  }
  ## 2.2 primary phyloglm (default start) -- candidate 1
  f1 <- capture_fit(phyloglm(y ~ x, data = d, phy = phy, method = "logistic_MPLE", btol = 10))
  fl1 <- if (is.null(f1)) character(0) else fit_flags(f1, attr(f1, "warnings"))
  err1 <- attr(f1, "error")
  candidates <- list(list(f = f1, flags = fl1, start = "default"))
  refit_detail <- ""
  ## 2.3 refit if flagged (or failed): 9 random starts; set.seed(123) immediately before draws
  if (length(fl1) > 0 || !is.null(err1)) {
    set.seed(123)
    for (r in 1:9) {
      sb <- runif(2, -2, 2)                       # coefficients first, uniform(-2, 2)
      sa <- exp(runif(1, alpha_lo, alpha_hi))     # then alpha: uniform on the log scale (search range is for log alpha)
      fr <- capture_fit(phyloglm(y ~ x, data = d, phy = phy, method = "logistic_MPLE",
                                 btol = 10, start.beta = sb, start.alpha = sa))
      if (is.null(fr)) next
      flr <- fit_flags(fr, attr(fr, "warnings"))
      candidates <- c(candidates, list(list(f = fr, flags = flr, start = sprintf("random%d", r))))
    }
    refit_detail <- paste0(sapply(candidates, function(cc) if (is.null(cc$f)) paste0(cc$start, ":error") else
      paste0(cc$start, ":code=", cc$f$convergence, ":penlogLik=", sprintf("%.3f", cc$f$penlogLik))), collapse = "; ")
  }
  ## keep highest penlogLik among converged fits with computable Wald
  ok <- Filter(function(cc) !is.null(cc$f) && cc$f$convergence == 0 &&
                 all(is.finite(stats::coefficients(cc$f))) &&
                 all(is.finite(summary(cc$f)$coefficients[, "StdErr"])), candidates)
  status <- "ok"; status_detail <- ""
  if (length(ok) == 0) {
    status <- "failed"
    status_detail <- if (!is.null(err1)) paste0("phyloglm error: ", err1) else "no starting value converged / Wald not computable"
    OR <- ORl <- ORu <- p <- alpha <- NA_real_; refitted <- "yes"; kept_flags <- NA_character_
  } else {
    pen <- sapply(ok, function(cc) cc$f$penlogLik)
    keep <- ok[[which.max(pen)]]
    f <- keep$f
    co <- summary(f)$coefficients
    est <- co["x", "Estimate"]; se <- co["x", "StdErr"]; p <- co["x", "p.value"]
    OR <- exp(est); ORl <- exp(est - 1.96 * se); ORu <- exp(est + 1.96 * se)
    alpha <- f$alpha
    refitted <- if (length(candidates) > 1) "yes" else "no"
    kept_flags <- if (length(keep$flags)) paste(keep$flags, collapse = "; ") else NA_character_
  }
  ## 2.4 phylolm lambda
  pl <- capture_fit(phylolm(y ~ x, data = d, phy = phy, model = "lambda",
                            lower.bound = 1e-7, upper.bound = 1))
  plerr <- attr(pl, "error"); plwarn <- attr(pl, "warnings")
  if (is.null(pl)) {
    r_phylo <- r_p <- lambda <- NA_real_
    lw <- paste0("phylolm error: ", plerr)
  } else {
    co2 <- summary(pl)$coefficients
    tval <- co2["x", "t.value"]; r_phylo <- tval / sqrt(tval^2 + n - 2); r_p <- co2["x", "p.value"]
    lambda <- pl$optpar
    lw <- if (length(plwarn)) paste(plwarn, collapse = "; ") else NA_character_
    if (abs(lambda) < 1e-6 || abs(lambda - 1) < 1e-6)
      lw <- paste(c(lw, sprintf("lambda within 1e-6 of %s", ifelse(lambda < 0.5, "0", "1"))), collapse = "; ")
  }
  ## 2.5 naive Pearson
  phi <- cor(x, y)
  writeLines(paste(c(pi_, pj, n,
    sum(x == 0 & y == 0), sum(x == 1 & y == 0), sum(x == 0 & y == 1), sum(x == 1 & y == 1),
    OR, ORl, ORu, p, alpha, refitted, kept_flags, refit_detail,
    r_phylo, r_p, lambda, lw, phi, if (is.null(pl)) plerr else NA_character_,
    if (length(plwarn)) paste(plwarn, collapse = "; ") else NA_character_,
    status, status_detail), collapse = "\t"), raw)
  flush(raw)
  if (k %% 25 == 0) cat(format(Sys.time(), "%H:%M:%S"), "pair", k, "/", npairs_total,
                        "elapsed", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
close(raw)
cat("pair loop done:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

##### assemble outputs #####
raw <- read.delim(rawfile, stringsAsFactors = FALSE)
tz <- raw[raw$status == "untested_zero_cell", ]
tt <- raw[raw$status %in% c("ok", "failed"), ]
tt$p <- as.numeric(tt$p); tt$r_phylo_p <- as.numeric(tt$r_phylo_p)

## 2.2 BH across tested pairs with non-missing p (failed fits excluded)
m_phy <- sum(!is.na(tt$p))
tt$q <- p.adjust(tt$p, method = "BH")
m_rph <- sum(!is.na(tt$r_phylo_p))
tt$r_phylo_q <- p.adjust(tt$r_phylo_p, method = "BH")

## disagreement flag (3.3)
sigL <- !is.na(tt$q) & tt$q < 0.05; sigR <- !is.na(tt$r_phylo_q) & tt$r_phylo_q < 0.05
opp  <- sigL & sigR & ((!is.na(tt$OR) & tt$OR > 1 & tt$r_phylo < 0) | (!is.na(tt$OR) & tt$OR < 1 & tt$r_phylo > 0))
tt$disagreement <- ifelse((sigL != sigR) | opp, "yes", "no")

## 3.1 pair table
pt <- tt[, c("predictor","response","n","n00","n10","n01","n11","OR","OR_lower","OR_upper",
             "p","q","alpha","r_phylo","r_phylo_p","r_phylo_q","lambda","phi","refitted",
             "fit_flags_kept","disagreement")]
names(pt)[names(pt) == "fit_flags_kept"] <- "fit_flag_persisted"
pt$any_flags <- ifelse(apply(pt[, c("fit_flag_persisted","disagreement")], 1,
                             function(r) any(!is.na(r) & r != "no")), "yes", "no")
write.table(pt, file.path(outdir, "pair_table.tsv"), sep = "\t", quote = FALSE, na = "", row.names = FALSE)

## 3.2 matrices (n x n, tested modifications; cell [predictor, response] = upper triangle filled)
mkmat <- function(valcol) {
  A <- matrix(NA_real_, length(mods_sorted), length(mods_sorted), dimnames = list(mods_sorted, mods_sorted))
  v <- tt[[valcol]]; v <- as.numeric(v)
  for (i in seq_len(nrow(tt))) A[tt$predictor[i], tt$response[i]] <- v[i]
  A
}
Rmat <- mkmat("r_phylo"); Pmat <- mkmat("phi"); Qmat <- mkmat("q")
write.table(Rmat, file.path(outdir, "r_phylo_matrix.tsv"), sep = "\t", quote = FALSE, na = "NA")
write.table(Pmat, file.path(outdir, "phi_matrix.tsv"),       sep = "\t", quote = FALSE, na = "NA")
write.table(Qmat, file.path(outdir, "phyloglm_q_matrix.tsv"),sep = "\t", quote = FALSE, na = "NA")

## 3.3 flags file: all 630 pairs
fp <- data.frame(predictor = all_pairs[1, ], response = all_pairs[2, ], stringsAsFactors = FALSE)
fp$tested <- ifelse(paste(fp$predictor, fp$response, sep = "|") %in% paste(tt$predictor, tt$response, sep = "|"), "yes", "no")
fp$untested_reason <- NA_character_
fp$untested_reason[fp$tested == "no"] <- apply(fp[fp$tested == "no", c("predictor","response")], 1, function(r) {
  ex <- c(r[1], r[2]) %in% excluded
  zc <- tz$untested_reason[tz$predictor == r[1] & tz$response == r[2]]
  parts <- c(if (any(ex)) paste0("modification excluded by presence filter: ",
             paste(c(r[1], r[2])[ex], collapse = ", ")), if (length(zc)) zc)
  paste(parts, collapse = "; ")
})
fp$lambda_bound <- ifelse(fp$tested == "yes" & !is.na(tt$lambda_warn[match(paste(fp$predictor, fp$response), paste(tt$predictor, tt$response))]) &
                            grepl("lambda within", tt$lambda_warn[match(paste(fp$predictor, fp$response), paste(tt$predictor, tt$response))]), "yes", "no")
mi <- match(paste(fp$predictor, fp$response), paste(tt$predictor, tt$response))
fp$disagreement      <- ifelse(fp$tested == "yes", tt$disagreement[mi], NA)
fp$fit_flag_persisted<- ifelse(fp$tested == "yes" & !is.na(tt$fit_flags_kept[mi]) & tt$fit_flags_kept[mi] != "", "yes", "no")
fp$logistic_failed   <- ifelse(fp$tested == "yes" & tt$status[mi] == "failed", "yes", "no")
fp$phylom_flag       <- ifelse(fp$tested == "yes" & (!is.na(tt$phylom_warn[mi]) | !is.na(tt$phylom_error[mi])), "yes", "no")
fp$flag_detail <- apply(fp, 1, function(r) {
  parts <- c()
  if (!is.na(r["untested_reason"])) parts <- c(parts, r["untested_reason"])
  if (!is.na(r["fit_flag_persisted"]) && r["fit_flag_persisted"] == "yes")
    parts <- c(parts, paste0("fit flags: ", tt$fit_flags_kept[mi][which(paste(fp$predictor, fp$response) == paste(r["predictor"], r["response"]))]))
  if (!is.na(r["phylom_flag"]) && r["phylom_flag"] == "yes")
    parts <- c(parts, paste0("phylolm: ", na.omit(c(tt$phylom_error[mi], tt$phylom_warn[mi]))[1]))
  if (length(parts)) paste(parts, collapse = "; ") else NA_character_
})
fp$flags <- apply(fp[, c("untested_reason","lambda_bound","disagreement","fit_flag_persisted","logistic_failed","phylom_flag")], 1, function(r) {
  f <- c(); if (!is.na(r[1])) f <- c(f, "untested"); if (!is.na(r[2]) && r[2] == "yes") f <- c(f, "lambda_bound")
  if (!is.na(r[3]) && r[3] == "yes") f <- c(f, "disagreement"); if (!is.na(r[4]) && r[4] == "yes") f <- c(f, "fit_flag_persisted")
  if (!is.na(r[5]) && r[5] == "yes") f <- c(f, "logistic_failed"); if (!is.na(r[6]) && r[6] == "yes") f <- c(f, "phylolm_flag")
  if (length(f)) paste(f, collapse = ";") else NA_character_
})
write.table(fp[, c("predictor","response","tested","flags","untested_reason","lambda_bound",
                   "disagreement","fit_flag_persisted","logistic_failed","phylom_flag","flag_detail")],
            file.path(outdir, "flags.tsv"), sep = "\t", quote = FALSE, na = "", row.names = FALSE)

## 3.5 versions + summary
vinfo <- data.frame(package = c("R", "phylolm", "ape", "phytools", "ggplot2"),
                    version = c(paste(R.version$major, R.version$minor, sep = "."),
                                as.character(packageVersion("phylolm")), as.character(packageVersion("ape")),
                                as.character(packageVersion("phytools")), as.character(packageVersion("ggplot2"))))
write.table(vinfo, file.path(outdir, "package_versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(outdir, "summary.txt"))
cat("Modification-modification phylogenetic correlation (Aquificota outgroup rooting)\n")
cat("Presence filter: fraction rule (< 5% or > 95% of", n, "genomes)\n")
cat("Modifications passing:", length(mods_sorted), "| excluded:", paste(excluded, collapse = ", "), "\n")
cat("Pairs total:", npairs_total, "| zero-cell untested:", nrow(tz), "| tested:", nrow(tt),
    "| failed:", sum(tt$status == "failed"), "\n")
cat("BH family sizes: phyloglm m =", m_phy, "| r_phylo m =", m_rph, "\n")
cat("Significant (q < 0.05): phyloglm", sum(tt$q < 0.05, na.rm = TRUE), "| r_phylo", sum(tt$r_phylo_q < 0.05, na.rm = TRUE), "\n")
cat("Refitted:", sum(tt$refitted == "yes", na.rm = TRUE), "| fit flag persisted:", sum(pt$fit_flag_persisted == "yes", na.rm = TRUE),
    "| disagreements:", sum(tt$disagreement == "yes"), "\n")
sink()

##### 3.2 heatmap #####
sym <- Rmat; sym[lower.tri(sym)] <- t(sym)[lower.tri(sym)]  # mirror for display/clustering
Rs <- sym; diag(Rs) <- 1
D <- as.dist(1 - abs(Rs)); D[is.na(D)] <- 1
hc <- hclust(D, method = "complete")
ord <- hc$labels[hc$order]
disp <- Rs[ord, ord]
qm <- Qmat[ord, ord]; qm[lower.tri(qm)] <- t(qm)[lower.tri(qm)]
library(ggplot2)
df <- expand.grid(row = ord, col = ord, stringsAsFactors = FALSE)
df$r <- mapply(function(a, b) disp[a, b], df$row, df$col)
df$sig <- mapply(function(a, b) { v <- qm[a, b]; !is.na(v) && v < 0.05 }, df$row, df$col)
df$row <- factor(df$row, levels = rev(ord)); df$col <- factor(df$col, levels = ord)
p <- ggplot(df, aes(col, row, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(data = subset(df, sig), color = "black", linewidth = 0.7, fill = NA) +
  scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400", midpoint = 0,
                       na.value = "grey90", name = expression(r[phylo]), limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_family = "Liberation Sans") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank(),
        axis.title = element_blank()) +
  coord_fixed()
ggsave(file.path(outdir, "r_phylo_heatmap.png"), p, width = 9, height = 8, dpi = 200)
ggsave(file.path(outdir, "r_phylo_heatmap.svg"), p, width = 9, height = 8)
cat("heatmap done\n")
cat("ALL DONE\n")
