## Sensitivity table + forest plot: phyloglm-significant pairs under midpoint vs Aquificota rooting
suppressMessages(library(ggplot2))
mid <- read.delim("/workspace/mod_mod_output/pair_table.tsv")
aqu <- read.delim("/workspace/mod_mod_output_aquificota/pair_table.tsv")
fmt <- function(d) data.frame(pair = paste0(d$predictor, " \u2192 ", d$response),
                              OR = d$OR, lo = d$OR_lower, hi = d$OR_upper, q = d$q,
                              r = d$r_phylo, rq = d$r_phylo_q, stringsAsFactors = FALSE)
m <- fmt(mid); a <- fmt(aqu)
sigm <- m[!is.na(m$q) & m$q < 0.05, ]; siga <- a[!is.na(a$q) & a$q < 0.05, ]
pairs_u <- union(sigm$pair, siga$pair)
mm <- m[match(pairs_u, m$pair), ]; aa <- a[match(pairs_u, a$pair), ]
tab <- data.frame(
  pair = pairs_u,
  OR_mid = sprintf("%.2f (%.2f\u2013%.2f)", mm$OR, mm$lo, mm$hi),
  q_mid = sprintf("%.2e", mm$q),
  OR_aqu = sprintf("%.2f (%.2f\u2013%.2f)", aa$OR, aa$lo, aa$hi),
  q_aqu = sprintf("%.2e", aa$q),
  r_phylo_mid = sprintf("%.3f", mm$r), r_phylo_aqu = sprintf("%.3f", aa$r),
  verdict = ifelse(!is.na(mm$q) & mm$q < 0.05 & !is.na(aa$q) & aa$q < 0.05,
                   "robust (both rootings)", "rooting-sensitive"),
  stringsAsFactors = FALSE)
tab <- tab[order(tab$verdict, tab$pair), ]
write.table(tab, "/mnt/results/mod_mod_correlation/sensitivity_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(tab, row.names = FALSE)

## forest plot: OR with 95% CI, log scale, two rootings side by side
pd <- rbind(cbind(mm, rooting = "Midpoint"), cbind(aa, rooting = "Aquificota"))
pd$rooting <- factor(pd$rooting, c("Midpoint", "Aquificota"))
pd$pair <- factor(pd$pair, levels = tab$pair)
pd$sig <- ifelse(pd$q < 0.05, "q < 0.05", "q \u2265 0.05")
pd$y <- as.numeric(pd$pair) + ifelse(pd$rooting == "Midpoint", -0.18, 0.18)  # manual y-dodge
p <- ggplot(pd, aes(OR, y, color = rooting, shape = sig, xmin = lo, xmax = hi)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_pointrange(linewidth = 0.55, fatten = 3) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq_along(tab$pair), labels = tab$pair, expand = expansion(add = 0.6)) +
  scale_color_manual(values = c(Midpoint = "#0279EE", Aquificota = "#FF9400"), name = "Rooting") +
  scale_shape_manual(values = c("q < 0.05" = 16, "q \u2265 0.05" = 1), name = "phyloglm") +
  labs(x = "Odds ratio (log scale, 95% CI)", y = NULL) +
  theme_minimal(base_family = "Liberation Sans", base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        axis.text.y = element_text(size = 11))
ggsave("/mnt/results/mod_mod_correlation/sensitivity_forest.png", p, width = 7.5, height = 4, dpi = 200)
ggsave("/mnt/results/mod_mod_correlation/sensitivity_forest.svg", p, width = 7.5, height = 4)
cat("figure done\n")
