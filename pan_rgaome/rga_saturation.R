#!/usr/bin/env Rscript
################################################################################
# rga_saturation.R
# Pan-RGAome saturation: how many assemblies are needed to recover RGA
# orthogroups (OGs). Assemblies are added in random order (permutations) and
# the cumulative number of OGs is recorded at each step.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Construction of Pan-RGAome (saturation analysis, 500 permutations)
# Figure     : Fig. 2b
#
# Input
#   --clusters  MMseqs2 cluster table, tab-separated, no header:
#               sample  gene_id  gene_id  cluster_id  cluster_size
#
# Options (defaults)
#   --min-samples 3      keep OGs present in >= this many samples
#   --perm        500    number of random assembly orders
#   --target      95     % of OGs used to report the saturation point
#   --thin        5      x-axis step for the thinned plot
#   --seed        42
#   --out         rga    output prefix
#
# Example
#   Rscript rga_saturation.R --clusters pisum85_cluster.tsv --perm 500 --out Pisum85_RGA
#
# Outputs
#   <out>_saturation_permutations.tsv.gz  perm, step, cum_OGs, cum_pct, remaining
#   <out>_saturation_summary.tsv          per step: median / 2.5% / 97.5% of cum_pct
#   <out>_saturation_all.pdf              boxplot at every step
#   <out>_saturation_thin.pdf             boxplot every --thin steps (Fig. 2b)
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

#---------------------------- arguments ----------------------------------------
defaults <- list(clusters = NA, `min-samples` = 3, perm = 500, target = 95, thin = 5,
                 seed = 42, out = "rga")
a <- commandArgs(trailingOnly = TRUE)
if (length(a) %% 2 != 0) stop("Arguments must be --name value pairs. See header of this script.")
for (i in seq(1, length(a), 2)) {
  k <- sub("^--", "", a[i]); if (!k %in% names(defaults)) stop("Unknown option --", k)
  defaults[[k]] <- a[i + 1]
}
opt <- defaults
if (is.na(opt$clusters)) stop("--clusters is required")
min_s  <- as.integer(opt$`min-samples`); n_perm <- as.integer(opt$perm)
target <- as.numeric(opt$target);        thin   <- as.integer(opt$thin)
out    <- opt$out
set.seed(as.integer(opt$seed))

#---------------------------- 1. load and filter OGs ---------------------------
dt <- fread(opt$clusters, header = FALSE, sep = "\t")
dt <- unique(dt[, .(sample = V1, og = V4)])
og_n <- dt[, .(n_samples = uniqueN(sample)), by = og]
dt <- dt[og %in% og_n[n_samples >= min_s, og]]

samples <- sort(unique(dt$sample))
n_samp  <- length(samples)
n_og    <- uniqueN(dt$og)
cat(sprintf("Assemblies: %d | OGs in >= %d assemblies: %d\n", n_samp, min_s, n_og))

#---------------------------- 2. permutations ----------------------------------
# OG membership per sample as integer indices -> fast cumulative union
dt[, og_idx := as.integer(factor(og))]
og_by_sample <- split(dt$og_idx, dt$sample)[samples]

one_perm <- function() {
  ord <- sample.int(n_samp)
  seen <- logical(n_og); cum <- integer(n_samp)
  for (i in seq_len(n_samp)) {
    seen[og_by_sample[[ord[i]]]] <- TRUE
    cum[i] <- sum(seen)
  }
  cum
}
res <- rbindlist(lapply(seq_len(n_perm), function(p)
  data.table(perm = p, step = seq_len(n_samp), cum_OGs = one_perm())))
res[, `:=`(cum_pct = 100 * cum_OGs / n_og, remaining = n_og - cum_OGs)]
fwrite(res, paste0(out, "_saturation_permutations.tsv.gz"), sep = "\t")

summ <- res[, .(median_pct = median(cum_pct), lo_pct = quantile(cum_pct, 0.025),
                hi_pct = quantile(cum_pct, 0.975), median_remaining = median(remaining)), by = step]
fwrite(summ, paste0(out, "_saturation_summary.tsv"), sep = "\t")
step_t <- summ[median_pct >= target, step][1]
cat(sprintf("Median %.0f%% of OGs recovered with %s assemblies (%d permutations)\n",
            target, ifelse(is.na(step_t), "NA", step_t), n_perm))

#---------------------------- 3. plots -----------------------------------------
scale2 <- 100 / n_og
sat_plot <- function(d, breaks, xlab) {
  ggplot(d, aes(step, group = step)) +
    geom_boxplot(aes(y = cum_pct), width = 0.6, colour = "black", fill = NA, outlier.size = 0.2) +
    geom_boxplot(aes(y = remaining * scale2), width = 0.35, colour = "#2C6FB7", fill = NA, outlier.size = 0.2) +
    geom_hline(yintercept = target, linetype = "dashed", linewidth = 0.6) +
    { if (!is.na(step_t)) geom_vline(xintercept = step_t, linetype = "dashed", linewidth = 0.6) } +
    { if (!is.na(step_t)) annotate("text", x = step_t, y = 5, hjust = -0.1,
                                   label = paste0(target, "% at ", step_t, " assemblies")) } +
    scale_y_continuous(name = "Orthogroups recovered (%)",
                       sec.axis = sec_axis(~ . / scale2, name = "Orthogroups remaining")) +
    scale_x_continuous(breaks = breaks, limits = c(0.5, n_samp + 0.5), expand = expansion(mult = 0.01)) +
    labs(x = xlab) +
    theme_classic(base_size = 13) +
    theme(axis.title.y.right = element_text(colour = "#2C6FB7"), axis.text.y.right = element_text(colour = "#2C6FB7"))
}
keep <- sort(unique(c(1, seq(thin, n_samp, by = thin), n_samp)))
ggsave(paste0(out, "_saturation_all.pdf"),  sat_plot(res, pretty(c(1, n_samp), 10), "Number of assemblies"), width = 7.5, height = 6)
ggsave(paste0(out, "_saturation_thin.pdf"), sat_plot(res[step %in% keep], keep, "Number of assemblies"), width = 7.5, height = 6)
cat("Done. Outputs written with prefix:", out, "\n")