#!/usr/bin/env Rscript
################################################################################
# rga_core_shell_cloud.R
# Pan-RGAome: orthogroup (OG) prevalence, core/shell/cloud classes, RGA-class
# assignment, class-by-category test and P. sativum-specific / shared OGs.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Construction of Pan-RGAome
# Figures    : Fig. 2c (prevalence density + chi-square test), Fig. 2d (taxon sharing)
#
# Input
#   --clusters  MMseqs2 cluster table, tab-separated, no header:
#               sample  gene_id  gene_id  cluster_id  cluster_size
#   --rga       RGA annotation, tab-separated, no header:
#               gene_id  length  RGA_class  ...
#   --taxa      (optional) sample<TAB>taxon; needed for Fig. 2d
#
# Options (defaults)
#   --min-samples 3        keep OGs present in >= this many samples
#   --cloud-max   18       cloud = min-samples .. cloud-max
#   --shell-max   63       shell = cloud-max+1 .. shell-max; core = > shell-max
#   --min-prop    0.5      OG gets an RGA class if >= this fraction of its genes have it
#   --psat-taxon  "P. sativum"   taxon label(s) counted as P. sativum (comma separated)
#   --B           100000   Monte Carlo replicates for the chi-square test
#   --seed        42
#   --out         rga      output prefix
#
# Example
#   Rscript rga_core_shell_cloud.R --clusters pisum85_cluster.tsv --rga pisum85_RGA.tsv \
#       --taxa pisum85_taxa.tsv --out Pisum85_RGA
#
# Outputs
#   <out>_OG_summary.tsv          one row per OG: n_samples, category, RGA class (or Unassigned)
#   <out>_OG_class_assignments.tsv OG-class pairs with proportion >= min-prop
#   <out>_category_by_class.tsv   counts used for the chi-square test + test result
#   <out>_taxon_sharing.tsv       Psat-specific / shared / non-Psat OGs by class (with --taxa)
#   <out>_prevalence_density.pdf  Fig. 2c
#   <out>_taxon_sharing.pdf       Fig. 2d (with --taxa)
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

#---------------------------- arguments ----------------------------------------
defaults <- list(clusters = NA, rga = NA, taxa = NA, `min-samples` = 3, `cloud-max` = 18,
                 `shell-max` = 63, `min-prop` = 0.5, `psat-taxon` = "P. sativum",
                 B = 1e5, seed = 42, out = "rga")
parse_args <- function(a, d) {
  if (length(a) %% 2 != 0 || !all(startsWith(a[c(TRUE, FALSE)], "--")))
    stop("Arguments must be --name value pairs. See header of this script.")
  for (i in seq(1, length(a), 2)) {
    k <- sub("^--", "", a[i]); if (!k %in% names(d)) stop("Unknown option --", k)
    d[[k]] <- a[i + 1]
  }
  d
}
opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
if (is.na(opt$clusters) || is.na(opt$rga)) stop("--clusters and --rga are required")
min_s     <- as.integer(opt$`min-samples`)
cloud_max <- as.integer(opt$`cloud-max`)
shell_max <- as.integer(opt$`shell-max`)
min_prop  <- as.numeric(opt$`min-prop`)
psat_lab  <- trimws(strsplit(opt$`psat-taxon`, ",")[[1]])
set.seed(as.integer(opt$seed))
out <- opt$out
log <- function(...) cat(sprintf(...), "\n")

#---------------------------- 1. clusters and OG prevalence --------------------
cl <- fread(opt$clusters, header = FALSE, sep = "\t")
cl <- unique(cl[, .(sample = V1, gene = V2, og = V4)])
n_total <- uniqueN(cl$sample)
log("Samples: %d | OGs (all): %d", n_total, uniqueN(cl$og))

og_tab <- cl[, .(n_samples = uniqueN(sample), n_genes = .N), by = og][n_samples >= min_s]
cl <- cl[og %in% og_tab$og]
og_tab[, category := fcase(n_samples <= cloud_max, "Cloud",
                       n_samples <= shell_max, "Shell",
                       default = "Core")]
og_tab[, category := factor(category, levels = c("Cloud", "Shell", "Core"))]
log("OGs in >= %d samples: %d (cloud %d | shell %d | core %d)", min_s, nrow(og_tab),
    og_tab[category == "Cloud", .N], og_tab[category == "Shell", .N], og_tab[category == "Core", .N])

#---------------------------- 2. RGA class per OG ------------------------------
rga <- fread(opt$rga, header = FALSE, sep = "\t")
rga <- unique(rga[, .(gene = V1, rga_class = trimws(as.character(V3)))])[nzchar(rga_class)]
rga[rga_class %in% c("NBS>TNL", "NBS>CNL", "NBS>RNL"), rga_class := sub("^NBS>", "", rga_class)]
log("RGA classes: %s", paste(sort(unique(rga$rga_class)), collapse = ", "))

x <- merge(cl, rga, by = "gene")
cls <- x[, .(n_class_genes = .N), by = .(og, rga_class)]
cls[, n_annotated_genes := sum(n_class_genes), by = og]
cls[, proportion := n_class_genes / n_annotated_genes]
assign <- merge(cls[proportion >= min_prop], og_tab, by = "og")
log("OGs assigned to more than one class (ties at %.0f%%): %d", 100 * min_prop,
    assign[, .N, by = og][N > 1, .N])

# One class per OG for summaries (highest proportion; ties -> alphabetical)
best <- assign[order(og, -proportion, rga_class)][, .SD[1], by = og]
og_sum <- merge(og_tab, best[, .(og, rga_class, proportion)], by = "og", all.x = TRUE)
og_sum[is.na(rga_class), rga_class := "Unassigned"]
fwrite(og_sum, paste0(out, "_OG_summary.tsv"), sep = "\t")
fwrite(assign, paste0(out, "_OG_class_assignments.tsv"), sep = "\t")

#---------------------------- 3. class composition vs category (Fig. 2c test) --
test_classes <- c("TNL", "CNL", "RNL", "RLK", "RLP")
tab <- table(droplevels(og_sum[rga_class %in% test_classes, .(rga_class, category)]))
chi <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = as.numeric(opt$B)))
log("Monte Carlo chi-square (B = %d): X2 = %.1f, P = %.2g", as.integer(opt$B), chi$statistic, chi$p.value)
tab_dt <- as.data.table(tab)
setnames(tab_dt, "N", "n_OGs")
tab_dt[, pct_within_class := round(100 * n_OGs / sum(n_OGs), 1), by = rga_class]
fwrite(tab_dt, paste0(out, "_category_by_class.tsv"), sep = "\t")
cat(sprintf("# Monte Carlo chi-square: X2=%.2f; P=%.3g; B=%d\n", chi$statistic, chi$p.value,
            as.integer(opt$B)), file = paste0(out, "_category_by_class.tsv"), append = TRUE)

#---------------------------- 4. prevalence density (Fig. 2c) ------------------
plot_classes <- c("All", test_classes)
dens_in <- rbind(og_tab[, .(rga_class = "All", n_samples)],
                 assign[rga_class %in% test_classes, .(rga_class, n_samples)])
bw <- stats::bw.nrd0(og_tab$n_samples); if (!is.finite(bw) || bw <= 0) bw <- 1
curves <- dens_in[, {
  if (.N < 2) list(x = numeric(0), y = numeric(0)) else {
    d <- density(n_samples, bw = bw, adjust = 1.5, n = 2048, from = min_s, to = n_total)
    list(x = d$x, y = d$y) }
}, by = rga_class]
n_lab <- dens_in[, .N, by = rga_class]
curves[, rga_class := factor(rga_class, levels = plot_classes,
         labels = paste0(plot_classes, " (", n_lab[match(plot_classes, rga_class), N], ")"))]
cols <- setNames(c("#333333", "#2878CF", "#F16A3A", "#1DB18B", "#9A54B4", "#D4A127"),
                 levels(curves$rga_class))
regions <- data.table(label = c("Cloud", "Shell", "Core"),
                      x = c((min_s + cloud_max) / 2, (cloud_max + shell_max) / 2, (shell_max + n_total) / 2))

p_dens <- ggplot(curves, aes(x, y, colour = rga_class)) +
  geom_line(data = curves[grepl("^All", rga_class)], linetype = "dashed", linewidth = 1) +
  geom_line(data = curves[!grepl("^All", rga_class)], linewidth = 1) +
  geom_vline(xintercept = c(cloud_max, shell_max) + 0.5, linetype = "dashed", colour = "grey50") +
  geom_text(data = regions, aes(x = x, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.5, size = 4.5) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = sort(unique(c(pretty(c(min_s, n_total), 8), n_total))),
                     expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(min_s, n_total)) +
  labs(x = "Number of assemblies with orthogroup", y = "Density", colour = NULL,
       caption = sprintf("Monte Carlo chi-square = %.1f, P = %.2g", chi$statistic, chi$p.value)) +
  theme_classic(base_size = 13) + theme(legend.position = c(0.5, 0.75))
ggsave(paste0(out, "_prevalence_density.pdf"), p_dens, width = 7, height = 5)

#---------------------------- 5. P. sativum-specific / shared (Fig. 2d) --------
if (!is.na(opt$taxa)) {
  tx <- fread(opt$taxa, header = FALSE, sep = "\t")[, .(sample = V1, taxon = trimws(V2))]
  miss <- setdiff(unique(cl$sample), tx$sample)
  if (length(miss)) stop("Samples missing from --taxa: ", paste(head(miss, 5), collapse = ", "))
  tx[, is_psat := taxon %in% psat_lab]
  share <- merge(unique(cl[, .(og, sample)]), tx, by = "sample")[
    , .(psat = any(is_psat), other = any(!is_psat)), by = og]
  share[, group := fcase(psat & !other, "P. sativum-specific",
                         psat & other, "Shared",
                         default = "Other Pisum-specific")]
  share <- merge(share, og_sum[, .(og, rga_class)], by = "og")
  sh_tab <- share[, .N, by = .(group, rga_class)][order(group, -N)]
  fwrite(sh_tab, paste0(out, "_taxon_sharing.tsv"), sep = "\t")
  log("P. sativum-specific: %d | shared: %d | other Pisum-specific: %d",
      share[group == "P. sativum-specific", .N], share[group == "Shared", .N],
      share[group == "Other Pisum-specific", .N])
  sh_tab[, group := factor(group, levels = c("P. sativum-specific", "Shared", "Other Pisum-specific"))]
  p_sh <- ggplot(sh_tab, aes(group, N, fill = rga_class)) +
    geom_col(width = 0.7) +
    geom_text(data = sh_tab[, .(N = sum(N)), by = group], aes(group, N, label = N),
              inherit.aes = FALSE, vjust = -0.4) +
    labs(x = NULL, y = "RGA orthogroups", fill = "RGA class") +
    theme_classic(base_size = 13)
  ggsave(paste0(out, "_taxon_sharing.pdf"), p_sh, width = 6, height = 5)
}
log("Done. Outputs written with prefix: %s", out)