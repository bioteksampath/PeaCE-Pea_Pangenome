#!/usr/bin/env Rscript
################################################################################
# r1r7_expression_plot.R
# Expression of the Er3-linked NLR cluster genes (R1-R7) before and after
# powdery mildew inoculation: log2(TPM + 1), mean +/- SD of biological replicates.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Figure     : Fig. 4f
#
# Options (defaults)
#   --samples    sample sheet from rnaseq_qc.R                    (required)
#   --genes      gene_id<TAB>label (e.g. Pful...g0123<TAB>R1), header  (required)
#   --genotypes  order of genotypes in the plot   [G4,G9,IN9,IN16,IN33,IN36]
#   --labels     display names                    [G4=WL2140,G9=Terno]
#   --out        output prefix                    [Fig4f_R1-R7_expression]
#
# Outputs: <out>.pdf, <out>_values.tsv (per sample), <out>_summary.tsv (mean, SD, n)
################################################################################
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "rnaseq_utils.R"))
suppressPackageStartupMessages(library(ggplot2))

opt <- parse_opts(list(samples = NA, genes = NA, `rsem-dir` = "", genotypes = "G4,G9,IN9,IN16,IN33,IN36",
                       labels = "G4=WL2140,G9=Terno", out = "Fig4f_R1-R7_expression"))
lab <- label_map(opt$labels); gts <- split_list(opt$genotypes)
genes <- read.delim(opt$genes, stringsAsFactors = FALSE)[, 1:2]; names(genes) <- c("gene_id", "gene")
samples <- read_samples(opt$samples, opt$`rsem-dir`) %>% filter(genotype %in% gts)

tpm <- bind_rows(lapply(seq_len(nrow(samples)), function(i) {
  r <- read.delim(samples$rsem_file[i], stringsAsFactors = FALSE, check.names = FALSE)
  r %>% transmute(gene_id, TPM) %>% filter(base_id(gene_id) %in% base_id(genes$gene_id)) %>%
    mutate(sample = samples$sample[i])
})) %>%
  mutate(id0 = base_id(gene_id)) %>%
  left_join(genes %>% mutate(id0 = base_id(gene_id)) %>% select(id0, gene), by = "id0") %>%
  left_join(samples %>% select(sample, genotype, time, replicate), by = "sample") %>%
  mutate(log2TPM = log2(TPM + 1),
         time = factor(recode(time, PM0hr = "0 h", PM24hr = "24 h", PM4dpi = "4 d"), levels = c("0 h", "24 h", "4 d")),
         genotype = factor(relabel(genotype, lab), levels = relabel(gts, lab)),
         gene = factor(gene, levels = genes$gene))
miss <- setdiff(genes$gene, unique(as.character(tpm$gene)))
if (length(miss)) warning("Genes not found in RSEM output: ", paste(miss, collapse = ", "))

summ <- tpm %>% group_by(gene, genotype, time) %>%
  summarise(mean = mean(log2TPM), sd = sd(log2TPM), n = n(), .groups = "drop")
write_tsv(tpm %>% select(gene, gene_id, sample, genotype, time, replicate, TPM, log2TPM), paste0(opt$out, "_values.tsv"))
write_tsv(summ, paste0(opt$out, "_summary.tsv"))

p <- ggplot(summ, aes(time, mean, colour = genotype, group = genotype)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.15, linewidth = 0.4) +
  facet_wrap(~ gene, nrow = 1) +
  labs(x = "Time after inoculation", y = expression(log[2](TPM + 1)), colour = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(paste0(opt$out, ".pdf"), p, width = 2 + 1.6 * nlevels(summ$gene), height = 3.8)
message("Saved ", opt$out, ".pdf")
