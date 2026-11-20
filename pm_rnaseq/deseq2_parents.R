#!/usr/bin/env Rscript
################################################################################
# deseq2_parents.R
# Differential expression between the resistant donor (G4, P. fulvum WL2140)
# and the susceptible recurrent parent (G9, P. sativum Terno).
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Powdery mildew expression analyses (parent comparison)
# Figure     : Fig. 4f context; candidate tables for the Er3 interval
#
# Part A  Within each time point: G4 vs G9 at 0 h, 24 h and 4 dpi
#         design ~ condition (genotype_time); genes kept if >= 10 reads in >= 4 samples
# Part B  Infection-response interaction: (G4 24h - G4 0h) - (G9 24h - G9 0h)
#         design ~ genotype * time on 0 h and 24 h samples; >= 10 reads in >= 2 samples
#         (same test as the introgression-line analysis, deseq2_introgression_lines.R)
#
# DE: Benjamini-Hochberg padj < 0.05 and |log2FC| >= 1
#
# Options (defaults)
#   --samples      sample sheet from rnaseq_qc.R                    (required)
#   --annotation   GeneID Chr Start End RgeneFamily (tab, header)      (required)
#   --outdir       [DESeq2_Er3]
#   --resistant G4  --susceptible G9  --labels G4=WL2140,G9=Terno
#   --times        time points for Part A                              [PM0hr,PM24hr,PM4dpi]
#   --region       candidate interval                                  [chr4:8414176-19685252]
#   --padj 0.05  --lfc 1  --top 100
#
# Example
#   Rscript deseq2_parents.R --samples QC_PM_RNAseq/rnaseq_samples.tsv \
#       --annotation G4_gene_Rgene.tsv --outdir DESeq2_Er3
################################################################################
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "rnaseq_utils.R"))
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap) })

opt <- parse_opts(list(samples = NA, annotation = NA, outdir = "DESeq2_Er3", `rsem-dir` = "",
                       resistant = "G4", susceptible = "G9", labels = "G4=WL2140,G9=Terno",
                       times = "PM0hr,PM24hr,PM4dpi", region = "chr4:8414176-19685252",
                       padj = 0.05, lfc = 1, top = 100))
padj_cut <- as.numeric(opt$padj); lfc_cut <- as.numeric(opt$lfc); top_n <- as.integer(opt$top)
R <- opt$resistant; S <- opt$susceptible; lab <- label_map(opt$labels)
lR <- relabel(R, lab); lS <- relabel(S, lab)
reg <- parse_region(opt$region); reg_name <- sprintf("%s_%.1f-%.1fMb", reg$chr, reg$start / 1e6, reg$end / 1e6)
times <- split_list(opt$times)
outdir <- file.path(opt$outdir, "parents"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- read_samples(opt$samples, opt$`rsem-dir`) %>% filter(genotype %in% c(R, S))
anno <- read_annotation(opt$annotation) %>% mutate(in_region = in_region(chr, start, end, reg))
add_anno <- function(x) x %>% mutate(id0 = base_id(gene_id)) %>%
  left_join(anno %>% select(-gene_id) %>% mutate(id0 = base_id(anno$gene_id)), by = "id0") %>%
  select(-id0) %>% mutate(is_Rgene = coalesce(is_Rgene, FALSE), in_region = coalesce(in_region, FALSE))

#==============================================================================
# Part A  G4 vs G9 within each time point
#==============================================================================
sa <- samples %>% filter(time %in% times) %>%
  mutate(genotype = factor(genotype, levels = c(S, R)), time = factor(time, levels = times),
         condition = factor(paste(genotype, time, sep = "_"),
                            levels = c(paste(S, times, sep = "_"), paste(R, times, sep = "_"))))
rownames(sa) <- sa$sample
print(table(sa$genotype, sa$time))
dds <- DESeqDataSetFromTximport(import_rsem(sa), colData = sa, design = ~ condition)
dds <- filter_low(dds, 10, 4)
dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "dds_parents_by_time.rds"))
write.csv(counts(dds, normalized = TRUE), file.path(outdir, "normalized_counts_parents.csv"))

resA <- bind_rows(lapply(times, function(tt) {
  nm <- paste0(R, "_vs_", S, "_", tt)
  out <- format_result(results(dds, contrast = c("condition", paste(R, tt, sep = "_"), paste(S, tt, sep = "_")),
                               alpha = padj_cut), nm, padj_cut, lfc_cut) %>%
    add_anno() %>% mutate(time = tt, label = paste(lR, "vs", lS, "at", tt))
  write_csv(out, file.path(outdir, paste0(nm, "_ALL.csv")))
  write_csv(filter(out, significant), file.path(outdir, paste0(nm, "_DE.csv")))
  write_csv(filter(out, in_region), file.path(outdir, paste0(nm, "_", reg_name, "_ALL.csv")))
  write_csv(filter(out, in_region, is_Rgene), file.path(outdir, paste0(nm, "_", reg_name, "_Rgenes.csv")))
  out
}))
write_csv(resA, file.path(outdir, paste0("MASTER_", R, "_vs_", S, "_by_time_long.csv")))

summA <- resA %>% group_by(contrast, label) %>%
  summarise(tested = n(), DE = sum(significant), Up = sum(direction == "Up"), Down = sum(direction == "Down"),
            DE_Rgenes = sum(significant & is_Rgene), region_genes = sum(in_region),
            region_DE = sum(in_region & significant), region_Rgenes = sum(in_region & is_Rgene),
            region_DE_Rgenes = sum(in_region & significant & is_Rgene), .groups = "drop")
write_csv(summA, file.path(outdir, "SUMMARY_DE_counts_by_time.csv")); print(summA)

# plots: volcano, positional log2FC in the interval, heatmaps
vsd <- vst(dds, blind = FALSE); vmat <- assay(vsd)
ann_col <- data.frame(Genotype = relabel(as.character(sa$genotype), lab), Time = sa$time, row.names = sa$sample)
vol <- resA %>% mutate(y = -log10(pmax(padj, 1e-300)),
                       class = case_when(significant & in_region & is_Rgene ~ "DE, interval R gene",
                                         significant & in_region ~ "DE, interval",
                                         significant & is_Rgene ~ "DE, R gene",
                                         significant ~ "DE", TRUE ~ "Not DE"))
ggsave(file.path(outdir, "VOLCANO_by_time.pdf"),
       ggplot(vol, aes(log2FoldChange, y, colour = class)) + geom_point(size = 0.7, alpha = 0.7) +
         geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey55") +
         geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey55") +
         facet_wrap(~ label, nrow = 1) +
         scale_colour_manual(values = c("DE, interval R gene" = "red3", "DE, interval" = "orange2",
                                        "DE, R gene" = "purple3", "DE" = "grey25", "Not DE" = "grey82")) +
         theme_bw() + labs(x = paste0("log2FC ", lR, " vs ", lS), y = "-log10 adjusted P", colour = NULL),
       width = 12, height = 4.5)
ggsave(file.path(outdir, paste0("POSITION_log2FC_", reg_name, ".pdf")),
       ggplot(filter(vol, in_region), aes(start / 1e6, log2FoldChange,
                                          colour = ifelse(is_Rgene, "R gene", "Other"), alpha = significant)) +
         geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey55") +
         geom_point(size = 1.3) + facet_wrap(~ label, ncol = 1) +
         scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.3)) +
         scale_colour_manual(values = c("R gene" = "red3", "Other" = "grey30")) +
         theme_bw() + labs(x = paste(reg$chr, "position (Mb)"), y = paste0("log2FC ", lR, " vs ", lS),
                           colour = NULL, alpha = "DE"),
       width = 10, height = 7)
heat <- function(ids, file, title) {
  ids <- unique(ids[ids %in% rownames(vmat)]); if (length(ids) < 2) return(invisible(NULL))
  m <- vmat[ids, , drop = FALSE]; m <- m - rowMeans(m)
  info <- resA %>% distinct(gene_id, .keep_all = TRUE)
  rownames(m) <- with(info[match(ids, info$gene_id), ],
                      paste0(gene_id, ifelse(is_Rgene, paste0(" | ", Rgene_family), "")))
  pdf(file.path(outdir, file), width = 10, height = max(6, min(24, 0.18 * nrow(m) + 4)))
  pheatmap(m, annotation_col = ann_col, cluster_cols = FALSE, fontsize_row = 5, main = title)
  dev.off()
}
top_ids <- function(d) d %>% filter(significant) %>% arrange(padj) %>% distinct(gene_id) %>% slice_head(n = top_n) %>% pull(gene_id)
heat(top_ids(resA), "HEATMAP_top_DE_genome.pdf", paste("Top", top_n, "DE genes"))
heat(top_ids(filter(resA, in_region)), paste0("HEATMAP_top_DE_", reg_name, ".pdf"), paste("Top DE genes,", reg_name))
heat(top_ids(filter(resA, in_region, is_Rgene)), paste0("HEATMAP_DE_Rgenes_", reg_name, ".pdf"), paste("DE R genes,", reg_name))

#==============================================================================
# Part B  Infection-response interaction (0 h and 24 h)
#==============================================================================
sb <- samples %>% filter(time %in% c("PM0hr", "PM24hr")) %>%
  mutate(genotype = factor(genotype, levels = c(S, R)), time = factor(time, levels = c("PM0hr", "PM24hr")))
rownames(sb) <- sb$sample
ddi <- DESeqDataSetFromTximport(import_rsem(sb), colData = sb, design = ~ genotype * time)
ddi <- filter_low(ddi, 10, 2)
ddi <- DESeq(ddi)
coef <- grep(paste0("genotype", R, ".*timePM24hr|timePM24hr.*genotype", R), resultsNames(ddi), value = TRUE)
if (length(coef) != 1) stop("Interaction coefficient not found; resultsNames: ", paste(resultsNames(ddi), collapse = ", "))
resB <- format_result(results(ddi, name = coef, alpha = padj_cut),
                      paste0(R, "_response_vs_", S, "_response"), padj_cut, lfc_cut)
write_csv(resB, file.path(opt$outdir, paste0(R, "_response_vs_", S, "_response_ALL.csv")))
write_csv(filter(resB, significant), file.path(opt$outdir, paste0(R, "_response_vs_", S, "_response_DE.csv")))
saveRDS(ddi, file.path(outdir, "dds_parents_interaction.rds"))
message(sprintf("Interaction %s: %d DE (%d up, %d down)", coef, sum(resB$significant),
                sum(resB$direction == "Up"), sum(resB$direction == "Down")))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))