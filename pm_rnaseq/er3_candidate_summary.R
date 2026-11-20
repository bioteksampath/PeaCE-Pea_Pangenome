#!/usr/bin/env Rscript
################################################################################
# er3_candidate_summary.R
# Combine parent and introgression-line DESeq2 results into R-gene candidate
# tables for the Er3-associated intervals, with LFC heatmaps and an Excel workbook.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Powdery mildew expression analyses
# Figure     : Fig. 4f context; candidate tables (Supplementary)
#
# Uses the *_response_vs_*_response_ALL.csv files from deseq2_parents.R and deseq2_introgression_lines.R.
# A gene has "support" in a comparison if padj < 0.05 and |log2FC| >= 1.
#
# Options (defaults)
#   --deseq-dir   folder with *_ALL.csv from the deseq2_*.R scripts                 [DESeq2_Er3]
#   --annotation  GeneID Chr Start End RgeneFamily                         (required)
#   --outdir      [Er3_candidate_summary]
#   --parent      parent interaction prefix                                [G4_response_vs_G9_response]
#   --reference   susceptible IL                                           [IN9]
#   --groups      [Chr4=IN16,IN33,IN36;Chr6=IN24,IN40]
#   --other       lines summarised on their own                            [IN32]
#   --regions     interval per group            [Chr4=chr4:8414176-19685252;Chr6=chr6:7750000-10250000]
#   --labels      [G4=WL2140,G9=Terno]
#   --padj 0.05 --lfc 1 --top 50
#
# Example
#   Rscript er3_candidate_summary.R --deseq-dir DESeq2_Er3 --annotation G4_gene_Rgene.tsv
################################################################################
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "rnaseq_utils.R"))
suppressPackageStartupMessages({ library(ggplot2); library(pheatmap); library(openxlsx) })

opt <- parse_opts(list(`deseq-dir` = "DESeq2_Er3", annotation = NA, outdir = "Er3_candidate_summary",
                       parent = "G4_response_vs_G9_response", reference = "IN9",
                       groups = "Chr4=IN16,IN33,IN36;Chr6=IN24,IN40", other = "IN32",
                       regions = "Chr4=chr4:8414176-19685252;Chr6=chr6:7750000-10250000",
                       labels = "G4=WL2140,G9=Terno", padj = 0.05, lfc = 1, top = 50))
padj_cut <- as.numeric(opt$padj); lfc_cut <- as.numeric(opt$lfc); top_n <- as.integer(opt$top)
kv <- function(x) { p <- strsplit(strsplit(x, ";")[[1]], "="); setNames(lapply(p, function(z) z[2]), vapply(p, function(z) trimws(z[1]), "")) }
groups  <- lapply(kv(opt$groups), split_list)
regions <- lapply(kv(opt$regions), parse_region)
other   <- if (nzchar(opt$other)) split_list(opt$other) else character(0)
ref <- opt$reference; lab <- label_map(opt$labels)
plot_dir <- file.path(opt$outdir, "plots"); dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
is_de <- function(l, p) !is.na(l) & !is.na(p) & p < padj_cut & abs(l) >= lfc_cut

#---------------------------- comparisons and inputs ---------------------------
lines_all <- c(unlist(groups), other)
cmp <- tibble(key = c("Parent", paste0(names(groups), "_group"), lines_all),
              file = c(opt$parent, paste0(names(groups), "_response_vs_", ref, "_response"),
                       paste0(lines_all, "_response_vs_", ref, "_response")),
              label = c(paste0("Parent\n", paste(relabel(strsplit(sub("_response_vs_", ",", sub("_response$", "", opt$parent)), ",")[[1]], lab), collapse = "/")),
                        paste0(names(groups), "\nmean"), lines_all))
anno <- read_annotation(opt$annotation)
for (g in names(regions)) anno[[paste0("in_", g)]] <- in_region(anno$chr, anno$start, anno$end, regions[[g]])

master <- anno
for (i in seq_len(nrow(cmp))) {
  f <- file.path(opt$`deseq-dir`, paste0(cmp$file[i], "_ALL.csv"))
  if (!file.exists(f)) stop("Missing DESeq2 result: ", f)
  x <- read_csv(f, show_col_types = FALSE) %>% transmute(gene_id, LFC = log2FoldChange, padj)
  names(x)[2:3] <- paste0(cmp$key[i], c("_LFC", "_padj"))
  master <- full_join(master, x, by = "gene_id")
}
master <- master %>% mutate(is_Rgene = coalesce(is_Rgene, FALSE))
for (k in cmp$key) master[[paste0(k, "_support")]] <- is_de(master[[paste0(k, "_LFC")]], master[[paste0(k, "_padj")]])
for (g in names(groups))
  master[[paste0(g, "_line_support_n")]] <- rowSums(as.matrix(master[paste0(groups[[g]], "_support")]))
master$resistant_line_support_n <- rowSums(as.matrix(master[paste0(lines_all, "_support")]))
master$total_support <- master$resistant_line_support_n + master$Parent_support +
  rowSums(as.matrix(master[paste0(names(groups), "_group_support")]))

#---------------------------- candidate tables ---------------------------------
tabs <- list()
for (g in names(groups)) {
  if (!g %in% names(regions)) next
  tabs[[paste0(g, "_R_candidates")]] <- master %>%
    filter(is_Rgene, .data[[paste0("in_", g)]]) %>%
    arrange(desc(.data[[paste0(g, "_group_support")]]), desc(.data[[paste0(g, "_line_support_n")]]),
            .data[[paste0(g, "_group_padj")]])
}
for (l in other)
  tabs[[paste0(l, "_R_candidates")]] <- master %>% filter(is_Rgene, .data[[paste0(l, "_support")]]) %>%
    mutate(score = abs(.data[[paste0(l, "_LFC")]]) * -log10(pmax(.data[[paste0(l, "_padj")]], 1e-300))) %>%
    arrange(desc(Parent_support), desc(score))
tabs$Genome_R_candidates <- master %>% filter(is_Rgene, total_support > 0) %>%
  arrange(desc(total_support), desc(resistant_line_support_n))
tabs$Multi_line_R <- master %>% filter(is_Rgene, resistant_line_support_n >= 2) %>%
  arrange(desc(resistant_line_support_n), desc(total_support))
tabs$Parent_IL_overlap <- master %>%
  filter(Parent_support, if_any(all_of(c(paste0(names(groups), "_group_support"), paste0(other, "_support"))))) %>%
  arrange(desc(is_Rgene), desc(total_support))

#---------------------------- summaries ----------------------------------------
DE_summary <- bind_rows(lapply(seq_len(nrow(cmp)), function(i) {
  k <- cmp$key[i]; l <- master[[paste0(k, "_LFC")]]; d <- master[[paste0(k, "_support")]]
  row <- tibble(Comparison = gsub("\n", " ", cmp$label[i]), Genes_tested = sum(!is.na(l)), DE_genes = sum(d),
                Up = sum(d & l > 0, na.rm = TRUE), Down = sum(d & l < 0, na.rm = TRUE), DE_R_genes = sum(d & master$is_Rgene))
  for (g in names(regions)) {
    row[[paste0(g, "_DE")]] <- sum(d & master[[paste0("in_", g)]], na.rm = TRUE)
    row[[paste0(g, "_DE_R")]] <- sum(d & master[[paste0("in_", g)]] & master$is_Rgene, na.rm = TRUE)
  }
  row
}))
overlap <- tibble(Category = "Annotated R genes", Count = sum(master$is_Rgene))
for (g in intersect(names(groups), names(regions))) {
  t <- tabs[[paste0(g, "_R_candidates")]]
  overlap <- bind_rows(overlap, tibble(
    Category = c(paste(g, "interval R genes"), paste(g, "group-supported R genes"),
                 paste(g, "R genes supported by >= 2 lines"), paste("Parent +", g, "group supported R genes")),
    Count = c(nrow(t), sum(t[[paste0(g, "_group_support")]]), sum(t[[paste0(g, "_line_support_n")]] >= 2),
              sum(t$Parent_support & t[[paste0(g, "_group_support")]]))))
}
overlap <- bind_rows(overlap, tibble(Category = c("R genes with any DE support", "R genes supported by >= 2 resistant lines"),
                                     Count = c(nrow(tabs$Genome_R_candidates), nrow(tabs$Multi_line_R))))
print(DE_summary, n = Inf); print(overlap, n = Inf)

#---------------------------- plots ---------------------------------------------
ggsave(file.path(plot_dir, "01_DE_counts.pdf"),
       ggplot(pivot_longer(DE_summary, c(Up, Down), names_to = "Direction", values_to = "Genes"),
              aes(reorder(Comparison, Genes), Genes, fill = Direction)) + geom_col() + coord_flip() +
         theme_bw(base_size = 12) + labs(x = NULL, y = "DE genes", subtitle = sprintf("padj < %g; |log2FC| >= %g", padj_cut, lfc_cut)),
       width = 8, height = 5)
lfc_heat <- function(d, file, title, max_genes = top_n) {
  d <- distinct(d, gene_id, .keep_all = TRUE) %>% slice_head(n = max_genes)
  m <- as.matrix(d[paste0(cmp$key, "_LFC")]); rownames(m) <- d$gene_id; colnames(m) <- cmp$label
  m <- m[rowSums(!is.na(m)) > 0, , drop = FALSE]; if (!nrow(m)) return(invisible(NULL))
  m[m > 5] <- 5; m[m < -5] <- -5                      # capped for display only
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") pdf(file.path(plot_dir, paste0(file, ".pdf")), width = 10, height = max(5, 2.5 + 0.22 * nrow(m)))
    else png(file.path(plot_dir, paste0(file, ".png")), width = 3000, height = max(1600, 80 * nrow(m) + 700), res = 300)
    pheatmap(m, cluster_rows = nrow(m) >= 2, cluster_cols = FALSE, border_color = NA, na_col = "grey90",
             fontsize_row = ifelse(nrow(m) <= 30, 8, 6), fontsize_col = 8, angle_col = 0,
             main = paste0(title, "\nDESeq2 log2 fold change"))
    dev.off()
  }
}
lfc_heat(tabs$Genome_R_candidates, "02_priority_Rgene_LFC_heatmap", "R genes with DE support")
for (g in intersect(names(groups), names(regions)))
  lfc_heat(tabs[[paste0(g, "_R_candidates")]], paste0("03_", g, "_interval_Rgene_LFC_heatmap"),
           paste(g, "interval R genes"), max_genes = 100)
for (l in other) lfc_heat(tabs[[paste0(l, "_R_candidates")]], paste0("04_", l, "_Rgene_LFC_heatmap"), paste(l, "DE R genes"), 30)

#---------------------------- outputs -------------------------------------------
write_csv(master, file.path(opt$outdir, "Er3_candidate_master.csv"))
write_csv(DE_summary, file.path(opt$outdir, "DE_summary.csv"))
write_csv(overlap, file.path(opt$outdir, "candidate_overlap_summary.csv"))
for (n in names(tabs)) write_csv(tabs[[n]], file.path(opt$outdir, paste0(n, ".csv")))

readme <- tibble(Item = c("DE threshold", names(regions), names(groups), "Reference line", "Parent test", "Note"),
                 Value = c(sprintf("padj < %g and |log2FC| >= %g", padj_cut, lfc_cut),
                           vapply(regions, function(r) sprintf("%s:%s-%s", r$chr, format(r$start, big.mark = ","), format(r$end, big.mark = ",")), ""),
                           vapply(groups, paste, "", collapse = ", "), ref,
                           "(resistant 24h - 0h) - (susceptible 24h - 0h)",
                           "Interaction tests compare infection responses with the susceptible line; expression support does not establish causality."))
wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", border = "bottom")
sheets <- c(list(README = readme, DE_summary = DE_summary, Overlap_summary = overlap), tabs, list(Candidate_master = master))
for (n in names(sheets)) {
  sn <- substr(n, 1, 31); addWorksheet(wb, sn)
  writeData(wb, sn, sheets[[n]], headerStyle = hs, withFilter = TRUE); freezePane(wb, sn, firstRow = TRUE)
}
pngs <- list.files(plot_dir, pattern = "\\.png$", full.names = TRUE)
if (length(pngs)) { addWorksheet(wb, "Plots")
  for (i in seq_along(pngs)) insertImage(wb, "Plots", pngs[i], startRow = 1 + 32 * ((i - 1) %/% 2),
                                         startCol = 1 + 11 * ((i - 1) %% 2), width = 8, height = 8) }
saveWorkbook(wb, file.path(opt$outdir, "Er3_DESeq2_summary.xlsx"), overwrite = TRUE)
message("Done: ", opt$outdir)