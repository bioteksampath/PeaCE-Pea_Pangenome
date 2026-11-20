#!/usr/bin/env Rscript
################################################################################
# deseq2_introgression_lines.R
# Infection response of P. fulvum introgression lines (ILs) compared with the
# susceptible IL (IN9), 0 h and 24 h after powdery mildew inoculation.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Powdery mildew expression analyses (introgression lines)
# Figure     : Fig. 4f context; Er3 candidate support (er3_candidate_summary.R)
#
# Part A  Each line, 24 h vs 0 h                      design ~ genotype_time
# Part B  Each line's response vs IN9 response        design ~ genotype * time
#           coefficient genotypeINx.timePM24hr = (INx 24h - INx 0h) - (IN9 24h - IN9 0h)
# Part C  Group means vs IN9 (numeric contrasts)      e.g. Chr4 = 1/3 (IN16 + IN33 + IN36)
# Part D  Group 1 vs group 2                          e.g. Chr4 mean - Chr6 mean
#
# Genes kept if >= 10 reads in >= 2 samples. DE: padj < 0.05 and |log2FC| >= 1.
#
# Options (defaults)
#   --samples    sample sheet from rnaseq_qc.R                                  (required)
#   --outdir     [DESeq2_Er3]
#   --reference  susceptible reference line                                       [IN9]
#   --groups     line groups, "name=line,line;name=line"     [Chr4=IN16,IN33,IN36;Chr6=IN24,IN40]
#   --other      other lines tested individually only                             [IN32]
#   --padj 0.05  --lfc 1
#
# Example
#   Rscript deseq2_introgression_lines.R --samples QC_PM_RNAseq/rnaseq_samples.tsv --outdir DESeq2_Er3
################################################################################
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "rnaseq_utils.R"))
suppressPackageStartupMessages(library(DESeq2))

opt <- parse_opts(list(samples = NA, outdir = "DESeq2_Er3", `rsem-dir` = "", reference = "IN9",
                       groups = "Chr4=IN16,IN33,IN36;Chr6=IN24,IN40", other = "IN32",
                       padj = 0.05, lfc = 1))
padj_cut <- as.numeric(opt$padj); lfc_cut <- as.numeric(opt$lfc)
ref <- opt$reference
groups <- lapply(strsplit(strsplit(opt$groups, ";")[[1]], "="), function(g) split_list(g[2]))
names(groups) <- vapply(strsplit(strsplit(opt$groups, ";")[[1]], "="), function(g) trimws(g[1]), "")
other <- if (nzchar(opt$other)) split_list(opt$other) else character(0)
lines <- unique(c(ref, unlist(groups), other))
outdir <- opt$outdir; dir.create(file.path(outdir, "introgression_lines"), showWarnings = FALSE, recursive = TRUE)
od <- file.path(outdir, "introgression_lines")

samples <- read_samples(opt$samples, opt$`rsem-dir`) %>%
  filter(genotype %in% lines, time %in% c("PM0hr", "PM24hr")) %>%
  mutate(genotype = factor(genotype, levels = lines), time = factor(time, levels = c("PM0hr", "PM24hr")))
rownames(samples) <- samples$sample
miss <- setdiff(lines, unique(as.character(samples$genotype)))
if (length(miss)) stop("No samples for: ", paste(miss, collapse = ", "))
print(table(samples$genotype, samples$time))
write_tsv(samples, file.path(od, "IL_samples_used.tsv"))
txi <- import_rsem(samples)

save_res <- function(res, name, dir = od) {
  out <- format_result(res, name, padj_cut, lfc_cut)
  write_csv(out, file.path(dir, paste0(name, "_ALL.csv")))
  write_csv(filter(out, significant), file.path(dir, paste0(name, "_DE.csv")))
  out
}

#==============================================================================
# Part A  each line: 24 h vs 0 h
#==============================================================================
samples$condition2 <- factor(paste(samples$genotype, samples$time, sep = "_"))
dds_c <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition2)
dds_c <- DESeq(filter_low(dds_c, 10, 2))
saveRDS(dds_c, file.path(od, "dds_IL_condition_model.rds"))
resA <- bind_rows(lapply(lines, function(g)
  save_res(results(dds_c, contrast = c("condition2", paste0(g, "_PM24hr"), paste0(g, "_PM0hr")), alpha = padj_cut),
           paste0(g, "_PM24hr_vs_PM0hr"))))

#==============================================================================
# Part B  each line's response vs the reference response (interaction)
#==============================================================================
dds_i <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ genotype + time + genotype:time)
mm <- model.matrix(design(dds_i), as.data.frame(colData(dds_i)))
if (qr(mm)$rank != ncol(mm)) stop("Interaction design is not full rank.")
dds_i <- DESeq(filter_low(dds_i, 10, 2))
saveRDS(dds_i, file.path(od, "dds_IL_genotype_time_interaction.rds"))
rn <- resultsNames(dds_i); writeLines(rn, file.path(od, "DESeq2_interaction_resultsNames.txt"))
coef_of <- function(g) {
  h <- grep(paste0("genotype", g, "\\.timePM24hr|timePM24hr\\.genotype", g), rn, value = TRUE)
  if (length(h) != 1) stop("Interaction coefficient for ", g, " not unique: ", paste(h, collapse = ", "))
  h
}
test_lines <- setdiff(lines, ref)
coefs <- setNames(vapply(test_lines, coef_of, ""), test_lines)
resB <- bind_rows(lapply(test_lines, function(g)
  save_res(results(dds_i, name = coefs[[g]], alpha = padj_cut),
           paste0(g, "_response_vs_", ref, "_response"), outdir)))

#==============================================================================
# Part C/D  group means vs reference; group 1 vs group 2
#==============================================================================
weights_vec <- function(w) { v <- setNames(rep(0, length(rn)), rn); v[names(w)] <- w; v }
group_w <- lapply(groups, function(g) setNames(rep(1 / length(g), length(g)), coefs[g]))
resC <- bind_rows(lapply(names(groups), function(nm) {
  message("Contrast ", nm, " mean: ", paste(sprintf("%.3f*%s", group_w[[nm]], names(group_w[[nm]])), collapse = " + "))
  save_res(results(dds_i, contrast = weights_vec(group_w[[nm]]), alpha = padj_cut),
           paste0(nm, "_response_vs_", ref, "_response"), outdir)
}))
if (length(groups) == 2) {
  g1 <- names(groups)[1]; g2 <- names(groups)[2]
  w <- c(group_w[[g1]], -group_w[[g2]])
  resC <- bind_rows(resC, save_res(results(dds_i, contrast = weights_vec(w), alpha = padj_cut),
                                   paste0(g1, "_response_vs_", g2, "_response"), outdir))
}

#==============================================================================
# Summaries and expression matrices
#==============================================================================
allres <- bind_rows(resA, resB, resC)
summ <- allres %>% group_by(contrast) %>%
  summarise(tested = n(), DE = sum(significant), Up = sum(direction == "Up"), Down = sum(direction == "Down"),
            .groups = "drop")
write_csv(summ, file.path(od, "SUMMARY_DE_counts.csv")); print(summ, n = Inf)
write_csv(bind_rows(resB, resC), file.path(od, "MASTER_interaction_results_long.csv"))
write_csv(resA %>% select(gene_id, contrast, log2FoldChange, padj) %>%
            pivot_wider(names_from = contrast, values_from = c(log2FoldChange, padj)),
          file.path(od, "MASTER_individual_response_wide.csv"))
write.csv(assay(vst(dds_c, blind = FALSE)), file.path(od, "VST_IL_expression_matrix.csv"))
write.csv(counts(dds_c, normalized = TRUE), file.path(od, "normalized_counts_IL.csv"))
capture.output(sessionInfo(), file = file.path(od, "sessionInfo.txt"))
message("Done: ", outdir)