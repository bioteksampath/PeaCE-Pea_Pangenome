#!/usr/bin/env Rscript
################################################################################
# rnaseq_qc.R
# Build the sample sheet from RSEM file names and run RNA-seq QC
# (library size, expressed genes, PCA, sample correlation/distance, replicates).
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : RNA isolation and analysis of introgression parents and population
#
# Sample names must look like:  G4-PM0hr-1, G9-PM24hr-3, CSPGR-IN16-PM24hr-2
#   genotype = G4 / G9 / IN<number>;  time = PM0hr / PM24hr / PM4dpi;  replicate = last number
#
# Options (defaults)
#   --rsem-dir     directory with <sample>.genes.results                      (required)
#   --outdir       output directory                                           [QC_PM_RNAseq]
#   --pattern      regex of sample names to keep                              [^(G4-PM|G9-PM|CSPGR-IN)]
#   --resistant    resistant parent                                           [G4]
#   --susceptible  susceptible parent                                         [G9]
#   --groups       IL groups, "group=line,line;group=line"
#                  [Chr4=IN16,IN33,IN36;Chr6=IN24,IN40;unknown=IN32;noEr3=IN9]
#   --susceptible-lines  ILs scored susceptible                               [IN9]
#   --labels       display names for plots                                    [G4=WL2140,G9=Terno]
#   --min-count 10 --min-samples 2   permissive QC gene filter
#
# Outputs (in --outdir)
#   rnaseq_samples.tsv  <- the sample sheet used by the deseq2_*.R scripts
#   QC_library_size_expressed_genes.tsv, QC_library_size.pdf
#   PCA_*.pdf (+ coordinates), QC_sample_correlation/distance heatmaps,
#   QC_replicate_correlations.tsv, QC_pairwise_replicate_correlations.tsv,
#   normalized_counts_QC.csv, VST_expression_matrix.csv
#
# Example
#   Rscript rnaseq_qc.R --rsem-dir pm_rsem/rsem --outdir QC_PM_RNAseq
################################################################################
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(script_dir, "rnaseq_utils.R"))
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel) })

opt <- parse_opts(list(`rsem-dir` = NA, outdir = "QC_PM_RNAseq", pattern = "^(G4-PM|G9-PM|CSPGR-IN)",
                       resistant = "G4", susceptible = "G9",
                       groups = "Chr4=IN16,IN33,IN36;Chr6=IN24,IN40;unknown=IN32;noEr3=IN9",
                       `susceptible-lines` = "IN9", labels = "G4=WL2140,G9=Terno",
                       `min-count` = 10, `min-samples` = 2))
outdir <- opt$outdir; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
labs_map <- label_map(opt$labels)

#---------------------------- 1. sample sheet from file names -------------------
files <- list.files(opt$`rsem-dir`, pattern = "\\.genes\\.results$", full.names = TRUE)
names(files) <- sub("\\.genes\\.results$", "", basename(files))
files <- files[grepl(opt$pattern, names(files))]
if (!length(files)) stop("No *.genes.results matching ", opt$pattern, " in ", opt$`rsem-dir`)

grp <- do.call(rbind, lapply(strsplit(strsplit(opt$groups, ";")[[1]], "="), function(g)
  data.frame(genotype = split_list(g[2]), er3_group = trimws(g[1]))))
susc_lines <- split_list(opt$`susceptible-lines`)

samples <- tibble(sample = names(files), rsem_file = unname(files)) %>%
  mutate(genotype = case_when(grepl(paste0("^", opt$resistant, "-"), sample) ~ opt$resistant,
                              grepl(paste0("^", opt$susceptible, "-"), sample) ~ opt$susceptible,
                              grepl("IN[0-9]+", sample) ~ sub(".*(IN[0-9]+).*", "\\1", sample),
                              TRUE ~ NA_character_),
         time = sub(".*(PM0hr|PM24hr|PM4dpi).*", "\\1", sample),
         time = ifelse(time == sample, NA, time),
         replicate = suppressWarnings(as.integer(sub(".*-([0-9]+)$", "\\1", sample))),
         experiment = ifelse(genotype %in% c(opt$resistant, opt$susceptible), "parents", "introgression_lines")) %>%
  left_join(grp, by = "genotype") %>%
  mutate(er3_group = case_when(genotype == opt$resistant ~ "donor",
                               genotype == opt$susceptible ~ "recurrent",
                               TRUE ~ er3_group),
         phenotype = ifelse(genotype == opt$susceptible | genotype %in% susc_lines, "S", "R"),
         condition = paste(genotype, time, sep = "_"),
         label = relabel(genotype, labs_map)) %>%
  select(sample, genotype, label, time, replicate, phenotype, er3_group, experiment, condition, rsem_file)
if (anyNA(samples$genotype) || anyNA(samples$time))
  stop("Could not parse genotype/time for: ",
       paste(samples$sample[is.na(samples$genotype) | is.na(samples$time)], collapse = ", "))
write_tsv(samples, file.path(outdir, "rnaseq_samples.tsv"))
print(table(samples$genotype, samples$time))

#---------------------------- 2. import and normalise ---------------------------
samples <- as.data.frame(samples); rownames(samples) <- samples$sample
txi <- import_rsem(samples)
cd <- samples %>% mutate(time = factor(time, levels = c("PM0hr", "PM24hr", "PM4dpi")))
rownames(cd) <- cd$sample
dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ 1)
dds <- filter_low(dds, as.integer(opt$`min-count`), as.integer(opt$`min-samples`))
dds <- estimateSizeFactors(dds)
raw <- counts(dds); norm <- counts(dds, normalized = TRUE)
write.csv(norm, file.path(outdir, "normalized_counts_QC.csv"))

#---------------------------- 3. library size -----------------------------------
lib <- tibble(sample = colnames(raw), total_expected_counts = colSums(raw),
              genes_ge1 = colSums(raw >= 1), genes_ge10 = colSums(raw >= 10)) %>%
  left_join(samples %>% select(sample, label, time, replicate), by = "sample") %>% arrange(total_expected_counts)
write_tsv(lib, file.path(outdir, "QC_library_size_expressed_genes.tsv"))
ggsave(file.path(outdir, "QC_library_size.pdf"),
       ggplot(lib, aes(reorder(sample, total_expected_counts), total_expected_counts / 1e6, fill = label)) +
         geom_col() + coord_flip() + theme_bw() +
         labs(x = NULL, y = "RSEM expected counts (millions)", fill = NULL),
       width = 8, height = max(5, 0.15 * nrow(lib) + 2))

#---------------------------- 4. PCA ---------------------------------------------
vsd <- vst(dds, blind = TRUE); vmat <- assay(vsd)
write.csv(vmat, file.path(outdir, "VST_expression_matrix.csv"))
make_pca <- function(keep, file, title) {
  keep <- intersect(keep, colnames(vmat)); if (length(keep) < 3) return(invisible(NULL))
  p <- prcomp(t(vmat[, keep]), center = TRUE); pv <- 100 * p$sdev^2 / sum(p$sdev^2)
  d <- as.data.frame(p$x[, 1:min(4, ncol(p$x))]) %>% rownames_to_column("sample") %>% left_join(samples, by = "sample")
  write_tsv(d, file.path(outdir, sub("\\.pdf$", "_coordinates.tsv", file)))
  ggsave(file.path(outdir, file),
         ggplot(d, aes(PC1, PC2, colour = label, shape = time)) + geom_point(size = 3.5) +
           geom_text_repel(aes(label = replicate), size = 2.5, max.overlaps = Inf, show.legend = FALSE) +
           theme_bw() + labs(title = title, colour = "Genotype", shape = "Time",
                             x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])),
         width = 8, height = 6)
}
make_pca(samples$sample, "PCA_all_samples.pdf", "All samples")
for (tt in unique(samples$time)) make_pca(samples$sample[samples$time == tt], paste0("PCA_", tt, ".pdf"), tt)
make_pca(samples$sample[samples$experiment == "parents"], "PCA_parents.pdf", "Parents")
make_pca(samples$sample[samples$experiment == "introgression_lines"], "PCA_introgression_lines.pdf", "Introgression lines")

#---------------------------- 5. correlation, distance, replicates -------------
ann <- samples[, c("label", "time", "phenotype", "er3_group")]
sc <- cor(vmat); write.csv(sc, file.path(outdir, "QC_sample_correlations.csv"))
sd_mat <- as.matrix(dist(t(vmat)))
for (x in list(list(sc, "QC_sample_correlation_heatmap.pdf", "Sample correlation (VST)"),
               list(sd_mat, "QC_sample_distance_heatmap.pdf", "Sample distance (VST)"))) {
  pdf(file.path(outdir, x[[2]]), width = 13, height = 13)
  pheatmap(x[[1]], annotation_col = ann, annotation_row = ann, border_color = NA,
           fontsize_row = 5, fontsize_col = 5, main = x[[3]])
  dev.off()
}
pairs <- bind_rows(lapply(split(samples$sample, paste(samples$genotype, samples$time)), function(ss) {
  if (length(ss) < 2) return(NULL)
  cb <- combn(ss, 2)
  tibble(sample1 = cb[1, ], sample2 = cb[2, ])
})) %>%
  left_join(samples %>% select(sample, genotype, time), by = c("sample1" = "sample")) %>%
  mutate(correlation = sc[cbind(sample1, sample2)]) %>% arrange(correlation)
write_tsv(pairs, file.path(outdir, "QC_pairwise_replicate_correlations.tsv"))
write_tsv(pairs %>% group_by(genotype, time) %>%
            summarise(n_pairs = n(), mean_r = mean(correlation), min_r = min(correlation), .groups = "drop"),
          file.path(outdir, "QC_replicate_correlations.tsv"))
message("QC done: ", ncol(dds), " samples, ", nrow(dds), " genes. Sample sheet: ",
        file.path(outdir, "rnaseq_samples.tsv"))