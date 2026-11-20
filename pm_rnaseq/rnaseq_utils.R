################################################################################
# rnaseq_utils.R
# Shared helpers for the powdery mildew RNA-seq scripts.
# Sourced by each script; not run on its own.
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
})

# --name value command-line options; unknown names stop with an error
parse_opts <- function(defaults, usage = "") {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0 || (length(a) && !all(startsWith(a[c(TRUE, FALSE)], "--"))))
    stop("Options must be --name value pairs.\n", usage, call. = FALSE)
  for (i in seq(1, length(a), by = 2)) {
    k <- sub("^--", "", a[i])
    if (!k %in% names(defaults)) stop("Unknown option --", k, "\n", usage, call. = FALSE)
    defaults[[k]] <- a[i + 1]
  }
  miss <- names(defaults)[vapply(defaults, function(x) identical(x, NA), logical(1))]
  if (length(miss)) stop("Required: ", paste0("--", miss, collapse = ", "), "\n", usage, call. = FALSE)
  defaults
}

split_list <- function(x) trimws(strsplit(x, ",")[[1]])

# "G4=WL2140,G9=Terno" -> named vector used for plot labels
label_map <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  kv <- strsplit(split_list(x), "=")
  setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
}
relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)

# "chr4:8414176-19685252" -> list(chr = "chr4", start, end)
parse_region <- function(x) {
  m <- regmatches(x, regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", x))[[1]]
  if (!length(m)) stop("Region must look like chr4:8414176-19685252, got: ", x)
  list(chr = tolower(m[2]), start = as.numeric(gsub(",", "", m[3])), end = as.numeric(gsub(",", "", m[4])))
}

# "Pful.WL2140.Chr4" / "Chr4" / "4" -> "chr4"
clean_chr <- function(x) {
  x <- sub(".*[._]", "", as.character(x))
  x <- tolower(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

# Sample sheet with at least: sample genotype time replicate rsem_file
read_samples <- function(file, rsem_dir = NULL) {
  s <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("sample", "genotype", "time", "replicate", "rsem_file")
  miss <- setdiff(need, colnames(s))
  if (length(miss)) stop("Sample sheet is missing columns: ", paste(miss, collapse = ", "))
  if (!is.null(rsem_dir) && nzchar(rsem_dir))
    s$rsem_file <- ifelse(file.exists(s$rsem_file), s$rsem_file, file.path(rsem_dir, basename(s$rsem_file)))
  bad <- s$rsem_file[!file.exists(s$rsem_file)]
  if (length(bad)) stop("Missing RSEM files:\n", paste(head(bad, 10), collapse = "\n"))
  rownames(s) <- s$sample
  s
}

# tximport of RSEM gene results; genes with zero/NA effective length in any sample removed
import_rsem <- function(samples) {
  suppressPackageStartupMessages(library(tximport))
  f <- setNames(samples$rsem_file, samples$sample)
  txi <- tximport(f, type = "rsem", txIn = FALSE, txOut = FALSE)
  ok <- rowSums(!is.na(txi$length) & txi$length > 0) == ncol(txi$length)
  message("Genes removed (invalid effective length): ", sum(!ok))
  for (k in c("abundance", "counts", "length")) txi[[k]] <- txi[[k]][ok, , drop = FALSE]
  txi
}

# keep genes with >= min_count reads in >= min_samples samples
filter_low <- function(dds, min_count, min_samples) {
  keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
  message(sprintf("Genes kept: %d of %d (>= %d reads in >= %d samples)",
                  sum(keep), length(keep), min_count, min_samples))
  dds[keep, ]
}

format_result <- function(res, contrast_name, padj_cutoff, lfc_cutoff) {
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(contrast = contrast_name,
           significant = !is.na(padj) & padj < padj_cutoff &
                         !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_cutoff,
           direction = case_when(significant & log2FoldChange > 0 ~ "Up",
                                 significant & log2FoldChange < 0 ~ "Down",
                                 TRUE ~ "Not_DE")) %>%
    arrange(padj, desc(abs(log2FoldChange)))
}

# Gene annotation: GeneID Chr Start End RgeneFamily (tab-separated, header)
read_annotation <- function(file) {
  a <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("GeneID", "Chr", "Start", "End", "RgeneFamily")
  miss <- setdiff(need, colnames(a))
  if (length(miss)) stop("Annotation is missing columns: ", paste(miss, collapse = ", "))
  a %>% transmute(gene_id = trimws(GeneID), chr = clean_chr(Chr),
                  start = as.numeric(Start), end = as.numeric(End),
                  Rgene_family = trimws(as.character(RgeneFamily))) %>%
    mutate(Rgene_family = ifelse(is.na(Rgene_family) | toupper(Rgene_family) == "NA", "", Rgene_family),
           is_Rgene = Rgene_family != "") %>%
    distinct(gene_id, .keep_all = TRUE)
}

in_region <- function(chr, start, end, r) chr == r$chr & end >= r$start & start <= r$end

# DESeq2 gene IDs may carry RSEM suffixes; match on the base ID
base_id <- function(x) sub("\\.[0-9]+$", "", sub("_.*$", "", x))
