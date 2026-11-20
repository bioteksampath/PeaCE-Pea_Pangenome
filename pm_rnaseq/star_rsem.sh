#!/usr/bin/env bash
################################################################################
# star_rsem.sh
# Powdery mildew RNA-seq: read trimming, STAR alignment and RSEM gene-level
# quantification against the P. fulvum WL2140 (G4) genome.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Powdery mildew expression analyses
# Figure     : Fig. 4f (input); DESeq2 in deseq2_*.R
#
# Flow
#   step 0  Build RSEM/STAR reference        rsem-prepare-reference --star
#   step 1  QC + trimming                    FastQC, cutadapt
#   step 2  STAR + RSEM per sample           rsem-calculate-expression --star
#   step 3  Gene matrices                    expected counts and TPM for all samples
#
# Usage
#   bash star_rsem.sh -g WL2140.fa -a WL2140.gff3 -s rnaseq_samples.tsv -o rsem_out [options]
#
#   -g  genome FASTA (P. fulvum WL2140)                          (required)
#   -a  gene annotation GTF/GFF3                                 (required for step 0)
#   -s  sample sheet: sample<TAB>R1.fq.gz<TAB>R2.fq.gz           (required)
#   -o  output directory                                         (required)
#   -t  threads                                   [16]
#   -p  steps to run, comma separated             [0,1,2,3]
#   -A  R1 adapter                                [AGATCGGAAGAGCACACGTCTGAACTCCAGTCA]
#   -B  R2 adapter                                [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
#
# Example
#   bash star_rsem.sh -g Pful.WL2140.fa -a Pful.WL2140.gff3 -s pm_rnaseq.tsv -o pm_rsem -t 32
#
# Software (Methods): FastQC v0.12.1, cutadapt v5.1, STAR v2.7.11b, RSEM v1.3.3
################################################################################
set -euo pipefail

THREADS=16
STEPS="0,1,2,3"
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"   # TruSeq (Novogene libraries)
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
TRIM_QUAL=20
MIN_LEN=30
# STAR options used in the production run: no chimeric output, end-to-end alignment (no soft clipping)
STAR_OPTS="--chimSegmentMin 100000 --alignEndsType EndToEnd"

usage() { sed -n '2,32p' "$0"; exit 1; }
while getopts "g:a:s:o:t:p:A:B:h" opt; do
  case $opt in
    g) GENOME=$OPTARG ;; a) ANNOT=$OPTARG ;; s) SAMPLES=$OPTARG ;; o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;; p) STEPS=$OPTARG ;; A) ADAPTER_R1=$OPTARG ;; B) ADAPTER_R2=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${GENOME:-}" && -n "${SAMPLES:-}" && -n "${OUT:-}" ]] || usage

log()      { echo "[$(date '+%F %T')] $*"; }
run_step() { [[ ",$STEPS," == *",$1,"* ]]; }

mkdir -p "$OUT"/{ref,fastqc,trimmed,rsem,matrix,logs}
exec > >(tee -a "$OUT/logs/star_rsem.log") 2>&1
grep -v -e '^#' -e '^[[:space:]]*$' -e '^sample[[:space:]]' "$SAMPLES" > "$OUT/samples.tsv"
REF="$OUT/ref/rsem_ref"
log "genome=$GENOME samples=$(wc -l < "$OUT/samples.tsv") out=$OUT"

#==============================================================================
# step 0  RSEM reference with STAR index
#==============================================================================
if run_step 0; then
  [[ -n "${ANNOT:-}" ]] || { echo "step 0 needs -a annotation (GTF/GFF3)"; exit 1; }
  log "step 0: rsem-prepare-reference --star"
  case "$ANNOT" in *.gtf) AOPT=(--gtf "$ANNOT") ;; *) AOPT=(--gff3 "$ANNOT") ;; esac
  [[ -s "$REF.grp" ]] || rsem-prepare-reference "${AOPT[@]}" --star -p "$THREADS" "$GENOME" "$REF" \
                           > "$OUT/logs/rsem_prepare_reference.log" 2>&1
fi

#==============================================================================
# step 1  FastQC and adapter/quality trimming
#==============================================================================
if run_step 1; then
  log "step 1: FastQC + cutadapt"
  while IFS=$'\t' read -r S R1 R2; do
    T1="$OUT/trimmed/${S}_R1.trim.fq.gz"; T2="$OUT/trimmed/${S}_R2.trim.fq.gz"
    [[ -s "$T2" ]] && continue
    fastqc -t 2 -o "$OUT/fastqc" "$R1" "$R2" > /dev/null 2>&1
    cutadapt -j "$THREADS" -a "$ADAPTER_R1" -A "$ADAPTER_R2" -q "$TRIM_QUAL" -m "$MIN_LEN" \
             -o "$T1" -p "$T2" "$R1" "$R2" > "$OUT/logs/$S.cutadapt.log"
  done < "$OUT/samples.tsv"
fi

#==============================================================================
# step 2  STAR alignment + RSEM quantification
#==============================================================================
if run_step 2; then
  log "step 2: RSEM (STAR)"
  while IFS=$'\t' read -r S _ _; do
    [[ -s "$OUT/rsem/$S.genes.results" ]] && continue
    log "  $S"
    rsem-calculate-expression --star --star-gzipped-read-file --paired-end \
        --star-options "$STAR_OPTS" --output-genome-bam -p "$THREADS" \
        "$OUT/trimmed/${S}_R1.trim.fq.gz" "$OUT/trimmed/${S}_R2.trim.fq.gz" \
        "$REF" "$OUT/rsem/$S" > "$OUT/logs/$S.rsem.log" 2>&1
  done < "$OUT/samples.tsv"
fi

#==============================================================================
# step 3  Gene x sample matrices (expected counts, TPM) and mapping rates
#==============================================================================
if run_step 3; then
  log "step 3: gene matrices"
  # all genes.results files list genes in the same order (same RSEM reference)
  for col in expected_count TPM; do
    c=$([[ $col == expected_count ]] && echo 5 || echo 6)
    awk -F'\t' -v c="$c" -v dir="$OUT/rsem" 'BEGIN{OFS="\t"}
        { s[NR]=$1 } END{
          printf "gene_id"; for (i=1; i<=NR; i++) printf "\t%s", s[i]; printf "\n"
          n=0; while ((getline line < (dir"/"s[1]".genes.results")) > 0) { n++; split(line, f, "\t"); id[n]=f[1] }
          close(dir"/"s[1]".genes.results")
          for (i=1; i<=NR; i++) { r=0; f2=dir"/"s[i]".genes.results"
            while ((getline line < f2) > 0) { r++; split(line, f, "\t"); v[r, i]=f[c] } close(f2) }
          for (r=2; r<=n; r++) { printf "%s", id[r]; for (i=1; i<=NR; i++) printf "\t%s", v[r, i]; printf "\n" } }' \
      "$OUT/samples.tsv" > "$OUT/matrix/genes.$col.tsv"
  done
  # STAR unique mapping rate per sample
  echo -e "sample\tinput_reads\tuniquely_mapped_pct\tmulti_mapped_pct" > "$OUT/matrix/star_mapping.tsv"
  while IFS=$'\t' read -r S _ _; do
    L=$(ls "$OUT/rsem/$S.temp/"*Log.final.out "$OUT/rsem/$S"*Log.final.out 2>/dev/null | head -n1 || true)
    [[ -n "$L" ]] || continue
    awk -F'|' -v s="$S" '/Number of input reads/{n=$2} /Uniquely mapped reads %/{u=$2} /% of reads mapped to multiple loci/{m=$2}
        END{gsub(/[ \t%]/,"",n); gsub(/[ \t%]/,"",u); gsub(/[ \t%]/,"",m); print s"\t"n"\t"u"\t"m}' "$L" >> "$OUT/matrix/star_mapping.tsv"
  done < "$OUT/samples.tsv"
  log "  matrices: $OUT/matrix/genes.expected_count.tsv, genes.TPM.tsv"
fi

log "finished steps: $STEPS"