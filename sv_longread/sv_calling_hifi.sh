#!/usr/bin/env bash
################################################################################
# sv_calling_hifi.sh
# Structural variant (SV) calling from PacBio HiFi reads with two callers
# (Sniffles2 + SVIM). Only consensus SVs found by BOTH callers are kept.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Long-read structural variant detection
# Figures    : Extended Data Fig. 5a-d; input for TIP discovery (tip_discovery)
#
# Flow
#   step 0  Index reference                     samtools faidx
#   step 1  Map HiFi reads, sort, mark dups     minimap2 -ax map-hifi -> samtools -> BAM
#   step 2  SV calling, caller 1                Sniffles2 (PASS calls only) + .snf file
#   step 3  SV calling, caller 2                SVIM
#   step 4  Consensus per sample                SURVIVOR merge (both callers must agree)
#   step 5  Per-sample statistics and BED       SURVIVOR stats / vcftobed
#   step 6  Pan-SV catalogue (all samples)      SURVIVOR merge
#   step 7  Summary table                       SV counts by type per sample
#
# Usage
#   bash sv_calling_hifi.sh -r ref.fa -s samples.tsv -o outdir [options]
#
#   -r  reference FASTA (e.g. CDC Amarillo)                       (required)
#   -s  sample sheet: sample<TAB>reads                            (required)
#       reads = HiFi FASTQ/FASTA(.gz), or an existing sorted BAM/CRAM
#       (BAM/CRAM input skips step 1 for that sample)
#   -o  output directory                                          (required)
#   -t  threads                                  [16]
#   -p  steps to run, comma separated            [0,1,2,3,4,5,6,7]
#   -n  cohort name (output prefix)              [cohort]
#   -l  minimum SV length (bp)                   [30]
#   -d  SVIM minimum read depth                  [10]
#
# Example
#   bash sv_calling_hifi.sh -r Psat.Amarillo.fa -s hifi.tsv -o sv_out -t 64 -n Pisum85
#   bash sv_calling_hifi.sh -r Psat.Amarillo.fa -s hifi.tsv -o sv_out -p 6,7
#
# Each step skips samples that are already finished (.done files).
#
# Software (Methods): minimap2 v2.28, SAMtools, Sniffles2 v2.6.2, SVIM v2.0.0,
#                     SURVIVOR v1.0.7, BCFtools
################################################################################
set -euo pipefail

#---------------------------- parameters (defaults) ----------------------------
THREADS=16
STEPS="0,1,2,3,4,5,6,7"
COHORT="cohort"
MIN_SV=30                  # minimum SV length (Sniffles2 --minsvlen, SVIM --min_sv_size)
SVIM_MIN_DEPTH=10          # SVIM --minimum_depth
SVIM_MIN_QUAL=0            # optional SVIM QUAL filter (0 = keep all, as in production run)
MAX_BED_LEN=10000000       # SURVIVOR vcftobed max SV length (10 Mb)

# SURVIVOR merge, per sample (consensus of the two callers)
#   max breakpoint distance, min callers, same type, same strand, estimate distance, min size
SAMPLE_MERGE="1000 2 1 0 0 0"
# SURVIVOR merge, across samples (pan-SV catalogue)
COHORT_MERGE="1000 1 1 0 0 30"

usage() { sed -n '2,42p' "$0"; exit 1; }

while getopts "r:s:o:t:p:n:l:d:h" opt; do
  case $opt in
    r) REF=$OPTARG ;; s) SAMPLES=$OPTARG ;; o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;; p) STEPS=$OPTARG ;; n) COHORT=$OPTARG ;;
    l) MIN_SV=$OPTARG ;; d) SVIM_MIN_DEPTH=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${REF:-}" && -n "${SAMPLES:-}" && -n "${OUT:-}" ]] || usage
REF=$(readlink -f "$REF"); SAMPLES=$(readlink -f "$SAMPLES")

#---------------------------- helpers ------------------------------------------
log()      { echo "[$(date '+%F %T')] $*"; }
run_step() { [[ ",$STEPS," == *",$1,"* ]]; }
done_ok()  { [[ -f "$OUT/.done/$1" ]]; }
mark()     { mkdir -p "$OUT/.done"; touch "$OUT/.done/$1"; }
aln_of()   {   # alignment file used for a sample
  case "$2" in *.bam|*.cram) readlink -f "$2" ;; *) echo "$OUT/01_align/$1.bam" ;; esac
}

mkdir -p "$OUT"/{01_align,02_sniffles,03_svim,04_consensus,05_stats,06_cohort,lists,logs}
exec > >(tee -a "$OUT/logs/sv_calling_hifi.log") 2>&1
log "reference=$REF samples=$SAMPLES out=$OUT threads=$THREADS steps=$STEPS min_sv=$MIN_SV"

grep -v -e '^#' -e '^[[:space:]]*$' -e '^sample[[:space:]]' "$SAMPLES" > "$OUT/lists/samples.tsv"
log "$(wc -l < "$OUT/lists/samples.tsv") samples"

#==============================================================================
# step 0  Index reference
#==============================================================================
if run_step 0; then
  log "step 0: index reference"
  [[ -s "$REF.fai" ]] || samtools faidx "$REF"
fi

#==============================================================================
# step 1  Map HiFi reads -> sorted, duplicate-marked BAM with read group
#==============================================================================
if run_step 1; then
  log "step 1: minimap2 map-hifi"
  while IFS=$'\t' read -r S READS; do
    case "$READS" in *.bam|*.cram) log "  $S: using existing alignment $READS"; continue ;; esac
    done_ok "map_$S" && continue
    TMP="$OUT/01_align/tmp_$S"; mkdir -p "$TMP"
    minimap2 -t "$THREADS" -ax map-hifi -R "@RG\tID:$S\tSM:$S\tPL:PACBIO" "$REF" "$READS" \
      2> "$OUT/logs/$S.minimap2.log" \
      | samtools sort -@ "$THREADS" -m 4G -T "$TMP/sort" -u - \
      | samtools markdup -@ "$THREADS" -T "$TMP/md" --write-index - "$OUT/01_align/$S.bam"
    rm -rf "$TMP"
    samtools flagstat -@ "$THREADS" "$OUT/01_align/$S.bam" > "$OUT/05_stats/$S.flagstat.txt"
    mark "map_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 2  Sniffles2 (keep PASS calls; .snf kept for later multi-sample / force calling)
#==============================================================================
if run_step 2; then
  log "step 2: Sniffles2"
  while IFS=$'\t' read -r S READS; do
    done_ok "sniffles_$S" && continue
    A=$(aln_of "$S" "$READS")
    sniffles -i "$A" --reference "$REF" -t "$THREADS" --minsvlen "$MIN_SV" \
             --sample-id "$S" -v "$OUT/02_sniffles/$S.sniffles.raw.vcf" \
             --snf "$OUT/02_sniffles/$S.snf" > "$OUT/logs/$S.sniffles.log" 2>&1
    bcftools view -f PASS "$OUT/02_sniffles/$S.sniffles.raw.vcf" -Ov -o "$OUT/02_sniffles/$S.sniffles.vcf"
    mark "sniffles_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 3  SVIM
#==============================================================================
if run_step 3; then
  log "step 3: SVIM"
  while IFS=$'\t' read -r S READS; do
    done_ok "svim_$S" && continue
    A=$(aln_of "$S" "$READS")
    svim alignment "$OUT/03_svim/$S" "$A" "$REF" --sample "$S" \
         --min_sv_size "$MIN_SV" --minimum_depth "$SVIM_MIN_DEPTH" > "$OUT/logs/$S.svim.log" 2>&1
    bcftools view -i "QUAL>=${SVIM_MIN_QUAL}" "$OUT/03_svim/$S/variants.vcf" -Ov -o "$OUT/03_svim/$S.svim.vcf"
    mark "svim_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 4  Consensus per sample: SV must be called by both Sniffles2 and SVIM
#         SURVIVOR merge <list> 1000 2 1 0 0 0
#==============================================================================
if run_step 4; then
  log "step 4: SURVIVOR consensus ($SAMPLE_MERGE)"
  while IFS=$'\t' read -r S _; do
    done_ok "consensus_$S" && continue
    printf '%s\n' "$OUT/02_sniffles/$S.sniffles.vcf" "$OUT/03_svim/$S.svim.vcf" > "$OUT/04_consensus/$S.vcfs.txt"
    # shellcheck disable=SC2086
    SURVIVOR merge "$OUT/04_consensus/$S.vcfs.txt" $SAMPLE_MERGE "$OUT/04_consensus/$S.consensus.vcf" \
      > "$OUT/logs/$S.survivor.log" 2>&1
    mark "consensus_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 5  Per-sample statistics and BED (30 bp - 10 Mb)
#==============================================================================
if run_step 5; then
  log "step 5: SURVIVOR stats / vcftobed"
  while IFS=$'\t' read -r S _; do
    V="$OUT/04_consensus/$S.consensus.vcf"
    SURVIVOR stats "$V" -1 -1 -1 "$OUT/05_stats/$S.consensus.stats" > /dev/null
    SURVIVOR vcftobed "$V" "$MIN_SV" "$MAX_BED_LEN" "$OUT/05_stats/$S.consensus.bedpe" > /dev/null
    # chrom, start, end, SV type
    awk 'BEGIN{OFS="\t"} {print $1, $2, $5, $11}' "$OUT/05_stats/$S.consensus.bedpe" > "$OUT/05_stats/$S.consensus.bed"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 6  Pan-SV catalogue: merge consensus SVs across all samples
#==============================================================================
if run_step 6; then
  log "step 6: SURVIVOR cohort merge ($COHORT_MERGE)"
  cut -f1 "$OUT/lists/samples.tsv" | sed "s|^|$OUT/04_consensus/|; s|$|.consensus.vcf|" > "$OUT/06_cohort/consensus_vcfs.txt"
  # shellcheck disable=SC2086
  SURVIVOR merge "$OUT/06_cohort/consensus_vcfs.txt" $COHORT_MERGE "$OUT/06_cohort/${COHORT}.panSV.vcf" \
    > "$OUT/logs/${COHORT}.survivor_cohort.log" 2>&1
  bcftools sort "$OUT/06_cohort/${COHORT}.panSV.vcf" -Oz -o "$OUT/06_cohort/${COHORT}.panSV.vcf.gz"
  bcftools index -f -t "$OUT/06_cohort/${COHORT}.panSV.vcf.gz"
  SURVIVOR stats "$OUT/06_cohort/${COHORT}.panSV.vcf" -1 -1 -1 "$OUT/06_cohort/${COHORT}.panSV.stats" > /dev/null
fi

#==============================================================================
# step 7  Summary: consensus SV counts by type per sample (Extended Data Fig. 5a)
#==============================================================================
if run_step 7; then
  log "step 7: summary table"
  T="$OUT/05_stats/${COHORT}.sv_counts.tsv"
  echo -e "sample\tsniffles_PASS\tsvim\tconsensus\tDEL\tINS\tDUP\tINV\tBND\tother" > "$T"
  while IFS=$'\t' read -r S _; do
    n_snf=$(grep -vc '^#' "$OUT/02_sniffles/$S.sniffles.vcf" || true)
    n_svm=$(grep -vc '^#' "$OUT/03_svim/$S.svim.vcf" || true)
    { grep -v '^#' "$OUT/04_consensus/$S.consensus.vcf" || true; } \
      | awk -v s="$S" -v a="$n_snf" -v b="$n_svm" 'BEGIN{OFS="\t"}
          { t="other"; if (match($8, /SVTYPE=[A-Z]+/)) t=substr($8, RSTART+7, RLENGTH-7)
            if (t ~ /^DUP/) t="DUP"; if (t=="TRA") t="BND"
            c[t]++; n++ }
          END{ print s, a, b, n+0, c["DEL"]+0, c["INS"]+0, c["DUP"]+0, c["INV"]+0, c["BND"]+0, c["other"]+0 }' >> "$T"
  done < "$OUT/lists/samples.tsv"
  log "summary: $T"
fi

log "finished steps: $STEPS"