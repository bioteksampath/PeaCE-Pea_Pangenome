#!/usr/bin/env bash
################################################################################
# snp_calling_qc.sh
# SNP calling and QC for short-read (paired-end) resequencing data
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : SNP calling, population structure and genetic diversity
# Figure     : Fig. 1a (PCA)
#
# Flow
#   step 0  Index reference                    bwa index, samtools faidx
#   step 1  Read trimming (R1 + R2)            cutadapt
#   step 2  Mapping, sorting, duplicate mark   bwa mem -> samtools fixmate/sort/markdup -> CRAM
#   step 3  Mapping statistics                 samtools flagstat / coverage
#   step 4  Joint SNP calling per region       bcftools mpileup | bcftools call -mv
#   step 5  Merge regions per chromosome       bcftools concat
#   step 6  SNP filtering                      biallelic SNPs, MAF, heterozygosity, missing rate
#   step 7  Imputation                         Beagle
#   step 8  Sample missingness + PCA           PLINK2
#
# Usage
#   bash snp_calling_qc.sh -r ref.fa -s samples.tsv -o outdir [options]
#
#   -r  reference FASTA                                   (required)
#   -s  sample sheet: sample<TAB>R1.fq.gz<TAB>R2.fq.gz    (required)
#   -o  output directory                                  (required)
#   -t  threads                                [8]
#   -c  chromosome name regex (from .fai)      ['chr[0-9]+$']
#   -w  calling region size (bp)               [5000000]
#   -p  steps to run, comma separated          [0,1,2,3,4,5,6,7,8]
#   -n  cohort name (output prefix)            [cohort]
#   -a  R1 adapter sequence                    [CTGTCTCTTATACACATCT]
#   -A  R2 adapter sequence                    [CTGTCTCTTATACACATCT]
#
# Sample sheet may start with a header line beginning with "sample"; it is skipped.
#
# Example
#   bash snp_calling_qc.sh -r Psat.Amarillo.fa -s cspgr.tsv -o snp_out -t 32 -n CSPGR
#   bash snp_calling_qc.sh -r Psat.Amarillo.fa -s cspgr.tsv -o snp_out -p 6,7,8
#
# Each step skips samples/regions that are already finished (.done files),
# so the script can be re-run after an interruption.
#
# Software (Methods): cutadapt v5.1, BWA-MEM v0.7.19, SAMtools, BCFtools v1.22,
#                     Beagle v5.4, PLINK v2.0
################################################################################
set -euo pipefail

#---------------------------- parameters (defaults) ----------------------------
THREADS=8
CHROM_REGEX='chr[0-9]+$'
WINDOW=5000000
STEPS="0,1,2,3,4,5,6,7,8"
COHORT="cohort"

# Trimming (settings from the production run)
ADAPTER_R1="CTGTCTCTTATACACATCT"   # Tn5/Nextera; override with -a
ADAPTER_R2="CTGTCTCTTATACACATCT"   # override with -A
TRIM_QUAL=15                       # cutadapt -q
NEXTSEQ_TRIM=20                    # --nextseq-trim (2-colour chemistry, NovaSeq X)
MIN_LEN="75:75"                    # -m R1:R2 minimum length after trimming

# Variant filters (Methods)
MIN_QUAL=30          # site QUAL >= 30 (production run)
MAX_QUAL=""          # production run used QUAL <= 100 (per-sample calls); leave empty for joint calls
MIN_DP=1             # site DP >= 1 (production run)
MAF_MIN=0.005        # MAF > 0.5%
HET_MAX=0.01         # heterozygous calls / called genotypes < 1%
MISS_MAX=0.30        # missing genotype rate <= 30%

# Imputation and PCA
BEAGLE="beagle"      # or: "java -Xmx100g -jar /path/to/beagle.jar"
BEAGLE_NE=1000
SAMPLE_MISS_MAX=0.5  # samples with > 50% missing genotypes excluded from PCA
N_PCS=10

usage() { sed -n '2,40p' "$0"; exit 1; }

while getopts "r:s:o:t:c:w:p:n:a:A:h" opt; do
  case $opt in
    r) REF=$OPTARG ;;   s) SAMPLES=$OPTARG ;;  o) OUT=$OPTARG ;;
    t) THREADS=$OPTARG ;; c) CHROM_REGEX=$OPTARG ;; w) WINDOW=$OPTARG ;;
    p) STEPS=$OPTARG ;; n) COHORT=$OPTARG ;;
    a) ADAPTER_R1=$OPTARG ;; A) ADAPTER_R2=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${REF:-}" && -n "${SAMPLES:-}" && -n "${OUT:-}" ]] || usage
REF=$(readlink -f "$REF"); SAMPLES=$(readlink -f "$SAMPLES")

#---------------------------- helpers ------------------------------------------
log()      { echo "[$(date '+%F %T')] $*"; }
run_step() { [[ ",$STEPS," == *",$1,"* ]]; }
done_ok()  { [[ -f "$OUT/.done/$1" ]]; }
mark()     { mkdir -p "$OUT/.done"; touch "$OUT/.done/$1"; }

mkdir -p "$OUT"/{01_trim,02_cram,03_stats,04_regions,05_chrom,06_filter,07_impute,08_pca,lists,logs}
exec > >(tee -a "$OUT/logs/snp_calling_qc.log") 2>&1
log "reference=$REF samples=$SAMPLES out=$OUT threads=$THREADS steps=$STEPS"

# Sample list (skip header, blank and comment lines)
grep -v -e '^#' -e '^[[:space:]]*$' -e '^sample[[:space:]]' "$SAMPLES" > "$OUT/lists/samples.tsv"
N_SAMPLES=$(wc -l < "$OUT/lists/samples.tsv")
log "$N_SAMPLES samples"

#==============================================================================
# step 0  Index reference
#==============================================================================
if run_step 0; then
  log "step 0: index reference"
  [[ -s "$REF.fai" ]] || samtools faidx "$REF"
  [[ -s "$REF.bwt" ]] || bwa index "$REF"
fi

# Chromosomes and calling regions are taken from the .fai
awk -v re="$CHROM_REGEX" '$1 ~ re {print $1}' "$REF.fai" > "$OUT/lists/chroms.txt"
awk -v re="$CHROM_REGEX" -v w="$WINDOW" 'BEGIN{OFS="\t"} $1 ~ re {
    n=0; for (s=1; s<=$2; s+=w) { e=s+w-1; if (e>$2) e=$2; n++
      printf "%s\t%d\t%d\t%s_%04d\n", $1, s, e, $1, n } }' "$REF.fai" > "$OUT/lists/regions.tsv"
[[ -s "$OUT/lists/chroms.txt" ]] || { echo "No chromosomes match: $CHROM_REGEX"; exit 1; }
log "$(wc -l < "$OUT/lists/chroms.txt") chromosomes, $(wc -l < "$OUT/lists/regions.tsv") regions of ${WINDOW} bp"

#==============================================================================
# step 1  Read trimming (paired-end)
#==============================================================================
if run_step 1; then
  log "step 1: cutadapt"
  while IFS=$'\t' read -r S R1 R2; do
    done_ok "trim_$S" && continue
    cutadapt -j "$THREADS" -a "$ADAPTER_R1" -A "$ADAPTER_R2" \
      -q "$TRIM_QUAL" --nextseq-trim "$NEXTSEQ_TRIM" -m "$MIN_LEN" \
      --json "$OUT/03_stats/$S.cutadapt.json" \
      -o "$OUT/01_trim/${S}_R1.trim.fq.gz" -p "$OUT/01_trim/${S}_R2.trim.fq.gz" \
      "$R1" "$R2" > "$OUT/03_stats/$S.cutadapt.log"
    mark "trim_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 2  Mapping, sorting, duplicate marking, read groups -> CRAM
#==============================================================================
if run_step 2; then
  log "step 2: bwa mem + samtools"
  while IFS=$'\t' read -r S _ _; do
    done_ok "map_$S" && continue
    TMP="$OUT/02_cram/tmp_$S"; mkdir -p "$TMP"
    bwa mem -t "$THREADS" -M -R "@RG\tID:$S\tSM:$S\tLB:$S\tPL:ILLUMINA" "$REF" \
        "$OUT/01_trim/${S}_R1.trim.fq.gz" "$OUT/01_trim/${S}_R2.trim.fq.gz" \
      | samtools fixmate -@ 2 -m -u - - \
      | samtools sort -@ "$THREADS" -m 2G -T "$TMP/sort" -u - \
      | samtools markdup -@ "$THREADS" -T "$TMP/md" -f "$OUT/03_stats/$S.markdup.txt" \
          --reference "$REF" -O cram --write-index - "$OUT/02_cram/$S.cram"
    rm -rf "$TMP"
    mark "map_$S"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 3  Mapping statistics
#==============================================================================
if run_step 3; then
  log "step 3: mapping statistics"
  echo -e "sample\tpct_mapped\tpct_properly_paired\tpct_duplicate\tmean_depth" > "$OUT/03_stats/${COHORT}.mapping_summary.tsv"
  while IFS=$'\t' read -r S _ _; do
    C="$OUT/02_cram/$S.cram"
    [[ -s "$OUT/03_stats/$S.flagstat.txt" ]] || samtools flagstat -@ "$THREADS" "$C" > "$OUT/03_stats/$S.flagstat.txt"
    [[ -s "$OUT/03_stats/$S.coverage.tsv" ]] || samtools coverage --reference "$REF" "$C" > "$OUT/03_stats/$S.coverage.tsv"
    MAP=$(grep -m1 ' mapped (' "$OUT/03_stats/$S.flagstat.txt" | sed -E 's/.*\(([0-9.]+)%.*/\1/')
    PP=$(grep -m1 'properly paired (' "$OUT/03_stats/$S.flagstat.txt" | sed -E 's/.*\(([0-9.]+)%.*/\1/')
    DUP=$(awk -F': ' '/^READ:/{r=$2} /^DUPLICATE TOTAL:/{d=$2} END{printf "%.2f", (r>0 ? 100*d/r : 0)}' "$OUT/03_stats/$S.markdup.txt")
    DEP=$(awk -v re="$CHROM_REGEX" 'NR>1 && $1~re {L=$3-$2+1; b+=L; d+=$7*L} END{printf "%.2f", (b>0 ? d/b : 0)}' "$OUT/03_stats/$S.coverage.tsv")
    echo -e "$S\t$MAP\t$PP\t$DUP\t$DEP" >> "$OUT/03_stats/${COHORT}.mapping_summary.tsv"
  done < "$OUT/lists/samples.tsv"
fi

#==============================================================================
# step 4  Joint SNP calling (all samples together), one region at a time
#         Regions are independent: they can be split across jobs if needed.
#==============================================================================
if run_step 4; then
  log "step 4: bcftools mpileup | call"
  cut -f1 "$OUT/lists/samples.tsv" | sed "s|^|$OUT/02_cram/|; s|$|.cram|" > "$OUT/lists/crams.txt"
  ulimit -n $(( N_SAMPLES + 1000 )) 2>/dev/null || log "WARNING: could not raise open-file limit (ulimit -n = $(ulimit -n))"
  while IFS=$'\t' read -r CHR START END ID; do
    done_ok "call_$ID" && continue
    bcftools mpileup --threads "$THREADS" -f "$REF" -r "$CHR:$START-$END" \
        -a FORMAT/AD,FORMAT/DP -d 20000 -b "$OUT/lists/crams.txt" -Ou \
      | bcftools call --threads "$THREADS" -m -v -Ob -o "$OUT/04_regions/$ID.bcf"
    bcftools index -f "$OUT/04_regions/$ID.bcf"
    mark "call_$ID"
  done < "$OUT/lists/regions.tsv"
fi

#==============================================================================
# step 5  Merge regions per chromosome
#==============================================================================
if run_step 5; then
  log "step 5: bcftools concat"
  while read -r CHR; do
    done_ok "concat_$CHR" && continue
    awk -v c="$CHR" -v d="$OUT/04_regions" '$1==c {print d"/"$4".bcf"}' "$OUT/lists/regions.tsv" > "$OUT/05_chrom/$CHR.parts.txt"
    bcftools concat --threads "$THREADS" -f "$OUT/05_chrom/$CHR.parts.txt" \
        -Ob -o "$OUT/05_chrom/${COHORT}.$CHR.raw.bcf" --write-index
    mark "concat_$CHR"
  done < "$OUT/lists/chroms.txt"
fi

#==============================================================================
# step 6  SNP filtering
#         biallelic SNPs; QUAL >= 30; DP >= 1; MAF > 0.5%;
#         het < 1% of called genotypes; missing genotypes <= 30%
#==============================================================================
if run_step 6; then
  log "step 6: SNP filtering"
  FILTER="QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DP} && MAF>${MAF_MIN} && F_MISSING<=${MISS_MAX} && COUNT(GT=\"het\")/(N_SAMPLES-N_MISSING)<${HET_MAX}"
  [[ -n "$MAX_QUAL" ]] && FILTER="$FILTER && QUAL<=${MAX_QUAL}"
  log "filter: $FILTER"
  echo -e "chrom\traw_variants\tfiltered_snps" > "$OUT/06_filter/${COHORT}.snp_counts.tsv"
  while read -r CHR; do
    IN="$OUT/05_chrom/${COHORT}.$CHR.raw.bcf"; FO="$OUT/06_filter/${COHORT}.$CHR.filtered.vcf.gz"
    if ! done_ok "filter_$CHR"; then
      bcftools view --threads "$THREADS" -m2 -M2 -v snps "$IN" -Ou \
        | bcftools view --threads "$THREADS" -i "$FILTER" -Ou \
        | bcftools annotate --threads "$THREADS" -x INFO,^FORMAT/GT \
            --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$FO" --write-index=tbi
      mark "filter_$CHR"
    fi
    echo -e "$CHR\t$(bcftools index -n "$IN")\t$(bcftools index -n "$FO")" >> "$OUT/06_filter/${COHORT}.snp_counts.tsv"
  done < "$OUT/lists/chroms.txt"
fi

#==============================================================================
# step 7  Imputation (Beagle), per chromosome
#==============================================================================
if run_step 7; then
  log "step 7: Beagle imputation"
  while read -r CHR; do
    done_ok "impute_$CHR" && continue
    $BEAGLE gt="$OUT/06_filter/${COHORT}.$CHR.filtered.vcf.gz" \
            out="$OUT/07_impute/${COHORT}.$CHR.imputed" ne="$BEAGLE_NE" nthreads="$THREADS"
    bcftools index -f -t "$OUT/07_impute/${COHORT}.$CHR.imputed.vcf.gz"
    mark "impute_$CHR"
  done < "$OUT/lists/chroms.txt"
fi

#==============================================================================
# step 8  Merge chromosomes, sample missingness, PCA
#         Missingness is measured before imputation (imputed data have none);
#         samples with > 50% missing genotypes are excluded from the PCA.
#==============================================================================
if run_step 8; then
  log "step 8: PLINK2 PCA"
  P="$OUT/08_pca/${COHORT}"
  sed "s|^|$OUT/06_filter/${COHORT}.|; s|$|.filtered.vcf.gz|" "$OUT/lists/chroms.txt" > "$P.filtered.list"
  sed "s|^|$OUT/07_impute/${COHORT}.|; s|$|.imputed.vcf.gz|"  "$OUT/lists/chroms.txt" > "$P.imputed.list"
  bcftools concat --threads "$THREADS" -f "$P.filtered.list" -Oz -o "$P.filtered.vcf.gz" --write-index=tbi
  bcftools concat --threads "$THREADS" -f "$P.imputed.list"  -Oz -o "$P.imputed.vcf.gz"  --write-index=tbi

  plink2 --vcf "$P.filtered.vcf.gz" --allow-extra-chr --missing sample-only \
         --threads "$THREADS" --out "$P.missing"
  awk -v m="$SAMPLE_MISS_MAX" 'NR==1{for(i=1;i<=NF;i++) if($i=="F_MISS") c=i; next} $c<=m {print $1}' \
      "$P.missing.smiss" > "$P.keep.txt"
  log "PCA on $(wc -l < "$P.keep.txt") of $N_SAMPLES samples"

  plink2 --vcf "$P.imputed.vcf.gz" --allow-extra-chr --keep "$P.keep.txt" \
         --max-alleles 2 --pca "$N_PCS" --threads "$THREADS" --out "$P.pca"
fi

log "finished steps: $STEPS"