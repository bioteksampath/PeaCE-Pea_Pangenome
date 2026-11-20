#!/usr/bin/env bash
################################################################################
# tip_discovery.sh
# Transposon insertion polymorphism (TIP) discovery from the pan-SV catalogue.
# Insertions/deletions are compared with a TE library by BLAST; SVs that are
# mostly TE sequence are called TIPs, classified to lineage (TEsorter), dated
# (LTR-retrotransposons) and written to a TIP VCF.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Transposon insertion polymorphism (TIP) discovery and genotyping
# Figures    : Fig. 3a-e; Extended Data Fig. 6a-f; Supplementary Table 6
# Based on   : https://github.com/Lcamdom/panTEvo/tree/main/04_pangenome_TIP_annotation
#
# Flow
#   step 1  SV sequences -> FASTA            DEL (REF allele) and INS (ALT allele), length filter
#   step 2  BLAST TE library vs SVs          blastn (identity >= 80%, TE coverage >= 80%)
#   step 3  Select TIPs                      alignment > 80% of SV length; best TE hit per SV
#   step 4  Lineage classification           TEsorter (rexdb-plant)
#   step 5  LTR insertion age                tip_ltr_age.py (self-BLAST + MAFFT + K2P)
#   step 6  TIP annotation table             SV + best TE hit + lineage + age
#   step 7  TIP VCF                          annotated, sorted, unique IDs (+ sites-only matrix)
#   step 8  Force-call TIPs (optional)       bcftools norm (left-normalize) + Sniffles2 --genotype-vcf
#
# Usage
#   bash tip_discovery.sh -v panSV.vcf -t TE_library.fa -o tip_out [options]
#
#   -v  pan-SV VCF with sequence-resolved DEL/INS (sv_longread, step 6)     (required)
#   -t  TE library FASTA (EDTA consensus + EDTA intact + DANTE_LTR)     (required)
#   -o  output directory                                                (required)
#   -r  reference FASTA (needed for step 8)
#   -b  alignments for step 8: sample<TAB>sorted.bam
#   -T  threads                                     [16]
#   -m  minimum SV length (bp)                      [200]
#   -M  maximum SV length (bp)                      [40000]
#   -i  minimum BLAST identity (%)                  [80]
#   -c  minimum fraction of SV covered by TE hit    [0.8]
#   -p  steps to run, comma separated               [1,2,3,4,5,6,7]
#   -n  output prefix                               [Pisum]
#
# Example
#   bash tip_discovery.sh -v Pisum85.panSV.vcf.gz -t Pisum_TE_library.fa -o tips -T 36
#   bash tip_discovery.sh -v Pisum85.panSV.vcf.gz -t Pisum_TE_library.fa -o tips \
#        -r Psat.Amarillo.fa -b hifi_bams.tsv -p 8
#
# Software: BCFtools, seqkit, BLAST+ v2.12.0, BEDTools, TEsorter v1.4.7, MAFFT,
#           Python 3 + Biopython, Sniffles2 v2.7.2 (step 8)
################################################################################
set -euo pipefail

#---------------------------- parameters (defaults) ----------------------------
THREADS=16
MIN_LEN=200          # Methods: insertions or deletions longer than 200 bp
MAX_LEN=40000
MIN_ID=80
MIN_SV_COV=0.8
STEPS="1,2,3,4,5,6,7"
PREFIX="Pisum"
TES_COV=10            # TEsorter -cov
TES_EVAL=1e-2         # TEsorter -eval
TES_DB="rexdb-plant"
LTR_RATE=1.3e-8       # substitutions/site/year for LTR age
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() { sed -n '2,48p' "$0"; exit 1; }
while getopts "v:t:o:r:b:T:m:M:i:c:p:n:h" opt; do
  case $opt in
    v) VCF=$OPTARG ;; t) TELIB=$OPTARG ;; o) OUT=$OPTARG ;; r) REF=$OPTARG ;; b) BAMS=$OPTARG ;;
    T) THREADS=$OPTARG ;; m) MIN_LEN=$OPTARG ;; M) MAX_LEN=$OPTARG ;; i) MIN_ID=$OPTARG ;;
    c) MIN_SV_COV=$OPTARG ;; p) STEPS=$OPTARG ;; n) PREFIX=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${VCF:-}" && -n "${TELIB:-}" && -n "${OUT:-}" ]] || usage
VCF=$(readlink -f "$VCF"); TELIB=$(readlink -f "$TELIB")

log()      { echo "[$(date '+%F %T')] $*"; }
run_step() { [[ ",$STEPS," == *",$1,"* ]]; }

mkdir -p "$OUT"/{01_sv_fasta,02_blast,03_tips,04_tesorter,05_ltr_age,06_annotation,07_vcf,08_forcecall,logs}
exec > >(tee -a "$OUT/logs/tip_discovery.log") 2>&1
log "vcf=$VCF te_library=$TELIB out=$OUT len=${MIN_LEN}-${MAX_LEN} identity>=${MIN_ID} sv_cov>${MIN_SV_COV}"

SVMAP="$OUT/01_sv_fasta/sv_map.tsv"            # key CHROM POS ID SVTYPE SVLEN
SVFA="$OUT/01_sv_fasta/${PREFIX}.SV.fa"
HITS="$OUT/02_blast/${PREFIX}.TE_vs_SV.blast.tsv"
TIPHITS="$OUT/03_tips/${PREFIX}.TIP_hits.tsv"
TIPLIST="$OUT/03_tips/${PREFIX}.TIP_keys.txt"
TIPFA="$OUT/03_tips/${PREFIX}.TIPs.fa"
CLS="$OUT/04_tesorter/${PREFIX}.TIPs.fa.${TES_DB}.cls.tsv"
AGE="$OUT/05_ltr_age/LTR_results.csv"
ANNOT="$OUT/06_annotation/${PREFIX}.TIP_annotation.tsv"

#==============================================================================
# step 1  SV sequences -> FASTA
#         DEL: deleted sequence = REF allele; INS: inserted sequence = ALT allele.
#         Each SV gets a unique key (sv0000001 ...) linked to CHROM/POS/ID in sv_map.tsv.
#==============================================================================
if run_step 1; then
  log "step 1: extract DEL/INS sequences"
  bcftools query -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' \
      -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%REF\t%ALT\n' "$VCF" \
    | awk -F'\t' -v map="$SVMAP" -v fa="$OUT/01_sv_fasta/all_SV.fa" 'BEGIN{OFS="\t"}
        { seq = ($4 == "DEL") ? $5 : $6
          k = sprintf("sv%07d", NR)
          if (seq ~ /^[ACGTNacgtn]+$/) { print ">" k "\n" seq > fa; print k, $1, $2, $3, $4, length(seq) > map; ok++ }
          else nos++ }
        END{ printf "  SVs with sequence: %d | symbolic/no sequence (skipped): %d\n", ok, nos > "/dev/stderr" }'
  seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" "$OUT/01_sv_fasta/all_SV.fa" > "$SVFA"
  log "  SVs ${MIN_LEN}-${MAX_LEN} bp: $(grep -c '>' "$SVFA")"
fi

#==============================================================================
# step 2  BLAST: TE library (query) vs SV sequences (database)
#==============================================================================
if run_step 2; then
  log "step 2: blastn TE library vs SVs"
  makeblastdb -in "$SVFA" -dbtype nucl -out "$OUT/02_blast/SVdb" > "$OUT/logs/makeblastdb.log"
  blastn -query "$TELIB" -db "$OUT/02_blast/SVdb" \
    -outfmt "6 qseqid sseqid qstart qend sstart send pident length qlen slen qcovs bitscore" \
    -perc_identity "$MIN_ID" -qcov_hsp_perc "$(awk -v c="$MIN_SV_COV" 'BEGIN{print c*100}')" \
    -num_threads "$THREADS" -out "$HITS"
  log "  hits: $(wc -l < "$HITS")"
fi

#==============================================================================
# step 3  TIPs: TE hit covers > 80% of the SV (and >= 80% of the TE, step 2)
#==============================================================================
if run_step 3; then
  log "step 3: select TIPs"
  awk -F'\t' -v c="$MIN_SV_COV" '$8 > c * $10' "$HITS" > "$TIPHITS"
  cut -f2 "$TIPHITS" | sort -u > "$TIPLIST"
  # hits per SV (BED, 0-based) and merged overlapping hits
  awk 'BEGIN{OFS="\t"} { s=$5; e=$6; st="+"; if (s > e) { t=s; s=e; e=t; st="-" }
       print $2, s-1, e, $1, $12, st }' "$TIPHITS" \
    | sort -k1,1 -k2,2n > "$OUT/03_tips/${PREFIX}.TIP_hits.bed"
  bedtools merge -i "$OUT/03_tips/${PREFIX}.TIP_hits.bed" -c 4,4 -o count,first \
    > "$OUT/03_tips/${PREFIX}.TIP_hits.merged.bed"
  seqkit grep -f "$TIPLIST" "$SVFA" > "$TIPFA"
  log "  TIPs: $(wc -l < "$TIPLIST") of $(grep -c '>' "$SVFA") SVs"
fi

#==============================================================================
# step 4  Lineage classification (TEsorter)
#==============================================================================
if run_step 4; then
  log "step 4: TEsorter"
  ( cd "$OUT/04_tesorter" && ln -sf "$TIPFA" "${PREFIX}.TIPs.fa" &&
    TEsorter "${PREFIX}.TIPs.fa" -p "$THREADS" -cov "$TES_COV" -eval "$TES_EVAL" -db "$TES_DB" \
      > "$OUT/logs/tesorter.log" 2>&1 )
  log "  classified: $(grep -vc '^#' "$CLS" || true)"
fi

#==============================================================================
# step 5  Insertion age of LTR-retrotransposon TIPs
#==============================================================================
if run_step 5; then
  log "step 5: LTR insertion age"
  awk -F'\t' '!/^#/ && $2 == "LTR" {print $1}' "$CLS" | sort -u > "$OUT/05_ltr_age/LTR_TIP_keys.txt"
  seqkit grep -f "$OUT/05_ltr_age/LTR_TIP_keys.txt" "$TIPFA" > "$OUT/05_ltr_age/${PREFIX}.LTR_TIPs.fa"
  log "  LTR TIPs: $(wc -l < "$OUT/05_ltr_age/LTR_TIP_keys.txt")"
  python3 "$SCRIPT_DIR/tip_ltr_age.py" "$OUT/05_ltr_age/${PREFIX}.LTR_TIPs.fa" "$OUT/05_ltr_age" \
      --cores "$THREADS" --rate "$LTR_RATE" --quiet
fi

#==============================================================================
# step 6  TIP annotation table
#         key CHROM POS ID SVTYPE SVLEN best_TE_hit pident TE_order superfamily clade K2P age_Myr
#==============================================================================
if run_step 6; then
  log "step 6: annotation table"
  awk -F'\t' 'BEGIN{OFS="\t"}
      FILENAME == ARGV[1] { if (!($2 in bs) || $12 > bs[$2]) { bs[$2]=$12; te[$2]=$1; id[$2]=$7 }; next }
      FILENAME == ARGV[2] { if ($0 !~ /^#/) { ord[$1]=$2; sf[$1]=$3; cl[$1]=$4 }; next }
      FILENAME == ARGV[3] { if (FNR > 1) { split($0, a, ","); k2[a[1]]=a[4]; ag[a[1]]=a[5] }; next }
      FNR == 1 { print "key","CHROM","POS","ID","SVTYPE","SVLEN","best_TE_hit","pident","TE_order","TE_superfamily","TE_clade","K2P","age_Myr" }
      ($1 in te) { print $1, $2, $3, $4, $5, $6, te[$1], id[$1],
                   ($1 in ord ? ord[$1] : "unknown"), ($1 in sf ? sf[$1] : "unknown"), ($1 in cl ? cl[$1] : "unknown"),
                   ($1 in k2 ? k2[$1] : "NA"), ($1 in ag ? ag[$1] : "NA") }' \
      "$TIPHITS" "$CLS" <( [[ -s "$AGE" ]] && cat "$AGE" || echo "ID" ) "$SVMAP" > "$ANNOT"
  log "  annotation: $ANNOT ($(($(wc -l < "$ANNOT") - 1)) TIPs)"
  # summary by TE lineage
  awk -F'\t' 'NR > 1 {n[$5"\t"$9"\t"$10"\t"$11]++} END{for (k in n) print k"\t"n[k]}' "$ANNOT" \
    | sort -k5,5nr | sed '1i SVTYPE\tTE_order\tTE_superfamily\tTE_clade\tn_TIPs' > "$OUT/06_annotation/${PREFIX}.TIP_lineage_counts.tsv"
fi

#==============================================================================
# step 7  TIP VCF: TIP records from the pan-SV VCF, TE annotation in INFO,
#         sorted, unique IDs; plus a sites-only matrix (no genotypes)
#==============================================================================
if run_step 7; then
  log "step 7: TIP VCF"
  V="$OUT/07_vcf/${PREFIX}.TIPs"
  {
    bcftools view -h "$VCF" | sed '$d'
    echo '##INFO=<ID=TE_ORDER,Number=1,Type=String,Description="TE order of the TIP (TEsorter)">'
    echo '##INFO=<ID=TE_LINEAGE,Number=1,Type=String,Description="TE superfamily_clade of the TIP (TEsorter rexdb-plant)">'
    echo '##INFO=<ID=TE_HIT,Number=1,Type=String,Description="Best TE library hit (BLASTN)">'
    echo '##INFO=<ID=TE_AGE,Number=1,Type=String,Description="LTR insertion age (Myr; K2P, mu=1.3e-8)">'
    bcftools view -h "$VCF" | tail -n1
    bcftools view -H -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS"' "$VCF" \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          FNR == NR { if (FNR > 1) { key=$2"\t"$3"\t"$4; ann[key]="TE_ORDER="$9";TE_LINEAGE="$10"_"$11";TE_HIT="$7";TE_AGE="$13 }; next }
          { key=$1"\t"$2"\t"$3; if (key in ann) { $8 = $8";"ann[key]; print } }' "$ANNOT" -
  } | bcftools sort -Ov -o "$V.tmp.vcf" -
  # unique IDs: <original ID>.SV<n>
  awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$3 = $3 ".SV" ++i; print}' "$V.tmp.vcf" \
    | bcftools view -Oz -o "$V.vcf.gz"
  bcftools index -f -t "$V.vcf.gz"
  rm -f "$V.tmp.vcf"
  bcftools view -G "$V.vcf.gz" -Oz -o "$V.sites.vcf.gz"
  bcftools index -f -t "$V.sites.vcf.gz"
  log "  TIP VCF: $V.vcf.gz ($(bcftools index -n "$V.vcf.gz") records); sites-only: $V.sites.vcf.gz"
fi

#==============================================================================
# step 8  (optional) Force-call TIPs in each long-read genome (Sniffles2)
#==============================================================================
if run_step 8; then
  [[ -n "${REF:-}" && -n "${BAMS:-}" ]] || { echo "step 8 needs -r reference and -b sample<TAB>bam list"; exit 1; }
  log "step 8: Sniffles2 force-calling"
  # site-only input VCF: left-normalized against the reference and position-sorted (Methods)
  SITES="$OUT/08_forcecall/${PREFIX}.TIPs.sites.norm.vcf.gz"
  if [[ ! -s "$SITES" ]]; then
    bcftools norm -f "$REF" -c w "$OUT/07_vcf/${PREFIX}.TIPs.sites.vcf.gz" -Ou 2> "$OUT/logs/bcftools_norm.log" \
      | bcftools sort -Oz -o "$SITES" -
    bcftools index -f -t "$SITES"
  fi
  : > "$OUT/08_forcecall/vcfs.txt"
  while IFS=$'\t' read -r S BAM; do
    [[ -z "$S" || "$S" == sample ]] && continue
    G="$OUT/08_forcecall/$S.TIPgt.vcf.gz"
    [[ -s "$G" ]] || sniffles --input "$BAM" --reference "$REF" --genotype-vcf "$SITES" \
                             --sample-id "$S" -t "$THREADS" --vcf "$G" > "$OUT/logs/$S.forcecall.log" 2>&1
    [[ -s "$G.tbi" ]] || bcftools index -f -t "$G"
    echo "$G" >> "$OUT/08_forcecall/vcfs.txt"
  done < "$BAMS"
  bcftools merge -l "$OUT/08_forcecall/vcfs.txt" -m none -Oz -o "$OUT/08_forcecall/${PREFIX}.TIPs.forcecalled.vcf.gz"
  bcftools index -f -t "$OUT/08_forcecall/${PREFIX}.TIPs.forcecalled.vcf.gz"
  log "  genotyped TIP VCF: $OUT/08_forcecall/${PREFIX}.TIPs.forcecalled.vcf.gz"
fi

log "finished steps: $STEPS"