#!/usr/bin/env bash
################################################################################
# ltr_age_pipeline.sh
# Full-length LTR retrotransposons from DANTE_LTR and EDTA, insertion age from
# 5' vs 3' LTR divergence, and one non-redundant LTR age table per genome and
# for all genomes.
#
# Manuscript : Pea pangenomics connects structural variation with disease
# Methods    : Identification and classification of full-length LTR
#              retrotransposons; Estimation of insertion ages
# Figures    : Fig. 3a (LTR-RT density); Extended Data Fig. 2b, 2d; Supplementary Fig. 3
#
# Flow (per genome, driven by one manifest)
#   step 1  EDTA table       EDTA *.pass.list + TEsorter *.cls.tsv -> chrom, start, end, strand, family, age (years)
#   step 2  DANTE table      DANTE_LTR GFF3 -> intact elements, unique IDs, strand, 5'/3' LTR coordinates
#   step 3  DANTE age        5' vs 3' LTR BLASTN divergence -> raw_T, K2P_T, JC69_T (years)
#   step 4  EDTA-unique      EDTA elements with < 50% overlap with DANTE; 5'/3' LTR positions by realignment
#   step 5  Merge            DANTE + EDTA-unique per genome; all-genome table; age in Mya; age group
#
# Usage
#   bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age [options]
#
#   -m  genome manifest (build_ltr_manifest.py)           (required)
#   -o  output directory                                           (required)
#   -g  genome IDs to run, comma separated (e.g. G1,G4)            [all in manifest]
#   -p  steps, comma separated                                     [1,2,3,4,5]
#   -t  threads                                                    [16]
#   -H  folder with the age helper scripts                         [<script dir>/lib]
#         make_ltr_inputs.py, ltr_divergence_dnateV2.py, edta_infer_ltr_positions.py
#
# Examples
#   python build_ltr_manifest.py --map1 code/map1_*.map --map2 code/map2_dante_table-age_*.txt \
#          --edta-dir-template /gpfs/.../data/EDTA/{GENOME_CODE} --out genome_manifest.tsv
#   bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -t 32               # all genomes, all steps
#   bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -g G4 -p 2,3,4      # one genome
#   bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -p 5                 # final merge only
#
# Software: BLAST+, BEDTools, Python 3 + Biopython (helpers), awk
################################################################################
set -euo pipefail

#---------------------------- parameters ---------------------------------------
THREADS=16
STEPS="1,2,3,4,5"
GIDS=""
HELPERS=""
# step 3: DANTE 5' vs 3' LTR BLAST (production values)
BLAST_TASK="blastn"; WORD_SIZE=7; EVALUE="1e-5"; MIN_ALIGN=50; MAX_OFF=120; END_WINDOW=300
# step 4: EDTA-unique elements
EDTA_DANTE_OVERLAP=0.5                 # removed if >= 50% of the EDTA element overlaps a DANTE element
INFER_WINDOW=4000; INFER_MIN_ALIGN=200; INFER_MAX_OFF=80; INFER_PIDENT_MIN=75
# step 5: age column and age groups (Mya) shared with all downstream figure scripts
AGE_COLUMN="K2P_T"

usage() { sed -n '2,38p' "$0"; exit 1; }
while getopts "m:o:g:p:t:H:h" opt; do
  case $opt in
    m) MANIFEST=$OPTARG ;; o) OUT=$OPTARG ;; g) GIDS=$OPTARG ;; p) STEPS=$OPTARG ;;
    t) THREADS=$OPTARG ;; H) HELPERS=$OPTARG ;; *) usage ;;
  esac
done
[[ -n "${MANIFEST:-}" && -n "${OUT:-}" ]] || usage
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
HELPERS=${HELPERS:-$SCRIPT_DIR/lib}
MANIFEST=$(readlink -f "$MANIFEST")
if grep -q $'\r' "$MANIFEST"; then          # manifest edited on Windows: drop carriage returns
  mkdir -p "$OUT"; tr -d '\r' < "$MANIFEST" > "$OUT/manifest.unix.tsv"; MANIFEST="$OUT/manifest.unix.tsv"
fi

log()      { echo "[$(date '+%F %T')] $*"; }
run_step() { [[ ",$STEPS," == *",$1,"* ]]; }
is_true()  { [[ "${1,,}" =~ ^(1|true|t|yes)$ ]]; }

mkdir -p "$OUT"/{01_edta_table,02_dante_table,03_dante_age,04_edta_unique,05_final,logs}
exec > >(tee -a "$OUT/logs/ltr_age_pipeline.log") 2>&1
log "manifest=$MANIFEST out=$OUT steps=$STEPS genomes=${GIDS:-all}"

# manifest columns are read by NAME
mapfile -t MCOLS < <(head -n1 "$MANIFEST" | tr '\t' '\n')
for c in GID GENOME_CODE SPECIES GENOME_FA DANTE_GFF EDTA_TABLE EDTA_DIR STRAND_COL ZERO_BASED MIU; do
  printf '%s\n' "${MCOLS[@]}" | grep -qx "$c" || { echo "Manifest is missing column $c"; exit 1; }
done
genome_rows() {   # prints manifest rows (no header) for the selected genomes, as key=value lines per genome
  awk -F'\t' -v sel=",$GIDS," 'NR==1 {for (i=1;i<=NF;i++) h[i]=$i; next}
      sel == ",," || index(sel, ","$1",") { line=""; for (i=1;i<=NF;i++) line=line h[i]"="$i"\t"; print line }' "$MANIFEST"
}
load_row() {      # sets GID GENOME_CODE ... from one key=value row
  local kv; IFS=$'\t' read -r -a kv <<< "$1"
  for x in "${kv[@]}"; do [[ -n "$x" ]] && printf -v "${x%%=*}" '%s' "${x#*=}"; done
}
N_GENOMES=$(genome_rows | wc -l)
(( N_GENOMES > 0 )) || { echo "No genomes selected (-g $GIDS)"; exit 1; }
log "$N_GENOMES genomes selected"

MASTER_HEADER="Chrom\tStart\tEnd\tLTR_ID\tstrand\tTSD\tTE_Class\tLTR_type\tSuperfamily\tLineage\tFinal_Classification\tRank\tLTR5_length\tLTR5_Start\tLTR5_End\tLTR3_length\tLTR3_Start\tLTR3_End\tSource"

#==============================================================================
# step 1  EDTA LTR table: pass.list (LTR_retriever) + TEsorter cls.tsv
#         columns are found by header name (Strand, Insertion_Time; Superfamily, Clade)
#==============================================================================
if run_step 1; then
  log "step 1: EDTA LTR tables"
  while read -r ROW; do
    load_row "$ROW"
    T="$OUT/01_edta_table/${GID}_EDTA_LTRs.tsv"
    if [[ "$EDTA_TABLE" != "NA" && -s "$EDTA_TABLE" && "$EDTA_DIR" == "NA" ]]; then
      cp "$EDTA_TABLE" "$T"; log "  $GID: using existing EDTA table $EDTA_TABLE"; continue
    fi
    [[ -d "$EDTA_DIR" ]] || { log "  WARNING $GID: no EDTA_DIR and no EDTA_TABLE -- skipped"; continue; }
    shopt -s nullglob; PASS=("$EDTA_DIR"/*.pass.list); CLS=("$EDTA_DIR"/*.cls.tsv); shopt -u nullglob
    (( ${#PASS[@]} && ${#CLS[@]} )) || { log "  WARNING $GID: pass.list/cls.tsv not found in $EDTA_DIR"; continue; }
    awk -F'\t' -v OFS='\t' -v gid="$GID" -v gcode="$GENOME_CODE" -v sp="$SPECIES" '
      function loc_parse(s, r,   i, j) {           # chr:start..end (chr may contain dots/colons)
        sub(/#.*/, "", s); sub(/\|.*/, "", s)
        i = match(s, /:[0-9]+\.\.[0-9]+$/); if (!i) return 0
        r["chr"] = substr(s, 1, i - 1); split(substr(s, i + 1), j, /\.\./); r["s"] = j[1]; r["e"] = j[2]; return 1
      }
      FNR == 1 { file++ }
      file <= ncls && /^#/ { for (i=1;i<=NF;i++) { if ($i=="Superfamily") sfc=i; if ($i=="Clade") clc=i }; next }
      file <= ncls { if (loc_parse($1, r)) fam[r["chr"]":"r["s"]".."r["e"]] = $(sfc ? sfc : 3) "_" $(clc ? clc : 4); next }
      /^#/ { for (i=1;i<=NF;i++) { if ($i=="Strand") stc=i; if ($i=="Insertion_Time") agc=i }; next }
      { if (!loc_parse($1, r)) next
        k = r["chr"]":"r["s"]".."r["e"]; n++; if (k in fam) hit++
        age = $(agc ? agc : 12); if (age == "") age = "NA"
        print gid "_" gcode, gid, sp, r["chr"], r["s"], r["e"], $(stc ? stc : 9), (k in fam ? fam[k] : "NA"), age }
      END { printf "  %s: %d EDTA LTRs, %d with TEsorter family\n", gid, n, hit > "/dev/stderr" }
    ' ncls="${#CLS[@]}" "${CLS[@]}" "${PASS[@]}" \
    | sed '1i Genome_code\tG_id\tSpecies\tChrom\tStart\tEnd\tStrand\tLTR_fam\tTE_age_years' > "$T"
  done < <(genome_rows)
fi

#==============================================================================
# step 2  DANTE_LTR master table
#         intact elements = GFF3 type transposable_element (not TE_partial);
#         5'/3' LTR coordinates from child long_terminal_repeat features (LTR=5LTR/3LTR),
#         or from LTR5_length / LTR3_length attributes if those are present instead.
#==============================================================================
if run_step 2; then
  log "step 2: DANTE master tables"
  while read -r ROW; do
    load_row "$ROW"
    G="$DANTE_GFF"
    if [[ -d "$G" ]]; then                       # per-chromosome GFF3 files -> one file
      shopt -s nullglob; parts=("$G"/*.gff3); shopt -u nullglob
      (( ${#parts[@]} )) || { log "  WARNING $GID: no *.gff3 in $G"; continue; }
      grep -hv '^#' "${parts[@]}" > "$OUT/02_dante_table/${GID}.DANTE_merged.gff3"; G="$OUT/02_dante_table/${GID}.DANTE_merged.gff3"
    fi
    [[ -s "$G" ]] || { log "  WARNING $GID: DANTE GFF3 not found ($DANTE_GFF)"; continue; }
    awk -F'\t' -v OFS='\t' -v gid="$GID" '
      function attr(a, k,   n, i, p, kv) { n = split(a, p, ";"); for (i=1;i<=n;i++) if (split(p[i], kv, "=") >= 2 && kv[1] == k) return kv[2]; return "NA" }
      /^#/ { next }
      $3 == "transposable_element" && $9 !~ /TE_partial/ {
        id = attr($9, "ID"); n++; te[n] = id
        chr[id]=$1; s[id]=($4<$5?$4:$5); e[id]=($4<$5?$5:$4); st[id]=$7; at[id]=$9; next }
      $3 == "long_terminal_repeat" {
        p = attr($9, "Parent"); w = attr($9, "LTR")
        if (w == "5LTR") { l5s[p]=$4; l5e[p]=$5 } else if (w == "3LTR") { l3s[p]=$4; l3e[p]=$5 } }
      END {
        for (k=1; k<=n; k++) { id=te[k]; print chr[id], s[id], e[id], id }
      }' "$G" | sort -k1,1 -k2,2n -k3,3n > "$OUT/02_dante_table/${GID}.order.tmp"
    awk -F'\t' -v OFS='\t' -v gid="$GID" -v hdr="$MASTER_HEADER" '
      function attr(a, k,   n, i, p, kv) { n = split(a, p, ";"); for (i=1;i<=n;i++) if (split(p[i], kv, "=") >= 2 && kv[1] == k) return kv[2]; return "NA" }
      FNR == NR { ord[++m] = $4; next }
      /^#/ { next }
      $3 == "transposable_element" && $9 !~ /TE_partial/ {
        id = attr($9, "ID"); chr[id]=$1; s[id]=($4<$5?$4:$5); e[id]=($4<$5?$5:$4); st[id]=$7; at[id]=$9; next }
      $3 == "long_terminal_repeat" {
        p = attr($9, "Parent"); w = attr($9, "LTR")
        if (w == "5LTR") { l5s[p]=$4; l5e[p]=$5 } else if (w == "3LTR") { l3s[p]=$4; l3e[p]=$5 } }
      END {
        gsub(/\\t/, "\t", hdr); print hdr
        for (k=1; k<=m; k++) {
          id = ord[k]; a = at[id]
          fc = attr(a, "Final_Classification"); nt = split(fc, tx, "|")
          sf = "NA"; for (i=1;i<=nt;i++) if (tx[i] ~ /^Ty[13]\//) { sf = tx[i]; sub(/^Ty[13]\//, "", sf); sf = toupper(substr(sf,1,1)) substr(sf,2) }
          lin = (nt ? tx[nt] : "NA")
          if (id in l5s) { a5s=l5s[id]; a5e=l5e[id]; a3s=l3s[id]; a3e=l3e[id] }
          else {                                        # fallback: lengths in attributes
            L5 = attr(a, "LTR5_length"); L3 = attr(a, "LTR3_length")
            if (L5 == "NA" || L3 == "NA") { a5s=a5e=a3s=a3e="NA" }
            else if (st[id] == "-") { a3s=s[id]; a3e=s[id]+L3-1; a5e=e[id]; a5s=e[id]-L5+1 }
            else                    { a5s=s[id]; a5e=s[id]+L5-1; a3e=e[id]; a3s=e[id]-L3+1 }
          }
          l5 = (a5s == "NA" ? "NA" : a5e - a5s + 1); l3 = (a3s == "NA" ? "NA" : a3e - a3s + 1)
          if (a5s == "NA") noltr++
          print chr[id], s[id], e[id], sprintf("%s_LTR_%06d", gid, k), st[id], attr(a, "TSD"),
                (nt>=1?tx[1]:"NA"), (nt>=2?tx[2]:"NA"), sf, lin, fc, attr(a, "Rank"),
                l5, a5s, a5e, l3, a3s, a3e, "DANTE_LTR"
        }
        printf "  %s: %d intact DANTE_LTR elements, %d without LTR coordinates\n", gid, m, noltr > "/dev/stderr"
      }' "$OUT/02_dante_table/${GID}.order.tmp" "$G" > "$OUT/02_dante_table/${GID}.master.tsv"
    rm -f "$OUT/02_dante_table/${GID}.order.tmp"
  done < <(genome_rows)
fi

#==============================================================================
# step 3  DANTE age: 5' vs 3' LTR divergence (helper scripts), joined by LTR_ID
#==============================================================================
if run_step 3; then
  log "step 3: DANTE 5' vs 3' LTR age"
  for h in make_ltr_inputs.py ltr_divergence_dnateV2.py; do
    [[ -s "$HELPERS/$h" ]] || { echo "Missing helper $HELPERS/$h (use -H)"; exit 1; }
  done
  while read -r ROW; do
    load_row "$ROW"
    M="$OUT/02_dante_table/${GID}.master.tsv"; P="$OUT/03_dante_age/$GID/$GID"; mkdir -p "$(dirname "$P")"
    [[ -s "$M" ]] || { log "  WARNING $GID: run step 2 first"; continue; }
    [[ -s "$P.master_with_age.tsv" ]] && continue
    ZB=(); is_true "$ZERO_BASED" && ZB=(--zero_based)
    python3 "$HELPERS/make_ltr_inputs.py" -i "$M" -o "$P" --strand_col "$STRAND_COL" "${ZB[@]}"
    bedtools getfasta -fi "$GENOME_FA" -bed "$P.LTR5.bed" -s -name | sed 's/::.*$//' > "$P.LTR5.fa"
    bedtools getfasta -fi "$GENOME_FA" -bed "$P.LTR3.bed" -s -name | sed 's/::.*$//' > "$P.LTR3.fa"
    python3 "$HELPERS/ltr_divergence_dnateV2.py" --ltr5 "$P.LTR5.fa" --ltr3 "$P.LTR3.fa" --concat_spacer 0 \
        -o "$P.LTRpair_age.tsv" -t "$THREADS" --miu "$MIU" --blast_task "$BLAST_TASK" --word_size "$WORD_SIZE" \
        --evalue "$EVALUE" --end_window "$END_WINDOW" --max_off "$MAX_OFF" --min_align "$MIN_ALIGN" \
        > "$OUT/logs/$GID.dante_age.log" 2>&1
    # join raw_T / K2P_T / JC69_T (found by header name) onto every master row by LTR_ID
    awk -F'\t' -v OFS='\t' '
      FNR == NR { if (FNR == 1) { for (i=1;i<=NF;i++) c[$i]=i; if (!c["raw_T"] || !c["K2P_T"] || !c["JC69_T"]) { print "age table lacks raw_T/K2P_T/JC69_T" > "/dev/stderr"; exit 1 }; next }
                  r[$1]=$(c["raw_T"]); k[$1]=$(c["K2P_T"]); j[$1]=$(c["JC69_T"]); next }
      FNR == 1 { print $0, "raw_T", "K2P_T", "JC69_T"; next }
      { id=$4; print $0, (id in r ? r[id] : "NA"), (id in k ? k[id] : "NA"), (id in j ? j[id] : "NA"); n++; if (id in k) a++ }
      END { printf "  %s: age for %d of %d DANTE elements\n", gid, a, n > "/dev/stderr" }' gid="$GID" \
      "$P.LTRpair_age.tsv" "$M" > "$P.master_with_age.tsv"
  done < <(genome_rows)
fi

#==============================================================================
# step 4  EDTA-unique elements (< 50% overlap with DANTE), LTR positions, age
#         EDTA age (LTR_retriever, years) is written to raw_T, K2P_T and JC69_T alike.
#==============================================================================
if run_step 4; then
  log "step 4: EDTA-unique LTRs"
  [[ -s "$HELPERS/edta_infer_ltr_positions.py" ]] || { echo "Missing helper $HELPERS/edta_infer_ltr_positions.py (use -H)"; exit 1; }
  while read -r ROW; do
    load_row "$ROW"
    D="$OUT/03_dante_age/$GID/$GID.master_with_age.tsv"; E="$OUT/01_edta_table/${GID}_EDTA_LTRs.tsv"
    W="$OUT/04_edta_unique/$GID"; mkdir -p "$W"
    [[ -s "$D" && -s "$E" ]] || { log "  WARNING $GID: needs steps 1 and 3 -- skipped"; continue; }
    tail -n +2 "$D" | cut -f1-4 | sort -k1,1 -k2,2n > "$W/dante.4col.bed"
    tail -n +2 "$E" | sort -k4,4 -k5,5n -k6,6n | awk -F'\t' -v OFS='\t' -v gid="$GID" -v b4="$W/edta.4col.bed" '
      { id = sprintf("%s_ed_%06d", gid, NR); print $4, $5, $6, id > b4; print $4, $5, $6, id, $7, $8 "#" $9 }' > "$W/edta.6col.bed"
    bedtools intersect -a "$W/edta.4col.bed" -b "$W/dante.4col.bed" -v -f "$EDTA_DANTE_OVERLAP" > "$W/edta.uniq.4col.bed"
    awk 'NR==FNR {k[$4]; next} ($4 in k)' "$W/edta.uniq.4col.bed" "$W/edta.6col.bed" > "$W/edta.uniq.6col.bed"
    log "  $GID: $(wc -l < "$W/edta.uniq.4col.bed") of $(wc -l < "$W/edta.4col.bed") EDTA LTRs are EDTA-unique"
    if [[ -s "$W/edta.uniq.6col.bed" ]]; then
      python3 "$HELPERS/edta_infer_ltr_positions.py" -b "$W/edta.uniq.6col.bed" -g "$GENOME_FA" -o "$W/$GID" \
          --window "$INFER_WINDOW" --min_align "$INFER_MIN_ALIGN" --max_off "$INFER_MAX_OFF" \
          --pident_min "$INFER_PIDENT_MIN" --threads "$THREADS" > "$OUT/logs/$GID.edta_infer.log" 2>&1
    fi
    awk -F'\t' -v OFS='\t' -v hdr="$MASTER_HEADER" '
      BEGIN { gsub(/\\t/, "\t", hdr); print hdr, "raw_T", "K2P_T", "JC69_T" }
      FNR == NR { if (FNR > 1) { a[$1]=$2; b[$1]=$3; c[$1]=$4; d[$1]=$5 }; next }
      { split($6, nf, "#"); fam = nf[1]; age = nf[2]; split(fam, sc, "_")
        if ($4 in a) { s5=a[$4]; e5=b[$4]; s3=c[$4]; e3=d[$4]; l5=e5-s5+1; l3=e3-s3+1 } else { s5=e5=s3=e3=l5=l3="NA" }
        print $1, $2, $3, $4, $5, "NA", "NA", "NA", (sc[1] != "" ? sc[1] : "NA"), (sc[2] != "" ? sc[2] : "NA"), fam, "NA",
              l5, s5, e5, l3, s3, e3, "EDTA", age, age, age }' \
      <( [[ -s "$W/$GID.ltr_ends.tsv" ]] && cat "$W/$GID.ltr_ends.tsv" || echo "id" ) "$W/edta.uniq.6col.bed" \
      > "$W/$GID.edta_unique_with_age.tsv"
  done < <(genome_rows)
fi

#==============================================================================
# step 5  Final tables: DANTE + EDTA-unique, age in years -> Mya (once), age group
#==============================================================================
if run_step 5; then
  log "step 5: final LTR age tables (age column: $AGE_COLUMN)"
  FINAL="$OUT/05_final/Pisum_all_LTR_age.tsv"
  HDR="Genome_code\tG_id\tSpecies\tLTR_loc\tLTR_start\tLTR_end\tStrand\tSuperfamily\tLineage\tLTR_fam\tLTR_ageMya\tLTR_ageYears\tTE_ageGroup\tLTR_ID\tSource\tTSD\tRank\tLTR5_Start\tLTR5_End\tLTR3_Start\tLTR3_End\tFinal_Classification"
  printf "%b\n" "$HDR" > "$FINAL.tmp"
  while read -r ROW; do
    load_row "$ROW"
    D="$OUT/03_dante_age/$GID/$GID.master_with_age.tsv"; E="$OUT/04_edta_unique/$GID/$GID.edta_unique_with_age.tsv"
    [[ -s "$D" || -s "$E" ]] || { log "  WARNING $GID: no age tables -- skipped"; continue; }
    { [[ -s "$D" ]] && tail -n +2 "$D"; [[ -s "$E" ]] && tail -n +2 "$E"; } | sort -k1,1 -k2,2n -k3,3n \
    | awk -F'\t' -v OFS='\t' -v gid="$GID" -v gcode="$GENOME_CODE" -v sp="$SPECIES" -v agecol="$AGE_COLUMN" '
        function grp(m) { if (m < 0.25) return "<0.25"; if (m < 0.5) return "0.25-0.5"; if (m < 1) return "0.5-1"
                          if (m < 2) return "1-2"; if (m < 3) return "2-3"; if (m < 4) return "3-4"; return ">4" }
        BEGIN { col = (agecol == "raw_T" ? 20 : agecol == "JC69_T" ? 22 : 21) }
        { y = $col; if (y == "" || y == "NA") { skip[$19]++; next }
          fam = ($9 != "NA" && $10 != "NA") ? $9 "_" $10 : $11
          loc = (index($1, gcode ".") == 1 ? $1 : gcode "." $1)
          print gid "_" gcode, gid, sp, loc, $2, $3, $5, $9, $10, fam, y / 1e6, y, grp(y / 1e6),
                $4, $19, $6, $12, $14, $15, $17, $18, $11; keep[$19]++ }
        END { printf "  %s: kept DANTE %d, EDTA %d; no age: DANTE %d, EDTA %d\n", gid,
                     keep["DANTE_LTR"], keep["EDTA"], skip["DANTE_LTR"], skip["EDTA"] > "/dev/stderr" }' \
    > "$OUT/05_final/$GID.body"
    { printf "%b\n" "$HDR"; cat "$OUT/05_final/$GID.body"; } > "$OUT/05_final/$GID.master_LTR_age.tsv"
    cat "$OUT/05_final/$GID.body" >> "$FINAL.tmp"; rm -f "$OUT/05_final/$GID.body"
  done < <(genome_rows)
  { head -n1 "$FINAL.tmp"; tail -n +2 "$FINAL.tmp" | sort -t$'\t' -k1,1 -k4,4 -k5,5n; } > "$FINAL"; rm -f "$FINAL.tmp"
  # summary: elements by genome, source and age group
  awk -F'\t' 'NR > 1 { n[$2"\t"$15"\t"$13]++ } END { for (k in n) print k"\t"n[k] }' "$FINAL" \
    | sort -k1,1V -k2,2 -k3,3 | sed '1i G_id\tSource\tTE_ageGroup\tn' > "$OUT/05_final/LTR_age_summary.tsv"
  log "  $(($(wc -l < "$FINAL") - 1)) LTRs -> $FINAL"
fi

log "finished steps: $STEPS"