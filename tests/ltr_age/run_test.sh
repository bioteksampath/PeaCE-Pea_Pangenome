#!/usr/bin/env bash
# End-to-end test of ltr_age/ltr_age_pipeline.sh on a simulated genome with known answers.
# Needs: python3 + Biopython, blastn, bedtools.   Usage: bash tests/ltr_age/run_test.sh
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT

python3 "$HERE/simulate_ltr_genome.py" "$W/sim"
bash "$ROOT/ltr_age/ltr_age_pipeline.sh" -m "$W/sim/manifest.tsv" -o "$W/out" -t 2 > "$W/log.txt" 2>&1 \
  || { cat "$W/log.txt"; echo "FAIL: pipeline exited with an error"; exit 1; }

F="$W/out/05_final/Pisum_all_LTR_age.tsv"
fail=0
check() { if eval "$2"; then echo "PASS  $1"; else echo "FAIL  $1"; fail=1; fi; }
row() { awk -F'\t' -v id="$1" '$14 == id' "$F"; }

check "3 elements in final table (partial and overlapping EDTA call removed)" \
      '[[ $(tail -n +2 "$F" | wc -l) -eq 3 ]]'
check "TE_1 (2% divergence): age 0.76-0.80 Mya, group 0.5-1" \
      'row G1_LTR_000001 | awk -F"\t" '"'"'{exit !($11 > 0.76 && $11 < 0.80 && $13 == "0.5-1")}'"'"''
check "TE_2 (identical LTRs): age 0" \
      'row G1_LTR_000002 | awk -F"\t" '"'"'{exit !($11 == 0)}'"'"''
check "EDTA-only element keeps LTR_retriever age 1.15 Mya" \
      'tail -n +2 "$F" | awk -F"\t" '"'"'$15=="EDTA" && $11==1.15 {f=1} END{exit !f}'"'"''
check "EDTA-only element: 5'/3' LTR coordinates within 10 bp (minus strand)" \
      'paste <(awk -F"\t" '"'"'$15=="EDTA" {print $18"\t"$19"\t"$20"\t"$21}'"'"' "$F") "$W/sim/expected_edta_ltrs.txt" \
         | awk '"'"'{for (i=1;i<=4;i++) { d=$i-$(i+4); if (d<0) d=-d; if (d>10) exit 1 }}'"'"''
check "lineage parsed from DANTE classification (Gypsy_Ogre, Gypsy_Tekay)" \
      '[[ $(cut -f10 "$F" | grep -c -E "^Gypsy_(Ogre|Tekay)$") -eq 2 ]]'

(( fail == 0 )) && echo "All LTR pipeline tests passed." || { echo "Some tests failed."; exit 1; }
