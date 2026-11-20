#!/usr/bin/env bash
# Smoke test of pan_rgaome/*.R on simulated MMseqs2 clusters (85 assemblies, 400 orthogroups).
# Needs: python3, R with data.table and ggplot2.   Usage: bash tests/pan_rgaome/run_test.sh
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
cd "$W"

python3 - <<'EOF'
import random
random.seed(1)
classes = ["TNL", "CNL", "RNL", "RLK", "RLP", "TM-CC"]
samples = [f"G{i}" for i in range(1, 86)]
with open("clusters.tsv", "w") as c, open("rga.tsv", "w") as r, open("taxa.tsv", "w") as t:
    for s in samples:
        t.write(f"{s}\t{'P. sativum' if int(s[1:]) < 70 else 'P. fulvum'}\n")
    g = 0
    for og in range(400):
        n = random.choice([random.randint(1, 18), random.randint(19, 63), random.randint(64, 85)])
        k = random.choice(classes)
        for s in random.sample(samples, n):
            g += 1
            c.write(f"{s}\t{s}_g{g}\t{s}_g{g}\tOG{og}\t{n}\n")
            r.write(f"{s}_g{g}\t900\t{k}\n")
EOF

Rscript "$ROOT/pan_rgaome/rga_core_shell_cloud.R" --clusters clusters.tsv --rga rga.tsv --taxa taxa.tsv --B 2000 --out t
Rscript "$ROOT/pan_rgaome/rga_saturation.R" --clusters clusters.tsv --perm 20 --out t
for f in t_OG_summary.tsv t_category_by_class.tsv t_taxon_sharing.tsv t_prevalence_density.pdf \
         t_taxon_sharing.pdf t_saturation_summary.tsv t_saturation_thin.pdf; do
  [[ -s "$f" ]] || { echo "FAIL  missing $f"; exit 1; }
done
# last saturation step must recover 100% of orthogroups
awk -F'\t' 'NR > 1 {last = $2} END {exit !(last == 100)}' t_saturation_summary.tsv \
  || { echo "FAIL  saturation does not reach 100%"; exit 1; }
echo "All pan-RGAome tests passed."
