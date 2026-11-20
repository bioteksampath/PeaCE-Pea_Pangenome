# pan_rgaome

Pan-RGAome analyses from MMseqs2 clusters of filtered resistance gene analogues (Fig. 2b–d).

| Script | Purpose | Figure |
|---|---|---|
| `rga_core_shell_cloud.R` | Core/shell/cloud classes, RGA class per orthogroup, Monte Carlo chi-square, *P. sativum*-specific vs shared orthogroups | Fig. 2c, 2d |
| `rga_saturation.R` | Orthogroup accumulation over 500 random assembly orders | Fig. 2b |

```bash
Rscript rga_core_shell_cloud.R --clusters pisum85_cluster.tsv --rga pisum85_RGA.tsv --taxa pisum85_taxa.tsv --out Pisum85_RGA
Rscript rga_saturation.R --clusters pisum85_cluster.tsv --perm 500 --out Pisum85_RGA
```

- **Inputs.** MMseqs2 cluster table (`sample gene_id gene_id cluster_id cluster_size`, no header); RGA annotation (`gene_id length RGA_class`, RGAugury classes after NBS and DeepTMHMM filtering); optional `sample<TAB>taxon`.
- **Categories.** Orthogroups in ≥ 3 assemblies; cloud 3–18, shell 19–63, core 64–85. An orthogroup takes an RGA class if ≥ 50% of its annotated genes have it.
- **Test.** TNL/CNL/RNL/RLK/RLP × category, Monte Carlo chi-square (B = 100,000, seed 42).

Test: `bash tests/pan_rgaome/run_test.sh`.
