# ltr_age

Full-length LTR retrotransposons (FL-LTR-RTs) from DANTE_LTR and EDTA, and their insertion ages (Fig. 3a LTR-RT density; Extended Data Fig. 2b, 2d; Supplementary Fig. 3).

| Script | Purpose |
|---|---|
| `build_ltr_manifest.py` | Builds the genome manifest (one row per genome) from per-batch map files |
| `ltr_age_pipeline.sh` | Steps 1–5 below |
| `lib/make_ltr_inputs.py` | Step 3: 5′ and 3′ LTR BED files from the DANTE master table |
| `lib/ltr_divergence_dnateV2.py` | Step 3: 5′ vs 3′ LTR BLASTN divergence; p, K2P, JC69 distances and ages |
| `lib/edta_infer_ltr_positions.py` | Step 4: 5′/3′ LTR ends of EDTA-only elements |

```bash
python build_ltr_manifest.py --map1 map1_*.map --map2 map2_*.txt \
       --edta-dir-template /path/EDTA/{GENOME_CODE} --out genome_manifest.tsv
bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -t 32              # all genomes
bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -g G4 -p 1,2,3,4   # one genome (e.g. SLURM array task)
bash ltr_age_pipeline.sh -m genome_manifest.tsv -o ltr_age -p 5                # final merge, all genomes
```

**Manifest columns** (read by name): `GID GENOME_CODE SPECIES GENOME_FA DANTE_GFF EDTA_TABLE EDTA_DIR STRAND_COL ZERO_BASED MIU`.

| Step | What |
|---|---|
| 1 | EDTA table: LTR_retriever `*.pass.list` + TEsorter `*.cls.tsv` (rexdb-plant) → position, strand, lineage, age (years) |
| 2 | DANTE table: intact `transposable_element` features (not `TE_partial`) from the DANTE_LTR GFF3, unique IDs, 5′/3′ LTR coordinates |
| 3 | DANTE age: LTRs extracted with `bedtools getfasta -s`, aligned with BLASTN (`-word_size 7 -evalue 1e-5`; HSP ≥ 50 bp, within 300 bp of both LTR starts, diagonal offset ≤ 120 bp); K2P distance; T = K / (2μ), μ = 1.3 × 10⁻⁸ |
| 4 | EDTA-unique elements (< 50% overlap with DANTE) added with their LTR_retriever age; LTR ends from aligning the first and last 4 kb (HSP ≥ 200 bp, ≥ 75% identity, within 80 bp of the ends) |
| 5 | Final tables per genome and for all genomes: age in years and Mya, age group `<0.25, 0.25-0.5, 0.5-1, 1-2, 2-3, 3-4, >4` Mya |

Output: `ltr_age/05_final/Pisum_all_LTR_age.tsv`, `ltr_age/05_final/LTR_age_summary.tsv`. Age column: `AGE_COLUMN` in the script (default `K2P_T`; `raw_T` and `JC69_T` also available).

Test: `bash tests/ltr_age/run_test.sh` (simulated genome with known ages).
