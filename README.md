# PeaCE-Pea_Pangenome

**Analysis code for the 85-genome pea (*Pisum*) pangenome**

[![checks](https://github.com/bioteksampath/PeaCE-Pea_Pangenome/actions/workflows/checks.yml/badge.svg)](https://github.com/bioteksampath/PeaCE-Pea_Pangenome/actions/workflows/checks.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Custom scripts for the paper:

> **Pea pangenomics connects structural variation with disease** (2026).
> Perumal S. *et al.* Full author list and journal reference will be added on publication.

Raw data: NCBI BioProject [PRJNA1495898](https://www.ncbi.nlm.nih.gov/bioproject/1495898).

Each folder holds one analysis. Each pipeline runs from start to end and writes the tables and plots used in the figures. Analyses that used published tools only are described, with their parameters, in the Methods of the paper.

---

## Analyses

| Folder | Analysis | Main tools | Figures |
|---|---|---|---|
| [`snp_population`](snp_population) | Short-read SNP calling, filtering, imputation, PCA | cutadapt, BWA-MEM, BCFtools, Beagle, PLINK2 | Fig. 1a |
| [`ltr_age`](ltr_age) | Full-length LTR retrotransposons (DANTE_LTR + EDTA) and insertion ages | BLAST+, BEDTools, Python | Fig. 3a; Ext. Data Fig. 2b, 2d; Supp. Fig. 3 |
| [`pan_rgaome`](pan_rgaome) | Pan-RGAome: core/shell/cloud classes, class test, saturation | R (data.table, ggplot2) | Fig. 2b–d |
| [`sv_longread`](sv_longread) | Read-based SVs: HiFi reads mapped to CDC Amarillo; pan-SV catalogue | minimap2, Sniffles2, SVIM, SURVIVOR | Ext. Data Fig. 5a–d |
| [`tip_discovery`](tip_discovery) | Transposon insertion polymorphisms: discovery, lineage, LTR age, force-calling | BLAST+, TEsorter, MAFFT, Sniffles2 | Fig. 3a–e; Ext. Data Fig. 6a–f; Supp. Table 6 |
| [`pm_rnaseq`](pm_rnaseq) | Powdery mildew RNA-seq and Er3 candidate genes | STAR, RSEM, DESeq2 | Fig. 4f |

Every folder has its own README with the run command, steps and settings. [`docs/script_index.tsv`](docs/script_index.tsv) lists every script with its Methods section, figures, inputs and outputs.

```
sequencing reads ─► snp_population ─► SNP VCF, PCA

genome + EDTA + DANTE_LTR ─► ltr_age ─► LTR age table

RGAugury → filtered RGAs → MMseqs2 clusters ─► pan_rgaome

HiFi reads ─► sv_longread ─► pan-SV VCF ─► tip_discovery ─► TIP VCF + ages
                                               ▲
                              pan-TE library ──┘

RNA-seq reads ─► pm_rnaseq ─► DE tables, Er3 candidates, Fig. 4f
```

## Getting started

```bash
git clone https://github.com/bioteksampath/PeaCE-Pea_Pangenome.git
cd PeaCE-Pea_Pangenome
conda env create -f envs/ltr_age.yml      # one environment per analysis, see envs/
conda activate pisum_ltr_age
bash ltr_age/ltr_age_pipeline.sh          # no arguments: prints all options
bash tests/ltr_age/run_test.sh            # test on a small simulated genome
```

| Environment | Used by | Tools (versions from the Methods) |
|---|---|---|
| `envs/snp.yml` | snp_population | cutadapt 5.1, BWA 0.7.19, SAMtools, BCFtools 1.22, Beagle 5.4, PLINK 2.0 |
| `envs/ltr_age.yml` | ltr_age | BLAST+, BEDTools, Python 3, Biopython |
| `envs/sv.yml` | sv_longread | minimap2 2.28, Sniffles2 2.6.2, SVIM 2.0.0, SURVIVOR 1.0.7 |
| `envs/tip.yml` | tip_discovery | BLAST+ 2.12.0, TEsorter 1.4.7, MAFFT, seqkit, BEDTools, Sniffles2 2.7.2 |
| `envs/rnaseq.yml` | pm_rnaseq | FastQC 0.12.1, cutadapt 5.1, STAR 2.7.11b, RSEM 1.3.3 |
| `envs/r_analysis.yml` | pan_rgaome, pm_rnaseq | R ≥ 4.3, data.table, ggplot2, DESeq2, tximport, tidyverse, pheatmap, ggrepel, openxlsx |

## How the pipelines behave

- **Help.** Run a Bash pipeline with no arguments (or `-h`) to see every option and its default. For R scripts, see the header (`head -40 script.R`).
- **Steps.** Pipelines are split into steps; `-p 2,3` runs only steps 2 and 3.
- **Restart.** Finished samples, regions or chromosomes are skipped when you run the same command again.
- **Logs.** Everything printed is also saved in `<outdir>/logs/`.
- **Sample sheets.** Tab-separated. A header line starting with `sample`, blank lines and `#` lines are ignored.
- **Scheduler.** Plain Bash, Python and R; wrap them in your own SLURM (or other) job script.

## Tests

| Test | What it checks |
|---|---|
| `tests/ltr_age/run_test.sh` | Simulated genome with known answers. 2% LTR divergence gives 0.78 Mya; identical LTRs give 0; partial elements and EDTA calls overlapping DANTE are removed; EDTA ages and LTR ends are kept correctly. |
| `tests/pan_rgaome/run_test.sh` | Simulated clusters for 85 assemblies: all outputs are written, and saturation reaches 100%. |

GitHub Actions runs these tests, plus syntax checks of all Bash, Python and R scripts, on every push ([`.github/workflows/checks.yml`](.github/workflows/checks.yml)).

## Related code

- panTEvo (solo-LTR detection, TIP annotation): <https://github.com/Lcamdom/panTEvo>
- PanG2P (pangenome-informed genotyping and association): <https://github.com/sivasubramanics/PanG2P>

## Citation

Please cite the paper and this repository ([`CITATION.cff`](CITATION.cff); GitHub shows a "Cite this repository" button). A Zenodo DOI will be added with the first release.

## Licence

MIT, see [`LICENSE`](LICENSE).

## Contact

Sampath Perumal, National Research Council Canada, Saskatoon, Canada. Please open a GitHub issue for questions or bugs.
