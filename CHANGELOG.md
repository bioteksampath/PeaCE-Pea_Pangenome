# Changelog

## [Unreleased]
- Analyses to be added: assembly, gene annotation, orthology, pangenome graph, synteny, TIP genotyping (Paragraph), TIP-eQTL, PanG2P GWAS, AES1 haplotypes.

## [1.0.0] - 2026-11-09
First version for the manuscript *Pea pangenomics connects structural variation with disease*.
- `snp_population`: short-read SNP calling, filtering, imputation, PCA (Fig. 1a)
- `ltr_age`: full-length LTR-RTs (DANTE_LTR + EDTA) and insertion ages (Fig. 3a; Ext. Data Fig. 2b, 2d; Supp. Fig. 3)
- `pan_rgaome`: core/shell/cloud classes and saturation (Fig. 2b-d)
- `sv_longread`: HiFi read-based SVs against CDC Amarillo; pan-SV catalogue (Ext. Data Fig. 5)
- `tip_discovery`: TIPs, lineages, LTR ages, force-calling (Fig. 3; Ext. Data Fig. 6)
- `pm_rnaseq`: powdery mildew RNA-seq and Er3 candidates (Fig. 4f)
- Tests on simulated data; GitHub Actions checks
