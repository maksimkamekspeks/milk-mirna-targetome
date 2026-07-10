# Milk Exosomal microRNA Targetome — Analysis Pipeline

Code to reproduce the analyses and figures of:

> **Conserved Core and Species-Specific Signatures in the Milk Exosomal microRNA Targetome: A Comparative In Silico Analysis of Human, Cow, Goat and Donkey Milk.**

This repository contains the **analysis code** plus a small `data/` folder of
**intermediate files** that reproduce every figure without the raw reads. The
raw sequencing data and the two target databases are large and/or third-party
and must be obtained separately only if you want to re-run the full pipeline
from scratch (see **Required external data** below). All results are
computational predictions intended as a hypothesis-generating resource.

---

## Repository layout

```
milk-mirna-targetome/
├── README.md                 # this file
├── LICENSE                   # MIT
├── CITATION.cff              # how to cite (software + article)
├── requirements.txt          # pip dependencies
├── environment.yml           # conda environment (analysis + Step 0 tools)
├── run_all.sh                # reproduce all figures from data/ in one command
├── data/                     # intermediate files -> reproduce figures WITHOUT raw reads
│   ├── README.md
│   ├── human/  cow/  goat/  donkey/    # per-sample miRNA spectra (CSV)
│   ├── FINAL_Integrated_miRNA_Targets_Report.xlsx
│   ├── Categorized_Genes.xlsx
│   ├── miRNA_gene_edges.csv.gz
│   └── combined_miRNA.xlsx
└── scripts/
    ├── 00_prepare_reads.sh             # Step 0: SRA -> cutadapt -> bowtie(miRBase) -> SAM
    ├── 01_build_count_matrix.py        # SAM alignments -> per-species miRNA count matrix + TOP-90%
    ├── 02_combine_spectra.py           # merge per-species spectra (Supplementary S3)
    ├── 03_build_integrated_targetome.py# mirDIP x abundance x TarBase -> scores + edge list
    ├── 04_categorize_genes.py          # cross-species gene categorization (Supplementary S4)
    ├── 05_select_top_hub_genes.py      # top genes per category (Supplementary S5)
    ├── 06_plot_regulatory_power.py     # category-mass figure
    ├── 07_plot_heatmap.py              # signature heatmap (Figure 3)
    ├── 08_enrichment_by_category.py    # weighted GO/KEGG per category (Supplementary S7)
    ├── 09_enrichment_by_milk.py        # weighted GO/KEGG per milk (Supplementary S8)
    ├── 10_robustness_permutation.py    # sensitivity + permutation + correlation (Figure 5, S9)
    ├── 11_sample_clustering_pca.py     # PCA + hierarchical clustering (Figure 6)
    └── 12_mirna_gene_network.py        # conserved-miRNA / core-gene network (Figure 7)
```

All top-level files sit in the repository root; only `data/` and `scripts/` are
folders.

## Installation

**Option A — conda (recommended; includes the Step 0 tools):**
```bash
conda env create -f environment.yml
conda activate milk-mirna-targetome
```

**Option B — pip (analysis scripts only; install sra-tools/cutadapt/bowtie separately for Step 0):**
```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```
Python 3.9+ recommended. Scripts 08–09 query the Enrichr web service and need
internet access.

---

## Required external data (NOT included in this archive)

Place these in the working directory before running the pipeline.

**1. Small-RNA sequencing reads (public; NCBI SRA / GEO).**
Download the FASTQ files for the accessions below, then trim adapters and align
each sample to the human miRBase reference (Step 0). Put the resulting SAM
files in one folder per species, named `clean_<SRR>.sam`:

| Species | Accessions (BioProject / GEO) |
|---------|-------------------------------|
| Human   | PRJNA147351 / GSE32253 ; PRJNA305166 / GSE75726 |
| Cow     | PRJNA771627 ; PRJNA1165141 |
| Goat    | PRJNA771627 ; PRJNA1130716 |
| Donkey  | PRJNA771627 ; PRJNA1130716 |

Expected folders: `human/`, `cow/`, `goat/`, `donkey/` (each with `*.sam`).
Per-sample provenance (milk type, lactation stage, breed, platform, EV
isolation) is given in Supplementary Table S10 of the paper.

**2. mirDIP unidirectional export (predicted targets).**
Obtain from https://ormbase.org/ (mirDIP v5.2), unidirectional search, and save as:
```
mirDIP_Unidirectional_search_v_5_2/mirDIP_Unidirectional_search_v.5.2.txt
```
Columns used: `GENE_SYMBOL, MICRORNA, INTEGRATED_RANK, SCORE_CLASS` (only
classes **V** and **H** are kept).

**3. TarBase v9 (experimentally validated interactions).**
Obtain the human file from https://dianalab.e-ce.uth.gr/tarbasev9 and save as:
```
Homo_sapiens.tsv
```
Columns used: `mirna_name, gene_name`.

---

## Step 0 — read preparation

`scripts/00_prepare_reads.sh` provides a ready example: download FASTQ from SRA
(`sra-tools`), trim the 3' adapter and length-filter (`cutadapt`), and align to
the human miRBase **mature** miRNAs (`bowtie`), writing `clean_<SRR>.sam` into
the species folder. Build the miRBase bowtie index once (commands are inside the
script), set `SPECIES`, `SRR_LIST` and `ADAPTER` per dataset, then:
```bash
bash scripts/00_prepare_reads.sh
```
**Note:** the 3' adapter differs between datasets (e.g. Illumina TruSeq small
RNA `TGGAATTCTCGGGTGCCAAGG` vs NEBNext `AGATCGGAAGAGCACACGTCT`); set it for each
dataset before running.

## Quick start — reproduce figures without raw reads

A `data/` folder is provided with the intermediate files (spectra, integrated
targetome, categorised genes, edge list). To reproduce all downstream results
and figures **without** downloading the ~30 GB of reads or the databases:
```bash
cp -r data/* .
python scripts/04_categorize_genes.py
python scripts/07_plot_heatmap.py            # Figure 3
python scripts/10_robustness_permutation.py  # Figure 5 (+ S9)
python scripts/11_sample_clustering_pca.py   # Figure 6
python scripts/12_mirna_gene_network.py      # Figure 7
# 05, 06 also run from here; 08, 09 additionally need internet (Enrichr)
```
See `data/README.md` for details.

## How to run from raw data (end to end)

```bash
python scripts/01_build_count_matrix.py          # -> <sp>/final_miRNA_CountMatrix_<sp>_full.csv and _TOP90.csv
python scripts/02_combine_spectra.py             # -> combined_miRNA.xlsx
python scripts/03_build_integrated_targetome.py  # -> FINAL_Integrated_miRNA_Targets_Report.xlsx, miRNA_gene_edges.csv.gz
python scripts/04_categorize_genes.py            # -> Categorized_Genes.xlsx
# (re-run 03 once more if you want the edge list restricted to categorised genes)
python scripts/05_select_top_hub_genes.py        # -> Top_Hub_Genes_For_Discussion.xlsx
python scripts/06_plot_regulatory_power.py       # -> Regulatory_Power_Distribution.png
python scripts/07_plot_heatmap.py                # -> Heatmap_Strict_Categories_Ranked_Scores.png
python scripts/08_enrichment_by_category.py      # -> Enrichment_Plots/*.png        (internet)
python scripts/09_enrichment_by_milk.py          # -> Enrichment_By_Milk/*.png       (internet)
python scripts/10_robustness_permutation.py      # -> Figure5_robustness.png, S9_robustness_analyses.xlsx
python scripts/11_sample_clustering_pca.py       # -> Figure6_sample_clustering.png
python scripts/12_mirna_gene_network.py          # -> Figure7_miRNA_gene_network.png
```

## Method summary

For each gene *g* and species *s* the weighted targeting score is
`score(g,s) = Σ_m a(m,s)·(r(m,g) + b(m,g))`, where *a* is the mean relative
exosomal abundance of miRNA *m* in *s*, *r* is the mirDIP integrated rank
(V/H classes) and *b* = 0.5 if the (*m*,*g*) pair is validated in TarBase v9.
Genes are filtered (noise cut-off) and assigned to cross-species categories by
a fold-change gap rule (defaults: 60th-percentile noise threshold, fold change
1.5, minimum summed score 10).

## Parameters you can change

`FOLD_CHANGE`, `PERCENTILE_THRESHOLD`, `MIN_SUM_SCORE` (in
`04_categorize_genes.py`) and `VALIDATION_BONUS` (in
`03_build_integrated_targetome.py`). The threshold-sensitivity analysis in
`10_robustness_permutation.py` sweeps these automatically.

## License

Released under the MIT License (add a LICENSE file before publishing).

## Citation

If you use this code, please cite the associated article.
