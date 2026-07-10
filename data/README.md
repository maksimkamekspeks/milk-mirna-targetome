# Provided intermediate data — reproduce without raw reads or databases

These files are the intermediate outputs of Steps 0–3 of the pipeline. They let
you reproduce **all downstream results and every figure without downloading the
~30 GB of raw sequencing reads and without the mirDIP / TarBase databases.**

## Contents

| File | Produced by | Lets you run |
|------|-------------|--------------|
| `<species>/final_miRNA_CountMatrix_<species>_TOP90.csv` and `_full.csv` | 01 | 02, 03, 11 |
| `combined_miRNA.xlsx` | 02 | — (Supplementary S3) |
| `FINAL_Integrated_miRNA_Targets_Report.xlsx` (sheet `summary` + per-species) | 03 | 04, 10 |
| `miRNA_gene_edges.csv.gz` | 03 | 12 |
| `Categorized_Genes.xlsx` | 04 | 05, 06, 07, 08, 09, 12 |

The per-species CSVs contain one column per sample plus `Avg_Percent` and
`Cumulative_Percent`; they are the top-90 % abundance spectra used in the study.
Sample provenance (which SRR belongs to which dataset/species, milk type,
lactation, breed, platform) is in Supplementary Table S10 of the paper.

## Quick start (no raw reads needed)

From the repository root:
```bash
cp -r data/* .          # bring the intermediate files into the working directory
python scripts/04_categorize_genes.py        # -> Categorized_Genes.xlsx  (Core 1809, Donkey 770, Human 430, Cow-Goat 296)
python scripts/05_select_top_hub_genes.py
python scripts/06_plot_regulatory_power.py
python scripts/07_plot_heatmap.py            # Figure 3
python scripts/10_robustness_permutation.py  # Figure 5 + S9  (silhouette/permutation)
python scripts/11_sample_clustering_pca.py   # Figure 6  (silhouette 0.82)
python scripts/12_mirna_gene_network.py      # Figure 7
# scripts 08 and 09 additionally need internet (Enrichr)
```

To reproduce Steps 0–3 themselves (from raw reads), see the main README and
download the reads (SRA/GEO) and databases (mirDIP, TarBase) listed there.
