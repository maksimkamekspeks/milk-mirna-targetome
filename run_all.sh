#!/usr/bin/env bash
# Reproduce all downstream results and figures from the intermediate files in
# data/ (no raw reads and no databases required). Run from the repository root.
# Scripts 08/09 (functional enrichment) additionally need internet (Enrichr) and
# are left commented out.
set -euo pipefail

cp -r data/* .

python scripts/04_categorize_genes.py        # -> Categorized_Genes.xlsx
python scripts/05_select_top_hub_genes.py     # -> Top_Hub_Genes_For_Discussion.xlsx
python scripts/06_plot_regulatory_power.py    # -> Regulatory_Power_Distribution.png
python scripts/07_plot_heatmap.py             # -> Heatmap_Strict_Categories_Ranked_Scores.png (Figure 3)
python scripts/10_robustness_permutation.py   # -> Figure5_robustness.png + S9_robustness_analyses.xlsx
python scripts/11_sample_clustering_pca.py    # -> Figure6_sample_clustering.png
python scripts/12_mirna_gene_network.py       # -> Figure7_miRNA_gene_network.png
# python scripts/08_enrichment_by_category.py # needs internet (Enrichr)
# python scripts/09_enrichment_by_milk.py     # needs internet (Enrichr)

echo "Done. See Figure3/5/6/7, Regulatory_Power_Distribution.png and S9_robustness_analyses.xlsx"
