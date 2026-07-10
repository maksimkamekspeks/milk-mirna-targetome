"""
03_build_integrated_targetome.py
================================
Build the abundance-weighted, multi-evidence targetome: for every gene and
species, integrate mirDIP predicted targeting with TarBase experimental
support, weighted by exosomal miRNA abundance. Also export the miRNA-gene edge
list used for the network figure.

WEIGHTED SCORE (per gene g, species s)
--------------------------------------
    score(g,s) = sum over miRNAs m in s of  a(m,s) * ( r(m,g) + b(m,g) )
where
    a(m,s) = mean relative exosomal abundance of m in s   (Avg_Percent, from step 01)
    r(m,g) = mirDIP INTEGRATED_RANK (V/H confidence classes only)
    b(m,g) = VALIDATION_BONUS (0.5) if the (m,g) pair is in TarBase, else 0

REQUIRED INPUT
--------------
1. Per-species TOP-90% matrices from 01_build_count_matrix.py:
       <species>/final_miRNA_CountMatrix_<species>_TOP90.csv
2. mirDIP unidirectional export (DOWNLOAD SEPARATELY, see README):
       mirDIP_Unidirectional_search_v_5_2/mirDIP_Unidirectional_search_v.5.2.txt
   (comma-separated; columns used: GENE_SYMBOL, MICRORNA, INTEGRATED_RANK, SCORE_CLASS)
3. TarBase v9 human interactions (DOWNLOAD SEPARATELY, see README):
       Homo_sapiens.tsv        (tab-separated; columns used: mirna_name, gene_name)

OUTPUT
------
    FINAL_Integrated_miRNA_Targets_Report.xlsx
        - sheet 'summary'  : per-gene score for each species + support statistics
                             (mirDIP_Unique_miRNAs, TarBase_Unique_miRNAs,
                              TarBase_Total_Experiments, Average_Targeting_Score)
        - one sheet per species with the per-species detail
    miRNA_gene_edges.csv.gz  : edge list (milk, miRNA, gene, mirdip_rank,
                               mirna_abundance_pct, TarBase_validated, edge_weight)

Note: the edge list is filtered to genes retained by 04_categorize_genes.py if
Categorized_Genes.xlsx already exists (keeps the file small); otherwise the full
edge list is written.

USAGE
-----
    python scripts/03_build_integrated_targetome.py
"""
import pandas as pd
import glob
import os

# --- CONFIG ---
VALIDATION_BONUS = 0.5
MIRDIP_FILE = "mirDIP_Unidirectional_search_v_5_2/mirDIP_Unidirectional_search_v.5.2.txt"
TARBASE_FILE = "Homo_sapiens.tsv"
REPORT_FILE = "FINAL_Integrated_miRNA_Targets_Report.xlsx"
EDGES_FILE = "miRNA_gene_edges.csv.gz"
# --------------

# 1. per-species miRNA abundance (TOP-90%)
print("Loading per-species miRNA spectra...")
top90 = glob.glob("*/final_miRNA_CountMatrix_*_TOP90.csv")
milk_data, milk_names = {}, []
unique_mirnas = set()
for f in top90:
    milk = os.path.basename(os.path.dirname(f))
    milk_names.append(milk)
    d = pd.read_csv(f)[["miRNA_ID", "Avg_Percent"]]
    milk_data[milk] = d
    unique_mirnas.update(d["miRNA_ID"].tolist())

# 2. databases
print("Loading mirDIP (V/H classes) and TarBase...")
mirdip = pd.concat(
    chunk[(chunk["SCORE_CLASS"].isin(["V", "H"])) & (chunk["MICRORNA"].isin(unique_mirnas))]
    for chunk in pd.read_csv(MIRDIP_FILE, sep=",", chunksize=10**6,
                             usecols=["GENE_SYMBOL", "MICRORNA", "INTEGRATED_RANK", "SCORE_CLASS"]))
tarbase = pd.read_csv(TARBASE_FILE, sep="\t", usecols=["mirna_name", "gene_name"])
tarbase = tarbase[tarbase["mirna_name"].isin(unique_mirnas)]
tarbase_pairs = set(tarbase["mirna_name"] + ":" + tarbase["gene_name"])
print(f"  TarBase experimental pairs in scope: {len(tarbase_pairs)}")

# 3. score per species (+ collect edges)
print("Scoring per species...")
all_genes = mirdip["GENE_SYMBOL"].unique()
summary = pd.DataFrame(index=all_genes); summary.index.name = "GENE_SYMBOL"
sheets, edge_frames = {}, []

for milk, dmir in milk_data.items():
    merged = pd.merge(mirdip, dmir, left_on="MICRORNA", right_on="miRNA_ID", how="inner")
    merged["Bonus"] = (merged["MICRORNA"] + ":" + merged["GENE_SYMBOL"]).isin(tarbase_pairs) * VALIDATION_BONUS
    merged["Weighted_Score"] = merged["Avg_Percent"] * (merged["INTEGRATED_RANK"] + merged["Bonus"])

    e = merged[["MICRORNA", "GENE_SYMBOL", "INTEGRATED_RANK", "Avg_Percent", "Bonus", "Weighted_Score"]].copy()
    e.insert(0, "milk", milk)
    e["TarBase_validated"] = (e["Bonus"] > 0).astype(int)
    edge_frames.append(e)

    per = pd.DataFrame(index=all_genes); per.index.name = "GENE_SYMBOL"
    per[f"Score_{milk}"] = merged.groupby("GENE_SYMBOL")["Weighted_Score"].sum()
    tb_milk = tarbase[tarbase["mirna_name"].isin(set(dmir["miRNA_ID"]))]
    per["mirDIP_Unique_miRNAs"] = merged.groupby("GENE_SYMBOL")["MICRORNA"].nunique()
    per["TarBase_Unique_miRNAs"] = tb_milk.groupby("gene_name")["mirna_name"].nunique()
    per["TarBase_Total_Experiments"] = tb_milk.groupby("gene_name")["mirna_name"].count()
    per = per.fillna(0)
    summary[milk] = per[f"Score_{milk}"]
    sheets[milk] = per.sort_values(f"Score_{milk}", ascending=False).round(3)

# 4. finalise summary
print("Finalising summary...")
summary = summary.fillna(0)
summary["mirDIP_Unique_miRNAs"] = mirdip.groupby("GENE_SYMBOL")["MICRORNA"].nunique()
summary["TarBase_Unique_miRNAs"] = tarbase.groupby("gene_name")["mirna_name"].nunique()
summary["TarBase_Total_Experiments"] = tarbase.groupby("gene_name")["mirna_name"].count()
for c in ["mirDIP_Unique_miRNAs", "TarBase_Unique_miRNAs", "TarBase_Total_Experiments"]:
    summary[c] = summary[c].fillna(0).astype(int)
summary["Average_Targeting_Score"] = summary[milk_names].mean(axis=1)
summary = summary.sort_values("Average_Targeting_Score", ascending=False).round(3)

with pd.ExcelWriter(REPORT_FILE, engine="openpyxl") as w:
    summary.to_excel(w, sheet_name="summary")
    for milk, d in sheets.items():
        d.to_excel(w, sheet_name=milk[:31])
print(f"Saved {REPORT_FILE} (sheets: summary, {', '.join(milk_names)})")

# 5. edge list (optionally restricted to categorised genes)
edges = pd.concat(edge_frames, ignore_index=True).rename(columns={
    "MICRORNA": "miRNA", "GENE_SYMBOL": "gene", "INTEGRATED_RANK": "mirdip_rank",
    "Avg_Percent": "mirna_abundance_pct", "Weighted_Score": "edge_weight"})
if os.path.exists("Categorized_Genes.xlsx"):
    keep = set(pd.read_excel("Categorized_Genes.xlsx")["GENE_SYMBOL"].astype(str))
    edges = edges[edges["gene"].isin(keep)]
    print(f"  edge list restricted to {len(keep)} categorised genes")
edges = edges[["milk", "miRNA", "gene", "mirdip_rank", "mirna_abundance_pct",
               "TarBase_validated", "edge_weight"]].round(4)
edges.to_csv(EDGES_FILE, index=False, compression="gzip")
print(f"Saved {EDGES_FILE} ({len(edges)} edges)")
