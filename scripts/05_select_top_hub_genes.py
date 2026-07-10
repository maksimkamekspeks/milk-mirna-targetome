"""
05_select_top_hub_genes.py
==========================
Select the top genes per category (Supplementary S5) for focused discussion.
For '*_Specific' categories genes are ranked by the mean score of the dominant
species and by the actual fold-change over the non-dominant species; for
Core_Pan_Milk / Variable_Gradient they are ranked by total score.

REQUIRED INPUT (from 04_categorize_genes.py)
--------------------------------------------
    Categorized_Genes.xlsx

OUTPUT
------
    Top_Hub_Genes_For_Discussion.xlsx

USAGE
-----
    python scripts/05_select_top_hub_genes.py
"""
import pandas as pd

INPUT_FILE = "Categorized_Genes.xlsx"
COLS = ["human", "cow", "goat", "donkey"]
MAX_TOP_GENES = 50
OUTPUT_FILE = "Top_Hub_Genes_For_Discussion.xlsx"

df = pd.read_excel(INPUT_FILE, index_col="GENE_SYMBOL").dropna(subset=["Strict_Category"])
out = []
for cat in df["Strict_Category"].unique():
    cat = str(cat)
    d = df[df["Strict_Category"] == cat].copy()
    if "Specific" in cat:
        dom = [s.lower() for s in cat.replace("_Specific", "").split("_")]
        non = [c for c in COLS if c not in dom]
        d["Dominant_Score"] = d[dom].mean(axis=1)
        if non:
            d["Non_Dominant_Max"] = d[non].max(axis=1)
            d["Actual_FC"] = d["Dominant_Score"] / (d["Non_Dominant_Max"] + 0.01)
            d = d.sort_values(["Dominant_Score", "Actual_FC"], ascending=[False, False])
        else:
            d = d.sort_values("Dominant_Score", ascending=False); d["Actual_FC"] = None
    else:
        d["Total_Score"] = d[COLS].sum(axis=1)
        d = d.sort_values("Total_Score", ascending=False)
        d["Dominant_Score"] = d["Total_Score"] / 4; d["Actual_FC"] = None
    top = d.head(min(MAX_TOP_GENES, len(d))).copy()
    top["Hub_Category"] = cat
    out.append(top)
    print(f"  {cat}: {len(top)}")

final = pd.concat(out)
keep = [c for c in COLS + ["Strict_Category", "Dominant_Score", "Actual_FC",
        "mirDIP_Unique_miRNAs", "TarBase_Unique_miRNAs", "TarBase_Total_Experiments"] if c in final.columns]
final[keep].to_excel(OUTPUT_FILE)
print(f"Saved {OUTPUT_FILE} ({len(final)} genes)")
