"""
08_enrichment_by_category.py
============================
Weighted functional enrichment (GO + KEGG) for each category, via Enrichr
(gseapy). Bars are ranked by a cumulative regulatory-impact score (sum of the
weighted gene scores of the pathway's genes); top driver genes are annotated.
Produces Supplementary S7.

REQUIRED INPUT (from 04_categorize_genes.py)
--------------------------------------------
    Categorized_Genes.xlsx
INTERNET ACCESS REQUIRED: gseapy queries the Enrichr web service.

OUTPUT
------
    Enrichment_Plots/Enrich_<Category>_<Library>.png

USAGE
-----
    python scripts/08_enrichment_by_category.py
"""
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
import os

INPUT_FILE = "Categorized_Genes.xlsx"
OUTPUT_DIR = "Enrichment_Plots/"
LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
COLS = ["human", "cow", "goat", "donkey"]
MIN_GENES = 10
os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_excel(INPUT_FILE)
df["GENE_SYMBOL"] = df["GENE_SYMBOL"].astype(str).str.upper()
df = df[df["Strict_Category"].notna() & (df["Strict_Category"] != "Low_Targeting_Noise")]

def wrap(t, width=45):
    return "\n".join(textwrap.wrap(t.split(" (GO")[0].replace("_", " ").strip().capitalize(), width))

for cat in df["Strict_Category"].unique():
    cat = str(cat)
    d = df[df["Strict_Category"] == cat].copy()
    genes = d["GENE_SYMBOL"].unique().tolist()
    if len(genes) < MIN_GENES:
        print(f"[skip] {cat}: {len(genes)} genes (< {MIN_GENES})"); continue
    if "Specific" in cat:
        dom = [s for s in cat.replace("_Specific", "").lower().split("_") if s in COLS]
        d["Weight"] = d[dom].mean(axis=1)
    else:
        d["Weight"] = d[COLS].mean(axis=1)
    wmap = dict(zip(d["GENE_SYMBOL"], d["Weight"]))
    for lib in LIBRARIES:
        try:
            res = gp.enrichr(gene_list=genes, gene_sets=[lib], organism="Human", outdir=None).results
            res = res[res["P-value"] < 0.05].head(15).copy()
            if res.empty:
                continue
            imp, gc, lab = [], [], []
            for gs in res["Genes"]:
                w = {g: wmap.get(g, 0) for g in gs.split(";") if g in wmap}
                imp.append(sum(w.values())); gc.append(len(w))
                lab.append(", ".join(f"{g} ({v:.1f})" for g, v in sorted(w.items(), key=lambda x: -x[1])[:3]))
            res["Impact_Score"], res["Gene_Count"], res["Top_3_Genes"] = imp, gc, lab
            res["Clean_Term"] = res["Term"].apply(wrap)
            res = res.sort_values("Impact_Score", ascending=False).head(15)
            palette = {"Human": "Greens", "Cow": "Blues", "Donkey": "Purples", "Goat": "Oranges"}
            pal = next((v for k, v in palette.items() if k in cat), "Reds")
            plt.figure(figsize=(12, 8)); sns.set_style("white")
            ax = sns.barplot(data=res, x="Impact_Score", y="Clean_Term", hue="Clean_Term",
                             palette=pal, legend=False, alpha=0.85)
            m = res["Impact_Score"].max(); ax.set_xlim(0, m * 1.35)
            for i, r in enumerate(res.itertuples()):
                ax.text(r.Impact_Score + m*0.02, i, f"p={r._4:.1e} (n={r.Gene_Count})",
                        va="center", ha="left", fontsize=10, fontweight="bold", color="#333")
                ax.text(m*0.01, i + 0.2, r.Top_3_Genes, fontsize=9.5, fontstyle="italic", va="top")
            plt.title(f"{cat.replace('_',' ')}\nWeighted {lib.split('_')[0]} Enrichment", fontsize=16, fontweight="bold", pad=20)
            plt.xlabel("Cumulative Regulatory Impact (Sum of Gene Scores)", fontsize=12); plt.ylabel("")
            sns.despine(left=True); plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"Enrich_{cat}_{lib.split('_')[0]}.png"), dpi=300, bbox_inches="tight")
            plt.close(); print(f"[ok] {cat} / {lib}")
        except Exception as e:
            print(f"[err] {cat} / {lib}: {e}")
print("Done.")
