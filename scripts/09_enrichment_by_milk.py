"""
09_enrichment_by_milk.py
========================
Weighted functional enrichment (GO + KEGG) for each milk analysed globally
(all genes with that milk's score >= 5), via Enrichr (gseapy). Supplementary S8.

REQUIRED INPUT (from 04_categorize_genes.py)
--------------------------------------------
    Categorized_Genes.xlsx
INTERNET ACCESS REQUIRED (Enrichr web service).

OUTPUT
------
    Enrichment_By_Milk/Supp_Enrich_GLOBAL_<Milk>_<Library>.png

USAGE
-----
    python scripts/09_enrichment_by_milk.py
"""
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
import os

INPUT_FILE = "Categorized_Genes.xlsx"
COLS = ["human", "cow", "goat", "donkey"]
LIBRARIES = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
OUTPUT_DIR = "Enrichment_By_Milk/"
MIN_SCORE = 5
os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_excel(INPUT_FILE)
if "GENE_SYMBOL" not in df.columns:
    df = df.reset_index().rename(columns={df.columns[0]: "GENE_SYMBOL"})

def wrap(t, width=45):
    return "\n".join(textwrap.wrap(t.split(" (GO")[0].replace("_", " ").strip().capitalize(), width))

palette = {"human": "Greens", "cow": "Blues", "goat": "Oranges", "donkey": "Purples"}
for milk in COLS:
    d = df[df[milk] >= MIN_SCORE].sort_values(milk, ascending=False).copy()
    genes = d["GENE_SYMBOL"].astype(str).str.upper().tolist()
    wmap = dict(zip(genes, d[milk]))
    print(f"[{milk}] genes: {len(genes)}")
    for lib in LIBRARIES:
        try:
            res = gp.enrichr(gene_list=genes, gene_sets=[lib], organism="Human", outdir=None).results
            res = res[res["P-value"] < 0.05].head(20).copy()
            if res.empty:
                continue
            imp, gc, lab = [], [], []
            for gs in res["Genes"]:
                w = {g: wmap.get(g, 0) for g in gs.split(";") if g in wmap}
                imp.append(sum(w.values())); gc.append(len(w))
                lab.append(", ".join(f"{g} ({v:.1f})" for g, v in sorted(w.items(), key=lambda x: -x[1])[:3]))
            res["Impact_Score"], res["Gene_Count"], res["Top_3_Genes"] = imp, gc, lab
            res["Clean_Term"] = res["Term"].apply(wrap)
            res = res.sort_values("Impact_Score", ascending=False).head(20)
            plt.figure(figsize=(12, 8)); sns.set_theme(style="white")
            ax = sns.barplot(data=res, x="Impact_Score", y="Clean_Term", hue="Clean_Term",
                             palette=palette[milk], legend=False, alpha=0.85)
            m = res["Impact_Score"].max(); ax.set_xlim(0, m * 1.35)
            for i, r in enumerate(res.itertuples()):
                ax.text(r.Impact_Score + m*0.02, i, f"p={r._4:.1e} (n={r.Gene_Count})",
                        va="center", ha="left", fontsize=10, fontweight="bold", color="#333")
                ax.text(m*0.01, i + 0.2, r.Top_3_Genes, fontsize=9.5, fontstyle="italic", va="top")
            plt.title(f"{milk.capitalize()} Milk - Global (n={len(genes)})\nWeighted {lib.split('_')[0]} Enrichment",
                      fontsize=16, fontweight="bold", pad=20)
            plt.xlabel("Cumulative Regulatory Impact (Sum of Gene Scores)", fontsize=12); plt.ylabel("")
            sns.despine(left=True); plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"Supp_Enrich_GLOBAL_{milk.capitalize()}_{lib.split('_')[0]}.png"),
                        dpi=300, bbox_inches="tight")
            plt.close(); print(f"[ok] {milk} / {lib}")
        except Exception as e:
            print(f"[err] {milk} / {lib}: {e}")
print("Done.")
