"""
07_plot_heatmap.py
==================
Figure: heatmap of the top genes per category, ranked by category weight.
Colour is row-normalised (row max = 1.0), i.e. it shows the cross-species
pattern of each gene, not its absolute score (per-category peak scores are
annotated on the right).

REQUIRED INPUT (from 04_categorize_genes.py)
--------------------------------------------
    Categorized_Genes.xlsx

OUTPUT
------
    Heatmap_Strict_Categories_Ranked_Scores.png

USAGE
-----
    python scripts/07_plot_heatmap.py
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_FILE = "Categorized_Genes.xlsx"
OUTPUT_FILE = "Heatmap_Strict_Categories_Ranked_Scores.png"
COLS = ["human", "cow", "goat", "donkey"]
TOP_N = 7

df = pd.read_excel(INPUT_FILE)
valid = df[df["Strict_Category"].notna() & (df["Strict_Category"] != "Low_Targeting_Noise")].copy()
valid["Gene_Total_Score"] = valid[COLS].sum(axis=1)
order = valid.groupby("Strict_Category")["Gene_Total_Score"].sum().sort_values(ascending=False).index.tolist()

blocks, sizes, peaks = [], [], []
for cat in order:
    d = valid[valid["Strict_Category"] == cat].copy()
    if d.empty:
        continue
    d["Max_Score"] = d[COLS].max(axis=1)
    top = d.sort_values("Max_Score", ascending=False).head(TOP_N)
    blocks.append(top); sizes.append(len(top)); peaks.append(top["Max_Score"].max())

final = pd.concat(blocks).set_index("GENE_SYMBOL")
data = final[COLS].astype(float)
scaled = data.div(data.max(axis=1), axis=0)

fig, (ax, cax) = plt.subplots(1, 2, figsize=(14, max(10, len(final)*0.28)),
                              gridspec_kw={"width_ratios": [1, 0.05], "wspace": 0.8})
sns.heatmap(scaled, cmap="viridis", ax=ax, cbar_ax=cax,
            cbar_kws={"label": "Relative Targeting Strength (Row Max = 1.0)"},
            linewidths=0.5, linecolor="lightgray")
y = 0
for cat, size, peak in zip(order, sizes, peaks):
    y += size
    if y < len(scaled):
        ax.axhline(y, color="black", lw=2.5)
    ax.text(len(COLS) + 0.2, y - size/2, f"{cat.replace('_',' ')}\n(Max Score: {peak:.1f})",
            va="center", ha="left", fontsize=11, fontweight="bold", color="#222")
ax.set_title("Signature miRNA Targets Across Milk Types\n(Categories Ranked by Total Weight)",
             fontsize=16, fontweight="bold", pad=20)
ax.set_ylabel("Target Genes", fontsize=14, fontweight="bold", labelpad=15)
ax.set_xlabel("Milk Type", fontsize=14, fontweight="bold", labelpad=10)
ax.set_xticklabels([c.capitalize() for c in COLS], rotation=45, ha="right", fontsize=13, fontweight="bold")
ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, fontstyle="italic")
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight"); plt.close()
print(f"Saved {OUTPUT_FILE}")
