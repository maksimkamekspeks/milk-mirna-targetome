"""
06_plot_regulatory_power.py
===========================
Figure: distribution of total regulatory weight (summed score) across the
categories, per species. Reproduces the category-mass overview figure.

REQUIRED INPUT (from 04_categorize_genes.py)
--------------------------------------------
    Categorized_Genes.xlsx

OUTPUT
------
    Regulatory_Power_Distribution.png

USAGE
-----
    python scripts/06_plot_regulatory_power.py
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_FILE = "Categorized_Genes.xlsx"
COLS = ["human", "cow", "goat", "donkey"]
OUTPUT_IMAGE = "Regulatory_Power_Distribution.png"

try:
    df = pd.read_excel(INPUT_FILE, index_col=0)
except Exception:
    df = pd.read_excel(INPUT_FILE)
if "Strict_Category" not in df.columns:
    raise SystemExit("Column 'Strict_Category' not found in input.")

agg = df.groupby("Strict_Category")[COLS].sum()
counts = df["Strict_Category"].value_counts()
agg["Total_Mass"] = agg[COLS].sum(axis=1)
global_mass = agg["Total_Mass"].sum()
agg = agg.sort_values("Total_Mass", ascending=True)

cats = agg.index.tolist()
y = np.arange(len(cats)); height = 0.2
fig, ax = plt.subplots(figsize=(14, 12))
sns.set_theme(style="whitegrid", rc={"grid.alpha": 0.5, "grid.linestyle": "--"})
colors = ["#4A90E2", "#F5A623", "#50E3C2", "#9013FE"]
for i, col in enumerate(COLS):
    offset = (i - len(COLS)/2 + 0.5) * height
    ax.barh(y + offset, agg[col], height=height, label=col.capitalize(),
            color=colors[i], edgecolor="black", linewidth=0.7, alpha=0.9)
ax.set_yticks(y)
ax.set_yticklabels([f"{c.replace('_',' ')}\n(n={counts[c]})" for c in cats], fontsize=11, fontweight="bold")
for i, c in enumerate(cats):
    pct = agg.loc[c, "Total_Mass"] / global_mass * 100
    x = agg.loc[c, COLS].max() + agg[COLS].max().max() * 0.02
    ax.text(x, y[i], f"{pct:.1f}%", va="center", ha="left", fontweight="bold", color="#B00020", fontsize=12)
ax.set_xlabel("Regulatory Weight (Score Sum)", fontweight="bold", fontsize=13)
ax.set_title("Regulatory Power Distribution Across Categories", fontweight="bold", fontsize=18, pad=20)
ax.legend(title="Milk Type", fontsize=12, title_fontsize=13, loc="lower right", frameon=True, shadow=True)
sns.despine(top=True, right=True)
plt.tight_layout(); plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches="tight"); plt.close()
print(f"Saved {OUTPUT_IMAGE}")
