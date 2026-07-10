"""
12_mirna_gene_network.py
========================
Conserved-miRNA / shared-core-gene network (manuscript Figure 7): the miRNAs
detected in all four milks and their strongest predicted targets among the
Core_Pan_Milk genes, distinguishing TarBase-validated from predicted-only edges.

REQUIRED INPUT
--------------
    miRNA_gene_edges.csv.gz   (from 03_build_integrated_targetome.py)
    Categorized_Genes.xlsx    (from 04_categorize_genes.py; provides Core_Pan_Milk set)

OUTPUT
------
    Figure7_miRNA_gene_network.png

USAGE
-----
    python scripts/12_mirna_gene_network.py
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

EDGES = "miRNA_gene_edges.csv.gz"
CATS = "Categorized_Genes.xlsx"
TOP_MIRNAS = 12
TB, PRED = "#1b7f79", "#b9b9b9"

e = pd.read_csv(EDGES)
cat = pd.read_excel(CATS)
core = set(cat.loc[cat["Strict_Category"] == "Core_Pan_Milk", "GENE_SYMBOL"].astype(str))
cons = e.groupby("miRNA")["milk"].nunique()
cons = set(cons[cons == 4].index)                        # present in all 4 species
ab = (e[e.miRNA.isin(cons)].groupby(["miRNA","milk"])["mirna_abundance_pct"].first()
      .groupby("miRNA").mean())
topmir = ab.sort_values(ascending=False).head(TOP_MIRNAS).index.tolist()
sub = e[e.miRNA.isin(topmir) & e.gene.isin(core)]
agg = sub.groupby(["miRNA","gene"]).agg(w=("edge_weight","mean"), rank=("mirdip_rank","max"),
                                        tb=("TarBase_validated","max")).reset_index()
rows = []
for m, grp in agg.groupby("miRNA"):
    rows.append(pd.concat([grp[grp.tb==1].sort_values("w", ascending=False).head(2),
                           grp[grp.tb==0].sort_values("rank", ascending=False).head(1)]))
pick = pd.concat(rows)
genes = list(pick["gene"].unique())
deg = pick.groupby("gene")["miRNA"].nunique().to_dict()
spread = lambda n: np.linspace(0.96, 0.04, n)
my = dict(zip(topmir, spread(len(topmir))))
gy_raw = {g: np.mean([my[m] for m in pick[pick.gene==g]["miRNA"]]) for g in genes}
gs = sorted(genes, key=lambda g: -gy_raw[g]); gy = dict(zip(gs, spread(len(gs))))
xm, xg = 0.06, 0.94
fig, ax = plt.subplots(figsize=(9.4, 8.0)); ax.axis("off"); ax.set_xlim(-0.17, 1.17); ax.set_ylim(-0.04, 1.03)
wmax = pick["w"].max()
for _, r in pick.iterrows():
    ax.plot([xm, xg], [my[r["miRNA"]], gy[r["gene"]]], color=(TB if r["tb"] else PRED),
            lw=(0.9+1.8*(r["w"]/wmax) if r["tb"] else 0.8), alpha=(0.85 if r["tb"] else 0.6),
            ls=("-" if r["tb"] else (0, (3, 2))), zorder=1)
amax = ab[topmir].max()
for m in topmir:
    ax.scatter(xm, my[m], s=170+900*(ab[m]/amax), color="#34495e", edgecolor="white", zorder=3)
    ax.text(xm-0.035, my[m], m.replace("hsa-", ""), ha="right", va="center", fontsize=8.4, fontweight="bold")
dmax = max(deg.values())
for g in gs:
    ax.scatter(xg, gy[g], s=95+430*(deg[g]/dmax), color="#e0a458", edgecolor="white", zorder=3)
    ax.text(xg+0.035, gy[g], g, ha="left", va="center", fontsize=7.7, fontweight=("bold" if deg[g]>=3 else "normal"))
ax.text(xm, 1.01, "Conserved miRNAs\n(shared by all four milks)", ha="center", va="bottom", fontsize=9.6, fontweight="bold")
ax.text(xg, 1.01, "Shared core target genes", ha="center", va="bottom", fontsize=9.6, fontweight="bold")
leg = [Line2D([0],[0],color=TB,lw=2.2,label="TarBase v9 validated"),
       Line2D([0],[0],color=PRED,lw=1.4,ls=(0,(3,2)),label="mirDIP predicted only"),
       Line2D([0],[0],marker="o",color="w",markerfacecolor="#34495e",markersize=9,label="miRNA (size proportional to abundance)"),
       Line2D([0],[0],marker="o",color="w",markerfacecolor="#e0a458",markersize=9,label="gene (size proportional to # conserved miRNAs)")]
ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5,-0.07), ncol=2, fontsize=7.9, frameon=False)
plt.tight_layout(); plt.savefig("Figure7_miRNA_gene_network.png", dpi=300, bbox_inches="tight")
print("Saved Figure7_miRNA_gene_network.png")
