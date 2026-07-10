"""
11_sample_clustering_pca.py
===========================
Sample-level comparability (manuscript Figure 6): PCA and Ward hierarchical
clustering of the per-sample miRNA composition, showing that samples group by
species despite heterogeneous dataset origins. Also prints the mean silhouette
width and the between/within-species distance ratio.

REQUIRED INPUT (from 01_build_count_matrix.py)
----------------------------------------------
    <species>/final_miRNA_CountMatrix_<species>_full.csv   for each species
    (these contain one column per sample plus Avg_Percent, Cumulative_Percent)

OUTPUT
------
    Figure6_sample_clustering.png

USAGE
-----
    python scripts/11_sample_clustering_pca.py
"""
import pandas as pd, numpy as np, glob, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

SPECIES = ["human", "cow", "goat", "donkey"]
COL = {"human":"#0072B2","cow":"#D55E00","goat":"#009E73","donkey":"#CC79A7"}
NON_SAMPLE = {"miRNA_ID", "Avg_Percent", "Cumulative_Percent"}

samples, sp_of = {}, {}
allmir = set()
for sp in SPECIES:
    fs = glob.glob(os.path.join(sp, f"final_miRNA_CountMatrix_{sp}_full.csv"))
    if not fs:
        print(f"[warn] no full matrix for {sp}"); continue
    d = pd.read_csv(fs[0])
    mircol = d.columns[0]
    scols = [c for c in d.columns if c not in NON_SAMPLE and c != mircol]
    for c in scols:
        samples[c] = dict(zip(d[mircol], d[c])); sp_of[c] = sp
    allmir.update(d[mircol].tolist())

allmir = sorted(allmir); sn = list(samples.keys())
M = np.array([[samples[s].get(m, 0) for m in allmir] for s in sn], float)
M = M / M.sum(axis=1, keepdims=True)          # relative abundance per sample
X = np.log10(M + 1e-6); X = X - X.mean(0)
lab = np.array([sp_of[s] for s in sn])

# silhouette + distance ratio
Dm = squareform(pdist(X)); sil = []
for i in range(len(sn)):
    same = [Dm[i,j] for j in range(len(sn)) if lab[j]==lab[i] and j!=i]
    a = np.mean(same) if same else 0
    b = min(np.mean([Dm[i,j] for j in range(len(sn)) if lab[j]==o]) for o in set(lab) if o!=lab[i])
    sil.append((b-a)/max(a,b))
wm = np.mean([Dm[i,j] for i in range(len(sn)) for j in range(i+1,len(sn)) if lab[i]==lab[j]])
bm = np.mean([Dm[i,j] for i in range(len(sn)) for j in range(i+1,len(sn)) if lab[i]!=lab[j]])
print(f"samples={len(sn)} miRNAs={len(allmir)} | silhouette={np.mean(sil):.3f} | between/within={bm/wm:.2f}")

U, S, Vt = np.linalg.svd(X, full_matrices=False)
PC = U[:, :2]*S[:2]; ev = (S**2)/np.sum(S**2)*100
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for sp in SPECIES:
    idx = [i for i,x in enumerate(lab) if x==sp]
    ax[0].scatter(PC[idx,0], PC[idx,1], s=70, c=COL[sp], label=f"{sp.capitalize()} (n={len(idx)})", edgecolor="white")
ax[0].set_xlabel(f"PC1 ({ev[0]:.0f}%)"); ax[0].set_ylabel(f"PC2 ({ev[1]:.0f}%)")
ax[0].set_title("A  PCA of samples (miRNA composition)", fontweight="bold", loc="left")
ax[0].legend(fontsize=8, frameon=False); ax[0].spines[["top","right"]].set_visible(False)
Z = linkage(X, method="ward")
dendrogram(Z, labels=[sp_of[s][:2] for s in sn], ax=ax[1], color_threshold=0, above_threshold_color="#888")
for t in ax[1].get_xmajorticklabels():
    t.set_color({"hu":COL["human"],"co":COL["cow"],"go":COL["goat"],"do":COL["donkey"]}.get(t.get_text(), "k"))
ax[1].set_title("B  Ward hierarchical clustering", fontweight="bold", loc="left"); ax[1].set_ylabel("Distance")
ax[1].spines[["top","right"]].set_visible(False)
plt.tight_layout(); plt.savefig("Figure6_sample_clustering.png", dpi=300, bbox_inches="tight")
print("Saved Figure6_sample_clustering.png")
