"""
10_robustness_permutation.py
============================
Robustness of the category structure (manuscript Figure 5 / Supplementary S9):
(A) sensitivity of category sizes to the fold-change and noise thresholds;
(B) permutation null test of the species-specific category sizes
    (per-gene shuffling of the four species scores, 1000 permutations);
(C) Spearman correlation between the per-species weighted scores.

REQUIRED INPUT (from 03_build_integrated_targetome.py)
------------------------------------------------------
    FINAL_Integrated_miRNA_Targets_Report.xlsx  (sheet 'summary')

OUTPUT
------
    Figure5_robustness.png
    S9_robustness_analyses.xlsx   (sensitivity, permutation, correlation tables)

USAGE
-----
    python scripts/10_robustness_permutation.py
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import spearmanr
from collections import Counter

INPUT_FILE = "FINAL_Integrated_miRNA_Targets_Report.xlsx"
COLS = ["human", "cow", "goat", "donkey"]
H, C, G, D = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
rng = np.random.default_rng(42)

df = pd.read_excel(INPUT_FILE, sheet_name="summary", index_col="GENE_SYMBOL")
M = df[COLS].to_numpy(float)

def categorize(SM, FC=1.5, pct=60, minsum=10):
    nz = SM[SM > 0]; T = np.percentile(nz, pct); out = []
    for sc in SM:
        if sc.max() < T or sc.sum() < minsum:
            out.append("Noise"); continue
        o = np.argsort(sc)[::-1]; sv = sc[o]; dom = None
        for i in range(3):
            nx = sv[i+1] if sv[i+1] > 0 else 0.01
            if sv[i]/nx >= FC:
                dom = sorted([COLS[o[j]] for j in range(i+1)]); break
        out.append(("_".join(c.capitalize() for c in dom) + "_Specific") if dom
                   else ("Core_Pan_Milk" if sc.min() >= T else "Variable_Gradient"))
    return np.array(out)

FCs = [1.3, 1.5, 2.0, 2.5, 3.0]
series = {k: [] for k in ["Core_Pan_Milk","Donkey_Specific","Human_Specific","Cow_Goat_Specific","Variable_Gradient"]}
for fc in FCs:
    cc = Counter(categorize(M, FC=fc))
    for k in series: series[k].append(cc.get(k, 0))

base = categorize(M); bc = Counter(base)
targets = ["Donkey_Specific", "Human_Specific", "Cow_Goat_Specific"]
obs = {t: bc.get(t, 0) for t in targets}
null = {t: [] for t in targets}
for _ in range(1000):
    cc = Counter(categorize(np.array([r[rng.permutation(4)] for r in M])))
    for t in targets: null[t].append(cc.get(t, 0))

Mr = M[base != "Noise"]
S = np.ones((4, 4))
for i in range(4):
    for j in range(4):
        if i != j: S[i, j] = spearmanr(Mr[:, i], Mr[:, j]).correlation

# ---- figure ----
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0))
cm = {"Core_Pan_Milk":"#444","Donkey_Specific":D,"Human_Specific":H,"Cow_Goat_Specific":C,"Variable_Gradient":"#999"}
lab = {"Core_Pan_Milk":"Core pan-milk","Donkey_Specific":"Donkey-specific","Human_Specific":"Human-specific","Cow_Goat_Specific":"Cow-goat","Variable_Gradient":"Variable gradient"}
for k, v in series.items(): ax[0].plot(FCs, v, marker="o", lw=1.8, color=cm[k], label=lab[k])
ax[0].axvline(1.5, ls="--", color="k", lw=0.8, alpha=0.6)
ax[0].set_xlabel("Fold-change threshold"); ax[0].set_ylabel("Number of genes")
ax[0].set_title("A  Sensitivity to FC threshold", fontweight="bold", loc="left"); ax[0].legend(fontsize=7.5, frameon=False)
ax[0].spines[["top","right"]].set_visible(False)
tcol = {"Donkey_Specific":D,"Human_Specific":H,"Cow_Goat_Specific":C}
for t in targets:
    ax[1].hist(null[t], bins=30, color=tcol[t], alpha=0.45); ax[1].axvline(obs[t], color=tcol[t], lw=2.2)
ax[1].set_xlabel("Genes per category (null vs observed)"); ax[1].set_ylabel("Permutations")
ax[1].set_title("B  Permutation null (1000 shuffles)", fontweight="bold", loc="left")
ax[1].spines[["top","right"]].set_visible(False)
Sm = S.copy()
for i in range(4):
    for j in range(4):
        if j < i: Sm[i, j] = np.nan
cmap = mpl.cm.YlGnBu.copy(); cmap.set_bad("white")
im = ax[2].imshow(Sm, vmin=0.6, vmax=1.0, cmap=cmap)
ax[2].set_xticks(range(4)); ax[2].set_yticks(range(4))
ax[2].set_xticklabels([c.capitalize() for c in COLS]); ax[2].set_yticklabels([c.capitalize() for c in COLS])
for i in range(4):
    for j in range(4):
        if j >= i: ax[2].text(j, i, f"{S[i,j]:.2f}", ha="center", va="center", fontsize=8.5,
                              color="white" if S[i,j] > 0.85 else "black")
for s in ["top","right","bottom","left"]: ax[2].spines[s].set_visible(False)
ax[2].set_title("C  Cross-species score correlation", fontweight="bold", loc="left")
fig.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04).set_label("Spearman rho", fontsize=8)
plt.tight_layout(); plt.savefig("Figure5_robustness.png", dpi=300, bbox_inches="tight")

# ---- S9 tables ----
sens = pd.DataFrame({"Fold_change": FCs, **{k: series[k] for k in series}})
perm = pd.DataFrame([{"Category": t, "Observed": obs[t], "Null_mean": round(np.mean(null[t]),1),
                      "Null_SD": round(np.std(null[t]),1),
                      "Empirical_p": (np.sum(np.array(null[t]) >= obs[t]) + 1)/(1000+1)} for t in targets])
corr = pd.DataFrame(S, index=[c.capitalize() for c in COLS], columns=[c.capitalize() for c in COLS])
with pd.ExcelWriter("S9_robustness_analyses.xlsx") as w:
    sens.to_excel(w, sheet_name="FC_sensitivity", index=False)
    perm.to_excel(w, sheet_name="Permutation_null", index=False)
    corr.to_excel(w, sheet_name="Spearman_correlation")
print("Saved Figure5_robustness.png and S9_robustness_analyses.xlsx")
