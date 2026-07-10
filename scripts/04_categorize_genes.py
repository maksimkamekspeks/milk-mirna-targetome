"""
04_categorize_genes.py
======================
Assign every gene to a cross-species category from its per-species weighted
scores, and remove low-targeting noise.

RULE
----
1. Noise filter: T = PERCENTILE_THRESHOLD-th percentile of all non-zero scores.
   A gene is 'Low_Targeting_Noise' (removed) if max(score) < T OR sum(score) < MIN_SUM_SCORE.
2. Sort the four species by score; scan for the first 'gap' where the higher
   score is >= FOLD_CHANGE times the next. Species above the gap form the
   dominant group -> "<Sp1>_<Sp2>_..._Specific".
3. No gap: 'Core_Pan_Milk' if all species >= T, else 'Variable_Gradient'.

REQUIRED INPUT (from 03_build_integrated_targetome.py)
------------------------------------------------------
    FINAL_Integrated_miRNA_Targets_Report.xlsx  (sheet 'summary')

OUTPUT
------
    Categorized_Genes.xlsx   (retained genes + 'Strict_Category')

USAGE
-----
    python scripts/04_categorize_genes.py
Change FOLD_CHANGE / PERCENTILE_THRESHOLD / MIN_SUM_SCORE to reproduce the
threshold-sensitivity analysis (see 10_robustness_permutation.py).
"""
import pandas as pd
import numpy as np

# --- CONFIG ---
INPUT_FILE = "FINAL_Integrated_miRNA_Targets_Report.xlsx"
SHEET_NAME = "summary"
COLS = ["human", "cow", "goat", "donkey"]
PERCENTILE_THRESHOLD = 60     # noise cut-off percentile
FOLD_CHANGE = 1.5             # gap ratio defining dominance
MIN_SUM_SCORE = 10            # minimum summed score to keep a gene
OUTPUT_FILE = "Categorized_Genes.xlsx"
# --------------

df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, index_col="GENE_SYMBOL")
non_zero = df[COLS].values.flatten()
T = np.percentile(non_zero[non_zero > 0], PERCENTILE_THRESHOLD)
print(f"Noise threshold T (p{PERCENTILE_THRESHOLD}) = {T:.3f} | FC = {FOLD_CHANGE} | min sum = {MIN_SUM_SCORE}")

def categorize(row):
    scores = {c: row[c] for c in COLS}
    if max(scores.values()) < T or sum(scores.values()) < MIN_SUM_SCORE:
        return "Low_Targeting_Noise"
    ordered = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    dominant = []
    for i in range(len(ordered) - 1):
        nxt = ordered[i + 1][1] if ordered[i + 1][1] > 0 else 0.01
        if ordered[i][1] / nxt >= FOLD_CHANGE:
            dominant = [s[0] for s in ordered[:i + 1]]
            break
    if dominant:
        dominant.sort()
        return "_".join(s.capitalize() for s in dominant) + "_Specific"
    return "Core_Pan_Milk" if min(scores.values()) >= T else "Variable_Gradient"

df["Strict_Category"] = df.apply(categorize, axis=1)
counts = df["Strict_Category"].value_counts()
for cat, n in counts.items():
    print(f"  {cat}: {n}")

clean = df[df["Strict_Category"] != "Low_Targeting_Noise"].copy()
print(f"Removed {len(df) - len(clean)} noise genes (max < T or sum < {MIN_SUM_SCORE}); kept {len(clean)}")
clean.to_excel(OUTPUT_FILE)
print(f"Saved {OUTPUT_FILE}")
