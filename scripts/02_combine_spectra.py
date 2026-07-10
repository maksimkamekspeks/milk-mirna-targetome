"""
02_combine_spectra.py
=====================
Merge the per-species mean abundance (Avg_Percent) of the top miRNAs into a
single cross-species table (Supplementary S3: pooled most-abundant miRNAs).

REQUIRED INPUT (produced by 01_build_count_matrix.py)
-----------------------------------------------------
    cow/final_miRNA_CountMatrix_cow_full.csv
    human/final_miRNA_CountMatrix_human_full.csv
    donkey/final_miRNA_CountMatrix_donkey_full.csv
    goat/final_miRNA_CountMatrix_goat_full.csv

OUTPUT
------
    combined_miRNA.xlsx   (miRNA_ID, Avg_Percent [mean across species], File_Count, Cumulative_Percent)

USAGE
-----
    python scripts/02_combine_spectra.py
"""
import pandas as pd
from functools import reduce

FILES = [
    "cow/final_miRNA_CountMatrix_cow_full.csv",
    "human/final_miRNA_CountMatrix_human_full.csv",
    "donkey/final_miRNA_CountMatrix_donkey_full.csv",
    "goat/final_miRNA_CountMatrix_goat_full.csv",
]

frames = []
for i, f in enumerate(FILES):
    d = pd.read_csv(f)[["miRNA_ID", "Avg_Percent"]].copy()
    d.rename(columns={"Avg_Percent": f"Avg_Percent_{i+1}"}, inplace=True)
    frames.append(d)

merged = reduce(lambda l, r: pd.merge(l, r, on="miRNA_ID", how="outer"), frames)
avg_cols = [c for c in merged.columns if c.startswith("Avg_Percent_")]
merged["File_Count"] = merged[avg_cols].notna().sum(axis=1)      # in how many species detected
merged["Avg_Percent"] = merged[avg_cols].mean(axis=1)            # mean over species where present

result = (merged[["miRNA_ID", "Avg_Percent", "File_Count"]]
          .sort_values("Avg_Percent", ascending=False)
          .reset_index(drop=True))
result["Cumulative_Percent"] = result["Avg_Percent"].cumsum()
result.to_excel("combined_miRNA.xlsx", index=False)
print(f"Saved combined_miRNA.xlsx ({len(result)} miRNAs)")
