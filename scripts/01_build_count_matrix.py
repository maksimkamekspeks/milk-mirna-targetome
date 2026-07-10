"""
01_build_count_matrix.py
========================
Quantify miRNA expression from aligned reads (SAM) into a per-species count
matrix, normalise each sample to relative abundance (%), and export the full
matrix and the TOP-90% subset used downstream.

WHAT IT DOES
------------
For each species it reads all per-sample SAM alignments, counts reads per
miRNA (column 3 = reference/miRNA name, isoforms summed), normalises every
sample to 100 %, averages the per-sample percentages, sorts by mean abundance,
computes the cumulative percentage, and writes:
  * <species>/final_miRNA_CountMatrix_<species>_full.csv   (all miRNAs; per-sample counts + Avg_Percent + Cumulative_Percent)
  * <species>/final_miRNA_CountMatrix_<species>_TOP90.csv  (miRNAs that jointly make up 90 % of the reads)

REQUIRED INPUT (you must provide these — NOT included in the code archive)
--------------------------------------------------------------------------
Per species, a folder named exactly as the species, containing the aligned
small-RNA reads as SAM files, one per sample:
    human/clean_SRR346516.sam, human/clean_SRR346517.sam, ...
    cow/clean_*.sam
    goat/clean_*.sam
    donkey/clean_*.sam
The SAM files are produced upstream (Step 0, see README): adapter trimming and
alignment of each sample's small-RNA reads to the human miRBase reference.

OUTPUT
------
The two CSV files per species listed above.

USAGE
-----
    python scripts/01_build_count_matrix.py
Edit SPECIES below to match your folder names.
"""
import pandas as pd
import glob
import os

# --- CONFIG ---
SPECIES = ["human", "cow", "goat", "donkey"]   # one folder per species, holding *.sam
TOP_CUMULATIVE = 90.0                          # keep miRNAs up to this cumulative %
# --------------

for sp in SPECIES:
    sam_files = glob.glob(os.path.join(sp, "*.sam"))
    if not sam_files:
        print(f"[skip] no .sam files found in ./{sp}/")
        continue
    print(f"[{sp}] reading {len(sam_files)} SAM file(s)...")

    all_counts = {}
    for f in sam_files:
        sample = os.path.basename(f).replace("clean_", "").replace(".sam", "")
        # SAM column 3 (index 2) = reference (miRNA) name; comment='@' skips headers
        d = pd.read_csv(f, sep="\t", comment="@", header=None, usecols=[2], names=["miRNA"])
        d = d[d["miRNA"] != "*"]               # drop unmapped reads
        all_counts[sample] = d["miRNA"].value_counts()

    matrix = pd.DataFrame(all_counts).fillna(0).astype(int)
    matrix.index.name = "miRNA_ID"

    # normalise each sample (column) to 100 %, then average across samples
    pct = matrix.div(matrix.sum(axis=0), axis=1) * 100
    matrix["Avg_Percent"] = pct.mean(axis=1)
    matrix = matrix.sort_values("Avg_Percent", ascending=False)
    matrix["Cumulative_Percent"] = matrix["Avg_Percent"].cumsum()

    matrix.to_csv(os.path.join(sp, f"final_miRNA_CountMatrix_{sp}_full.csv"))

    # keep the miRNAs that jointly reach TOP_CUMULATIVE % (including the crossing one)
    top = matrix[matrix["Cumulative_Percent"].shift(fill_value=0) < TOP_CUMULATIVE]
    top.to_csv(os.path.join(sp, f"final_miRNA_CountMatrix_{sp}_TOP90.csv"))

    print(f"[{sp}] miRNAs total: {len(matrix)} | kept (~{TOP_CUMULATIVE:.0f}%): {len(top)}")
print("Done.")
