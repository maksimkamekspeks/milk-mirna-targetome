#!/usr/bin/env bash
# ============================================================================
# 00_prepare_reads.sh  —  Step 0: from SRA accessions to per-sample miRBase
#                          alignments (SAM), consumed by 01_build_count_matrix.py
#
# Example small-RNA workflow:
#   SRA download (sra-tools) -> adapter trimming + length filter (cutadapt)
#   -> alignment to human miRBase MATURE miRNAs (bowtie) -> SAM.
# In the resulting SAM, column 3 (reference name) is the miRNA ID, which is what
# 01_build_count_matrix.py counts.
#
# Requirements (see environment.yml): sra-tools, cutadapt, bowtie
#
# IMPORTANT: the 3' adapter depends on the library kit and DIFFERS between
# datasets. Set ADAPTER for each dataset (common examples given below). Run this
# script once per species/dataset (edit SPECIES and SRR_LIST).
# ============================================================================
set -euo pipefail

# ----------------------- CONFIG (edit me) -----------------------
SPECIES="human"                                   # output folder: human/cow/goat/donkey
SRR_LIST=(SRR346516 SRR346517 SRR346518 SRR346519)# accessions for THIS dataset
ADAPTER="TGGAATTCTCGGGTGCCAAGG"                   # 3' adapter — CHANGE per dataset:
#   Illumina TruSeq small RNA : TGGAATTCTCGGGTGCCAAGG
#   NEBNext small RNA         : AGATCGGAAGAGCACACGTCT
MIN_LEN=18
MAX_LEN=30
THREADS=4
MIRBASE_INDEX="ref/hsa_mature"                    # bowtie index prefix (built once, below)
# ----------------------------------------------------------------

mkdir -p "$SPECIES" ref

# ---- Build the miRBase index ONCE (uncomment and run the first time) ----
# 1) Download mature.fa from https://www.mirbase.org/download/
# 2) Keep human entries and convert RNA (U) to DNA (T), then build the index:
#   awk '/^>hsa/{p=1;print;next} /^>/{p=0} p' mature.fa > ref/hsa_mature_rna.fa
#   sed '/^[^>]/ y/uU/tT/' ref/hsa_mature_rna.fa > ref/hsa_mature.fa
#   bowtie-build ref/hsa_mature.fa ref/hsa_mature

for SRR in "${SRR_LIST[@]}"; do
  echo ">> $SRR"

  # 1. download FASTQ (skip if ${SRR}.fastq already present)
  prefetch "$SRR"
  fasterq-dump "$SRR" -O . --split-files

  # small-RNA libraries are single-end; use ${SRR}.fastq (or ${SRR}_1.fastq)
  RAW="${SRR}.fastq"; [ -f "$RAW" ] || RAW="${SRR}_1.fastq"

  # 2. adapter trim + keep 18-30 nt reads
  cutadapt -a "$ADAPTER" -m "$MIN_LEN" -M "$MAX_LEN" -j "$THREADS" \
           -o "clean_${SRR}.fastq" "$RAW"

  # 3. align to miRBase mature (sense strand; miRNA name -> SAM column 3)
  bowtie -p "$THREADS" -S -v 1 -k 1 --norc \
         "$MIRBASE_INDEX" "clean_${SRR}.fastq" "${SPECIES}/clean_${SRR}.sam"
done

echo "Done -> ${SPECIES}/clean_*.sam  (now run scripts/01_build_count_matrix.py)"
