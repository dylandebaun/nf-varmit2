#!/usr/bin/env python3
"""
02_cons_stats_INDIV.py
Calculate missing data (N/n/- characters) per chromosome in a per-individual
consensus FASTA, and write a summary TSV.

Expected input filenames: <indiv>_<chromo>_cons.fa  (one record per file)
"""

import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--cons_fn_list", nargs="+", required=True,
                    help="Per-chromosome consensus FASTA files")
parser.add_argument("-i", "--individual", required=True,
                    help="Individual identifier")
parser.add_argument("-o", "--out_dir", default=".",
                    help="Output directory (default: current directory)")
args = parser.parse_args()

rows = []
for fa in args.cons_fn_list:
    chromo  = os.path.basename(fa).split("_")[1]
    record  = SeqIO.read(fa, "fasta")
    length  = len(record.seq)
    missing = record.seq.count("N") + record.seq.count("n") + record.seq.count("-")
    rows.append([chromo, length, missing, round(missing / length * 100, 2)])

df = pd.DataFrame(rows, columns=["Chromo", "Length", "Missing_data", "Missing_data_PERCENTAGE"])

total_len     = df["Length"].sum()
total_missing = df["Missing_data"].sum()
df = pd.concat([df, pd.DataFrame([[
    "Total", total_len, total_missing,
    round(total_missing / total_len * 100, 2)
]], columns=df.columns)], ignore_index=True)

out_fn = os.path.join(args.out_dir, f"{args.individual}_missing_data.tsv")
df.to_csv(out_fn, sep="\t", index=False)
