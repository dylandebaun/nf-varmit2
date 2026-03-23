#!/usr/bin/env python3
"""
00_cov_stats_INDIV.py
Concatenate per-chromosome samtools coverage TSVs into a single per-individual file.
"""

import os
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--chromo_cov_tsv_list", nargs="+", required=True,
                    help="Per-chromosome samtools coverage TSV files")
parser.add_argument("-i", "--individual", required=True,
                    help="Individual identifier")
parser.add_argument("-o", "--out_dir", default=".",
                    help="Output directory (default: current directory)")
args = parser.parse_args()

frames = [pd.read_csv(f, sep="\t", header=0) for f in args.chromo_cov_tsv_list]
indiv_df = pd.concat(frames)

out_fn = os.path.join(args.out_dir, f"{args.individual}_coverage.tsv")
indiv_df.to_csv(out_fn, sep="\t", index=False)
