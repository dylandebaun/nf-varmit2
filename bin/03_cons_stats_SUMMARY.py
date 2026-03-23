#!/usr/bin/env python3
"""
03_cons_stats_SUMMARY.py
Combine per-individual missing-data TSVs into a summary TSV and dot plot.
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--cons_stats_fn_list", nargs="+", required=True,
                    help="Per-individual missing-data TSV files (from 02_cons_stats_INDIV.py)")
parser.add_argument("-o", "--out_dir", default=".",
                    help="Output directory (default: current directory)")
args = parser.parse_args()

out_dir = args.out_dir

# ── Aggregate ─────────────────────────────────────────────────────────────────

rows = []
for fn in args.cons_stats_fn_list:
    indiv_id = os.path.basename(fn).split("_")[0]
    df       = pd.read_csv(fn, sep="\t", index_col=0)
    rows.append([indiv_id,
                 df.loc["Total", "Missing_data"],
                 df.loc["Total", "Missing_data_PERCENTAGE"]])

summary = pd.DataFrame(rows, columns=["Indiv", "Missing_data", "Missing_data_PERCENTAGE"])
summary.to_csv(os.path.join(out_dir, "all_indivs_missing_data.tsv"), sep="\t", index=False)

# ── Dot plot ──────────────────────────────────────────────────────────────────

sns.set_theme(style="whitegrid")
g = sns.PairGrid(
    summary.sort_values("Missing_data_PERCENTAGE", ascending=False),
    x_vars=["Missing_data_PERCENTAGE"], y_vars=["Indiv"],
    height=10,
)
g.map(sns.stripplot, size=10, orient="h", jitter=False,
      palette="flare_r", linewidth=1, edgecolor="w")
g.set(xlim=(0, 100))

for ax in g.axes.flat:
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.set_xlabel("Missing data (%)", fontweight="bold", fontsize=18, labelpad=20)
    ax.set_ylabel("Individual",       fontweight="bold", fontsize=18, labelpad=20)

sns.despine(left=True, bottom=True)
plt.savefig(os.path.join(out_dir, "all_indivs_missing_data.pdf"),
            bbox_inches="tight", dpi=300)
