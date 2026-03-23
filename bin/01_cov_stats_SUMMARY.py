#!/usr/bin/env python3
"""
01_cov_stats_SUMMARY.py
Combine per-individual coverage TSVs and produce a summary TSV and violin plot
showing mean sequencing depth per individual, split by chromosome type.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--indiv_cov_tsv_list", nargs="+", required=True,
                    help="Per-individual coverage TSV files (from 00_cov_stats_INDIV.py)")
parser.add_argument("-s", "--sex_chr_list", nargs="+", default=[],
                    help="Sex chromosome names")
parser.add_argument("-p", "--plots_per_panel", type=int, default=10,
                    help="Individuals per panel row (default: 10)")
parser.add_argument("-o", "--out_dir", default=".",
                    help="Output directory (default: current directory)")
args = parser.parse_args()

out_dir = args.out_dir

# Handle Nextflow list-passing artifacts
sex_chromos = {c.strip(" [],'\"") for c in args.sex_chr_list if c.strip(" [],'\"}")}

# ── Load and annotate ─────────────────────────────────────────────────────────

frames = []
for tsv in args.indiv_cov_tsv_list:
    indiv_id = os.path.basename(tsv).split("_coverage")[0]
    df = pd.read_csv(tsv, sep="\t", header=0)
    df["individual"] = indiv_id
    frames.append(df)

df = pd.concat(frames, ignore_index=True)
df["chromo_type"] = df["#rname"].apply(
    lambda c: "sex chromosome" if c in sex_chromos else "autosome"
)

indivs = df["individual"].unique()

# ── Summary TSV ───────────────────────────────────────────────────────────────

summary_frames = []
for indiv in indivs:
    sub = df[df["individual"] == indiv]
    grp = sub.groupby("chromo_type", as_index=False).mean(numeric_only=True)
    grp["individual"] = indiv
    summary_frames.append(grp[["individual", "chromo_type", "coverage", "meandepth"]])

df_summary = pd.concat(summary_frames)
df_summary.to_csv(os.path.join(out_dir, "all_indivs_summary.tsv"), sep="\t", index=False)

# ── Violin plot ───────────────────────────────────────────────────────────────

n_per_panel   = args.plots_per_panel
upper_depth   = int(df["meandepth"].max())
n_panels      = int(np.ceil(len(indivs) / n_per_panel))
indiv_chunks  = [indivs[i:i + n_per_panel] for i in range(0, len(indivs), n_per_panel)]

sns.set_style("whitegrid", {"axes.grid": False})
fig, axes = plt.subplots(nrows=max(n_panels, 1), sharex=False, figsize=(19.6, 19.6))
fig.subplots_adjust(hspace=0.2)

if n_panels == 1:
    axes = [axes]

PALETTE = {"autosome": "#9b59b6", "sex chromosome": "#3498db"}

for i, chunk in enumerate(indiv_chunks):
    ax  = axes[i]
    sub = df[df["individual"].isin(chunk)]
    sns.violinplot(
        data=sub, x="individual", y="meandepth", hue="chromo_type",
        split=True, inner="quart", fill=False, palette=PALETTE,
        ax=ax, legend=False,
    )
    ax.set(ylim=(0, upper_depth))
    ax.set_xticklabels(ax.get_xticklabels(), fontweight="bold", fontsize=12, rotation=45)
    ax.set_yticklabels(ax.get_yticklabels(), fontweight="bold", fontsize=12)
    ax.set_ylabel("Mean depth", fontweight="bold", fontsize=18, labelpad=20)
    ax.set_xlabel("")

axes[-1].set_xlabel("Individuals", fontweight="bold", fontsize=18, labelpad=20)

legend_patches = [mpatches.Patch(color=v, label=k) for k, v in PALETTE.items()]
axes[0].legend(handles=legend_patches, bbox_to_anchor=(1.05, 0.8), fontsize=14, frameon=False)

sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "all_indivs_summary.pdf"),
            format="pdf", bbox_inches="tight", dpi=300)
