#!/usr/bin/env python3
"""
00_calc_depth_thresholds.py
Calculate per-individual sequencing depth thresholds from samtools coverage output.

Thresholds are multiplier-based relative to mean depth across qualifying scaffolds,
computed separately for autosomes and sex chromosomes:
    min_cov = max(5, floor(mean * min_depth_multiplier))
    max_cov = ceil(mean * max_depth_multiplier)

This approach is more robust than mean ± 2SD because scaffold-level variance inflates
SD estimates in non-model genomes (especially near sex chromosomes and repeats).
"""

import os
import sys
import math
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--chromo_cov_tsv_list", nargs="+", required=True,
                    help="Per-chromosome samtools coverage TSV files")
parser.add_argument("-i", "--individual", required=True,
                    help="Individual identifier")
parser.add_argument("-s", "--sex_chromos", nargs="+", default=[],
                    help="Sex chromosome names (handles Nextflow bracket/comma formatting)")
parser.add_argument("-l", "--auto_chromos_min_length", type=int, default=0,
                    help="Minimum scaffold length (bp) for autosome depth calculation (default: 0)")
parser.add_argument("--min_depth_multiplier", type=float, default=0.5,
                    help="Lower threshold = mean × this value (default: 0.5)")
parser.add_argument("--max_depth_multiplier", type=float, default=2.0,
                    help="Upper threshold = mean × this value (default: 2.0)")
args = parser.parse_args()

indiv                   = args.individual
auto_chromos_min_length = args.auto_chromos_min_length
min_mult                = args.min_depth_multiplier
max_mult                = args.max_depth_multiplier

# Handle Nextflow list-passing artifacts (e.g. ['chrZ','chrW'])
sex_chromos = {t.strip(" [],'\"") for t in args.sex_chromos if t.strip(" [],'\"}")}

# ── Load per-chromosome coverage ──────────────────────────────────────────────

auto_depths = []
sex_depths  = []

for tsv_fn in args.chromo_cov_tsv_list:
    df = pd.read_csv(tsv_fn, sep="\t", header=0)
    for _, row in df.iterrows():
        chrom    = str(row["#rname"])
        depth    = float(row["meandepth"])
        scaf_len = int(row["endpos"])
        if chrom in sex_chromos:
            sex_depths.append(depth)
        elif scaf_len >= auto_chromos_min_length:
            auto_depths.append(depth)

# ── Calculate thresholds ──────────────────────────────────────────────────────

def calc_thresholds(depths):
    if not depths:
        return 5, 1000, float("nan")
    mean    = sum(depths) / len(depths)
    min_cov = max(5, math.floor(mean * min_mult))
    max_cov = math.ceil(mean * max_mult)
    return min_cov, max_cov, mean

auto_min, auto_max, auto_mean = calc_thresholds(auto_depths)
sex_min,  sex_max,  sex_mean  = calc_thresholds(sex_depths) if sex_depths \
                                 else (auto_min, auto_max, auto_mean)

# ── Write output ──────────────────────────────────────────────────────────────

out_fn = f"{indiv}_depth_thresholds.tsv"
with open(out_fn, "w") as fh:
    fh.write("type\tmin_cov\tmax_cov\tmean_cov\n")
    fh.write(f"auto\t{auto_min}\t{auto_max}\t"
             f"{auto_mean:.1f}\n" if auto_depths else f"auto\t{auto_min}\t{auto_max}\tNA\n")
    fh.write(f"sex\t{sex_min}\t{sex_max}\t"
             f"{sex_mean:.1f}\n" if sex_depths else f"sex\t{sex_min}\t{sex_max}\tNA\n")

# ── Print summary ─────────────────────────────────────────────────────────────

if auto_depths:
    print(f"[{indiv}] Autosome thresholds : min={auto_min}, max={auto_max}"
          f"  (n={len(auto_depths)} scaffolds >= {auto_chromos_min_length} bp,"
          f" mean={auto_mean:.1f}, multipliers={min_mult}x–{max_mult}x)")
else:
    print(f"[{indiv}] No autosome data (check auto_chromos_min_length={auto_chromos_min_length})")

if sex_depths:
    print(f"[{indiv}] Sex-chr thresholds  : min={sex_min}, max={sex_max}"
          f"  (n={len(sex_depths)} chromosomes, mean={sex_mean:.1f},"
          f" multipliers={min_mult}x–{max_mult}x)")
else:
    print(f"[{indiv}] No sex-chromosome data; using autosome thresholds")
