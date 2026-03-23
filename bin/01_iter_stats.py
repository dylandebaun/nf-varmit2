#!/usr/bin/env python3
"""
01_iter_stats.py
Compute per-round reference-bias statistics from per-chromosome bgzipped VCFs
and a BAM file.  Statistics are reported separately for autosomes and sex
chromosomes.  A formatted summary is printed to stdout and a TSV is written
to the working directory.

Key metrics:
  mean_ab_het      – Mean allele balance at het sites; converges toward 0.5
                     as reference bias is eliminated across iterations.
  n_dominant_snps  – SNPs where ALT is dominant (AF ≥ 0.5); decreases as
                     the per-individual reference converges to the true sequence.
  titv_ratio       – Transitions/Transversions; sanity check (expect 1.5–2.5).
  pct_mapped       – Fraction of reads mapped to the current reference (BAM-level).
"""

import os
import sys
import argparse
import math
import subprocess
import statistics
import re

# ── Arguments ─────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-v", "--vcf_list", nargs="+", required=True,
                    help="Per-chromosome bgzipped VCF files")
parser.add_argument("-b", "--bam", required=True,
                    help="Sorted, indexed BAM file")
parser.add_argument("-i", "--individual", required=True,
                    help="Individual identifier")
parser.add_argument("-r", "--round", type=int, required=True,
                    help="Iteration round number")
parser.add_argument("-s", "--sex_chromos", nargs="+", default=[],
                    help="Sex chromosome names")
parser.add_argument("-l", "--auto_chromos_min_length", type=int, default=0,
                    help="Minimum scaffold length (bp) to count as autosome (default: 0)")
args = parser.parse_args()

indiv                   = args.individual
rnd                     = args.round
bam_fn                  = args.bam
vcf_list                = args.vcf_list
auto_chromos_min_length = args.auto_chromos_min_length

sex_chromos = {t.strip(" [],'\"") for t in args.sex_chromos if t.strip(" [],'\"}")}

# ── Scaffold lengths from BAM header ─────────────────────────────────────────

scaffold_lengths = {}
try:
    result = subprocess.run(["samtools", "view", "-H", bam_fn],
                            capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            if "SN" in fields and "LN" in fields:
                scaffold_lengths[fields["SN"]] = int(fields["LN"])
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Warning: could not read BAM header: {e}\n")

def is_valid_autosome(chrom):
    return scaffold_lengths.get(chrom, 0) >= auto_chromos_min_length

# ── Helper functions ──────────────────────────────────────────────────────────

TRANSITIONS = {frozenset(("A", "G")), frozenset(("C", "T"))}

def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and alt not in (".", "*")

def is_transition(ref, alt):
    return frozenset((ref.upper(), alt.upper())) in TRANSITIONS

def get_info_value(info_str, key):
    m = re.search(r"(?:^|;)" + re.escape(key) + r"=([^;]+)", info_str)
    return m.group(1).split(",")[0] if m else None

def empty_counters():
    return dict(n_variants=0, n_snps=0, n_het_snps=0, n_hom_snps=0,
                n_dominant_snps=0, ti=0, tv=0, ab_het=[])

# ── Parse VCFs ────────────────────────────────────────────────────────────────

counters = {"auto": empty_counters(), "sex": empty_counters()}

for vcf_fn in vcf_list:
    proc = subprocess.Popen(["bgzip", "-d", "-c", vcf_fn],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    for line in proc.stdout:
        line = line.rstrip("\n")
        if line.startswith("#") or len(line.split("\t")) < 8:
            continue

        cols      = line.split("\t")
        chrom     = cols[0]
        ref       = cols[3]
        alt_field = cols[4]
        info_str  = cols[7]

        if chrom in sex_chromos:
            ctype = "sex"
        elif is_valid_autosome(chrom):
            ctype = "auto"
        else:
            continue  # scaffold too short

        c = counters[ctype]
        c["n_variants"] += 1
        alt = alt_field.split(",")[0]

        try:
            af = float(get_info_value(info_str, "AF") or 0.0)
        except ValueError:
            af = 0.0
        try:
            ab = float(get_info_value(info_str, "AB")) if get_info_value(info_str, "AB") else None
        except ValueError:
            ab = None

        if not is_snp(ref, alt):
            continue

        c["n_snps"] += 1
        if af >= 0.5:
            c["n_dominant_snps"] += 1
        if is_transition(ref, alt):
            c["ti"] += 1
        else:
            c["tv"] += 1
        if ab is not None:
            if 0.2 <= ab <= 0.8:
                c["n_het_snps"] += 1
                c["ab_het"].append(ab)
            else:
                c["n_hom_snps"] += 1
    proc.wait()

# ── BAM-level mapping stats ───────────────────────────────────────────────────

pct_mapped = pct_properly_paired = float("nan")
try:
    fs = subprocess.run(["samtools", "flagstat", bam_fn],
                        capture_output=True, text=True, check=True)
    for line in fs.stdout.splitlines():
        if re.search(r"\bmapped\b", line) and "primary" not in line and "%" in line:
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                pct_mapped = float(m.group(1))
        if "properly paired" in line and "%" in line:
            m = re.search(r"\(([0-9.]+)%", line)
            if m:
                pct_properly_paired = float(m.group(1))
except subprocess.CalledProcessError as e:
    sys.stderr.write(f"Warning: samtools flagstat failed: {e}\n")

# ── Derived stats ─────────────────────────────────────────────────────────────

def derive(c):
    ab, nhom, tv = c["ab_het"], c["n_hom_snps"], c["tv"]
    return dict(
        n_variants      = c["n_variants"],
        n_snps          = c["n_snps"],
        n_het_snps      = c["n_het_snps"],
        n_hom_snps      = nhom,
        het_hom_ratio   = c["n_het_snps"] / nhom if nhom > 0 else float("nan"),
        mean_ab_het     = statistics.mean(ab)   if ab else float("nan"),
        median_ab_het   = statistics.median(ab) if ab else float("nan"),
        sd_ab_het       = (statistics.stdev(ab) if len(ab) > 1 else 0.0) if ab else float("nan"),
        n_dominant_snps = c["n_dominant_snps"],
        titv_ratio      = c["ti"] / tv if tv > 0 else float("nan"),
    )

derived = {ct: derive(counters[ct]) for ct in ("auto", "sex")}

# ── Write TSV ─────────────────────────────────────────────────────────────────

def fmt(v, d=4):
    return "nan" if (isinstance(v, float) and math.isnan(v)) else f"{v:.{d}f}"

HEADER = (
    "individual\tround\tchrom_type\t"
    "n_variants\tn_snps\tn_het_snps\tn_hom_snps\thet_hom_ratio\t"
    "mean_ab_het\tmedian_ab_het\tsd_ab_het\t"
    "n_dominant_snps\ttitv_ratio\t"
    "pct_mapped\tpct_properly_paired\n"
)

def tsv_row(ctype, d):
    return (
        f"{indiv}\t{rnd}\t{ctype}\t"
        f"{d['n_variants']}\t{d['n_snps']}\t{d['n_het_snps']}\t{d['n_hom_snps']}\t"
        f"{fmt(d['het_hom_ratio'], 3)}\t"
        f"{fmt(d['mean_ab_het'])}\t{fmt(d['median_ab_het'])}\t{fmt(d['sd_ab_het'])}\t"
        f"{d['n_dominant_snps']}\t{fmt(d['titv_ratio'], 3)}\t"
        f"{fmt(pct_mapped, 2)}\t{fmt(pct_properly_paired, 2)}\n"
    )

with open(f"{indiv}_round{rnd}_iter_stats.tsv", "w") as fh:
    fh.write(HEADER)
    fh.write(tsv_row("auto", derived["auto"]))
    if sex_chromos:
        fh.write(tsv_row("sex", derived["sex"]))

# ── Print human-readable summary ──────────────────────────────────────────────

def disp(v, d=4):
    if isinstance(v, float) and math.isnan(v):
        return "N/A"
    return f"{v:.{d}f}" if isinstance(v, float) else f"{v:,}"

def print_section(label, d):
    print(f"  [{label}]")
    print(f"  Variants called           : {disp(d['n_variants']):>10}")
    print(f"  SNPs                      : {disp(d['n_snps']):>10}")
    print(f"  Het SNPs  (AB 0.2–0.8)   : {disp(d['n_het_snps']):>10}")
    print(f"  Hom SNPs  (AB <0.2/>0.8) : {disp(d['n_hom_snps']):>10}")
    print(f"  Het/Hom ratio             : {disp(d['het_hom_ratio'], 3):>10}")
    print(f"  Dominant SNPs (AF≥0.5)   : {disp(d['n_dominant_snps']):>10}"
          f"  (decreases as reference converges)")
    print(f"  {'─'*54}")
    print(f"  Allele balance at het sites (target → 0.5 = unbiased):")
    print(f"    Mean AB                 : {disp(d['mean_ab_het']):>10}")
    print(f"    Median AB               : {disp(d['median_ab_het']):>10}")
    print(f"    SD AB                   : {disp(d['sd_ab_het']):>10}")
    print(f"  {'─'*54}")
    print(f"  Ti/Tv ratio               : {disp(d['titv_ratio'], 3):>10}  (expected 1.5–2.5)")

print(f"\n{'='*60}")
print(f"  Reference bias stats  |  {indiv}  |  Round {rnd}"
      f"  |  auto min len: {auto_chromos_min_length} bp")
print(f"{'='*60}")
print_section("Autosomes", derived["auto"])
if sex_chromos:
    print(f"{'─'*60}")
    print_section("Sex chromosomes", derived["sex"])
print(f"{'─'*60}")
print(f"  Reads mapped (%)          : {disp(pct_mapped, 2):>10}  (BAM-level)")
print(f"  Properly paired (%)       : {disp(pct_properly_paired, 2):>10}  (BAM-level)")
print(f"{'='*60}\n")
