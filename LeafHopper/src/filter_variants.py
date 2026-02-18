#!/usr/bin/env python3
"""
Variant Filtering for LeafHopper Pipeline

Applies standardized filtering to input data before downstream analysis.
Outputs a filtered CSV containing only rows for variants that pass all
quality filters, ensuring consistent variant sets across analyses.

Filtering steps:
  1. Filter to specified donors
  2. Drop rows with null/empty mutid
  3. Compute per donor-variant-dilution stats (mean, std, count) on maf_col
  4. Remove groups with mean=0, std=0, or single observation
  5. Keep only variants with data from 2+ donors
  6. Output rows from original data matching surviving variants

Usage:
    python filter_variants.py --input data.csv --output results/ \
        --maf-col duplex_maf \
        --mutid-col mutid \
        --donor-col "Donor ID" \
        --dilution-col "Copy number" \
        --donors 417-1005 191-1055 Accugenomics
"""

import argparse
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter variants for downstream LeafHopper analyses"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input CSV file path"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--maf-col",
        type=str,
        required=True,
        help="Column name for MAF values used in quality filtering"
    )
    parser.add_argument(
        "--mutid-col",
        type=str,
        required=True,
        help="Column name for variant/mutation ID"
    )
    parser.add_argument(
        "--donor-col",
        type=str,
        required=True,
        help="Column name for donor/sample ID"
    )
    parser.add_argument(
        "--dilution-col",
        type=str,
        required=True,
        help="Column name for dilution level"
    )
    parser.add_argument(
        "--donors",
        type=str,
        nargs="+",
        required=True,
        help="Donors to include (required)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} rows")

    # Step 1: Filter to specified donors
    df = df[df[args.donor_col].isin(args.donors)]
    print(f"After donor filter: {len(df)} rows")

    # Step 2: Drop rows with null/empty mutid
    df = df[df[args.mutid_col].notna() & (df[args.mutid_col] != '')].copy()
    print(f"After mutid filter: {len(df)} rows ({df[args.mutid_col].nunique()} unique variants)")

    # Step 3: Compute per donor-variant-dilution stats on maf_col
    stats = df.groupby([args.donor_col, args.mutid_col, args.dilution_col]).agg(
        mean_val=(args.maf_col, 'mean'),
        std_val=(args.maf_col, 'std'),
        n_obs=(args.maf_col, 'count')
    ).reset_index()

    # Step 4: Remove groups with mean=0, std=0, or single observation
    stats = stats[(stats['mean_val'] > 0) & (stats['std_val'] > 0) & (stats['n_obs'] > 1)]
    remaining_variants = stats[args.mutid_col].unique()
    print(f"After quality filter: {len(remaining_variants)} variants with valid stats")

    # Step 5: Keep only variants with data from 2+ donors
    donors_per_variant = stats.groupby(args.mutid_col)[args.donor_col].nunique()
    multi_donor_variants = donors_per_variant[donors_per_variant >= 2].index.tolist()
    print(f"After multi-donor filter: {len(multi_donor_variants)} variants with 2+ donors")

    # Step 6: Filter original data to surviving variants
    df_filtered = df[df[args.mutid_col].isin(multi_donor_variants)].copy()
    print(f"Filtered output: {len(df_filtered)} rows")

    # Save filtered CSV
    filtered_path = args.output / "filtered_data.csv"
    df_filtered.to_csv(filtered_path, index=False)
    print(f"Saved: {filtered_path}")

    # Save variant list
    variants_path = args.output / "filtered_variants.txt"
    with open(variants_path, 'w') as f:
        for v in sorted(multi_donor_variants):
            f.write(f"{v}\n")
    print(f"Saved: {variants_path} ({len(multi_donor_variants)} variants)")


if __name__ == "__main__":
    main()
