#!/usr/bin/env python3
"""
Hit Rate LoD (Limit of Detection) Analysis

Calculates detection (hit) rate across dilution levels per variant and donor.
A "hit" is defined as alt_n_col >= threshold. Generates one plot per variant
showing hit rate vs dilution level, colored by Donor ID.

Usage:
    python hit_rate_lod.py --input data.csv --output results/ \
        --alt-n-col super_duplex_mutant_count \
        --mutid-col mutid \
        --donor-col "Donor ID" \
        --dilution-col "Copy number" \
        --donors 417-1005 191-1055 Accugenomics \
        --hit-threshold 1
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate hit rate LoD plots - one per variant"
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
        "--alt-n-col",
        type=str,
        required=True,
        help="Column name for allele count values (e.g. super_duplex_mutant_count)"
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
    parser.add_argument(
        "--hit-threshold",
        type=float,
        default=1,
        help="Minimum allele count to qualify as a hit (default: 1)"
    )
    return parser.parse_args()


def compute_hit_rate(df: pd.DataFrame, alt_n_col: str, mutid_col: str,
                     donor_col: str, dilution_col: str,
                     threshold: float) -> pd.DataFrame:
    """
    Compute hit rate for each variant within each donor and dilution level.

    A hit is defined as alt_n_col >= threshold. The denominator (n_reps) for
    each donor-variant group is taken from the dilution level with the most
    observations, since undetected replicates at low dilutions produce no row.
    Missing dilution levels are filled in with 0 hits.

    Returns DataFrame with columns: donor, mutid, dilution, n_reps, n_hits, hit_rate
    """
    df = df.copy()
    df['_hit'] = (df[alt_n_col] >= threshold).astype(int)

    # Count hits per donor-variant-dilution
    hits = df.groupby([donor_col, mutid_col, dilution_col]).agg(
        n_obs=('_hit', 'count'),
        n_hits=('_hit', 'sum')
    ).reset_index()

    hits.columns = ['donor', 'mutid', 'dilution', 'n_obs', 'n_hits']

    # For each donor-variant group, the true replicate count is the max
    # observed across all dilution levels (highest dilutions have all reps)
    n_reps = hits.groupby(['donor', 'mutid'])['n_obs'].max().reset_index()
    n_reps.columns = ['donor', 'mutid', 'n_reps']
    hits = hits.merge(n_reps, on=['donor', 'mutid'])

    # Ensure all dilution levels are represented for each donor-variant group
    all_dilutions = sorted(hits['dilution'].unique())
    full_index = pd.MultiIndex.from_product(
        [hits[['donor', 'mutid']].drop_duplicates().itertuples(index=False, name=None),
         all_dilutions],
        names=['donor_mutid', 'dilution']
    )
    # Build a proper full grid
    donor_mutid_pairs = hits[['donor', 'mutid']].drop_duplicates()
    rows = []
    for _, row in donor_mutid_pairs.iterrows():
        for dil in all_dilutions:
            rows.append({'donor': row['donor'], 'mutid': row['mutid'], 'dilution': dil})
    full_grid = pd.DataFrame(rows)

    stats = full_grid.merge(hits, on=['donor', 'mutid', 'dilution'], how='left')
    stats['n_hits'] = stats['n_hits'].fillna(0).astype(int)

    # Fill n_reps from the group-level value
    stats = stats.drop(columns=['n_obs', 'n_reps'], errors='ignore')
    stats = stats.merge(n_reps, on=['donor', 'mutid'])

    stats['hit_rate'] = stats['n_hits'] / stats['n_reps']

    return stats


def plot_variant_hit_rate(stats: pd.DataFrame, variant: str, output_path: Path,
                          donors: list, threshold: float):
    """Generate hit rate plot for a single variant."""

    variant_stats = stats[stats['mutid'] == variant].copy()

    if len(variant_stats) == 0:
        return False

    # Sort by dilution for line plotting
    variant_stats = variant_stats.sort_values('dilution')

    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {'417-1005': 'blue', '191-1055': 'green', 'Accugenomics': 'red'}
    default_colors = plt.cm.tab10(np.linspace(0, 1, len(donors)))

    for i, donor in enumerate(donors):
        donor_stats = variant_stats[variant_stats['donor'] == donor]
        if len(donor_stats) == 0:
            continue

        color = colors.get(donor, default_colors[i])

        ax.plot(
            donor_stats['dilution'],
            donor_stats['hit_rate'],
            marker='o',
            markersize=5,
            label=donor,
            color=color,
            alpha=0.8,
            linewidth=1.5
        )

    ax.set_xlabel('Dilution Level', fontsize=12)
    ax.set_ylabel('Hit Rate', fontsize=12)
    ax.set_title(f'Hit Rate LoD: {variant} (threshold >= {threshold})', fontsize=14)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(title='Donor ID')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return True


def main():
    args = parse_args()

    # Create output directories
    args.output.mkdir(parents=True, exist_ok=True)
    plots_dir = args.output / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Load pre-filtered data
    df = pd.read_csv(args.input)
    variants = df[args.mutid_col].unique()
    print(f"Loaded {len(df)} rows, {len(variants)} variants")

    # Ensure dilution column is numeric
    df[args.dilution_col] = pd.to_numeric(df[args.dilution_col], errors='coerce')
    df = df.dropna(subset=[args.dilution_col])

    # Compute hit rates
    stats = compute_hit_rate(df, args.alt_n_col, args.mutid_col,
                             args.donor_col, args.dilution_col,
                             args.hit_threshold)
    print(f"Computed hit rates for {len(stats)} donor-variant-dilution combinations")

    # Generate one plot per variant
    print(f"Generating plots...")
    plot_count = 0
    for variant in variants:
        safe_name = str(variant).replace(':', '_').replace('/', '_').replace(' ', '_')
        plot_path = plots_dir / f"lod_{safe_name}.png"

        if plot_variant_hit_rate(stats, variant, plot_path, args.donors, args.hit_threshold):
            plot_count += 1

    print(f"Saved {plot_count} plots to {plots_dir}")

    # Save stats
    stats_path = args.output / "hit_rate_stats.csv"
    stats.to_csv(stats_path, index=False)
    print(f"Saved: {stats_path}")


if __name__ == "__main__":
    main()
