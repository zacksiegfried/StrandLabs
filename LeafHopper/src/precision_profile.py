#!/usr/bin/env python3
"""
Precision Profile Plot: %CV vs Mean SDMC by Variant

Generates one scatter plot per variant showing coefficient of variation
vs mean count, colored by Donor ID. Each point represents a group level
(e.g., Copy number) within a donor.

Usage:
    python precision_profile.py --input data.csv --output results/ \
        --sdmc-col super_duplex_mutant_count \
        --mutid-col mutid \
        --donor-col "Donor ID" \
        --group-col "Copy number" \
        --donors 417-1005 191-1055 Accugenomics
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path


def power_func(x, a, b):
    """Power function: y = a * x^b"""
    return a * np.power(x, b)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate precision profile plots (%CV vs mean SDMC) - one per variant"
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
        "--sdmc-col",
        type=str,
        default="super_duplex_mutant_count",
        help="Column name for SDMC values (default: super_duplex_mutant_count)"
    )
    parser.add_argument(
        "--mutid-col",
        type=str,
        default="mutid",
        help="Column name for variant/mutation ID (default: mutid)"
    )
    parser.add_argument(
        "--donor-col",
        type=str,
        default="Donor ID",
        help="Column name for donor/sample ID (default: 'Donor ID')"
    )
    parser.add_argument(
        "--group-col",
        type=str,
        default="Copy number",
        help="Column name for grouping variable (default: 'Copy number')"
    )
    parser.add_argument(
        "--donors",
        type=str,
        nargs="+",
        required=True,
        help="Donors to include (required)"
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="precision",
        help="Output file prefix (default: precision)"
    )
    return parser.parse_args()


def compute_precision_stats(df: pd.DataFrame, sdmc_col: str, mutid_col: str,
                            donor_col: str, group_col: str) -> pd.DataFrame:
    """
    Compute mean and %CV for each variant within each donor and group level.

    Returns DataFrame with columns: donor, mutid, group, mean_sdmc, std_sdmc, cv_percent, n_obs
    """
    stats = df.groupby([donor_col, mutid_col, group_col]).agg({
        sdmc_col: ['mean', 'std', 'count']
    }).reset_index()

    stats.columns = ['donor', 'mutid', 'group', 'mean_value', 'std_value', 'n_obs']

    # Calculate %CV (coefficient of variation)
    stats['cv_percent'] = (stats['std_value'] / stats['mean_value']) * 100

    # Remove rows with zero mean, zero std, or single observation
    stats = stats[(stats['mean_value'] > 0) & (stats['std_value'] > 0) & (stats['n_obs'] > 1)].copy()

    return stats


def plot_variant_precision(stats: pd.DataFrame, variant: str, output_path: Path,
                           donors: list, quant_col: str):
    """Generate precision profile scatter plot for a single variant."""

    variant_stats = stats[stats['mutid'] == variant]

    if len(variant_stats) == 0:
        return False

    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {'417-1005': 'blue', '191-1055': 'green', 'Accugenomics': 'red'}
    default_colors = plt.cm.tab10(np.linspace(0, 1, len(donors)))

    # Get x range for fitting curves
    x_min = variant_stats['mean_value'].min()
    x_max = variant_stats['mean_value'].max()
    x_line = np.linspace(max(x_min, 0.1), x_max, 100)

    for i, donor in enumerate(donors):
        donor_stats = variant_stats[variant_stats['donor'] == donor]
        if len(donor_stats) == 0:
            continue

        color = colors.get(donor, default_colors[i])

        # Plot scatter points
        ax.scatter(
            donor_stats['mean_value'],
            donor_stats['cv_percent'],
            label=f"{donor} (n={len(donor_stats)})",
            alpha=0.7,
            s=20,
            color=color
        )

        # Fit and plot power curve y = a * x^b
        if len(donor_stats) >= 3:
            x_data = donor_stats['mean_value'].values
            y_data = donor_stats['cv_percent'].values

            # Filter out zeros/negatives for fitting
            mask = (x_data > 0) & (y_data > 0)
            if mask.sum() >= 3:
                try:
                    popt, _ = curve_fit(power_func, x_data[mask], y_data[mask],
                                        p0=[100, -0.5], maxfev=5000)
                    x_fit = np.linspace(x_data[mask].min(), x_data[mask].max(), 100)
                    y_fit = power_func(x_fit, *popt)
                    ax.plot(x_fit, y_fit, color=color, linestyle='--', alpha=0.6, linewidth=1.5)
                except (RuntimeError, ValueError):
                    pass  # Skip curve if fitting fails

    # Format column name for display (replace underscores with spaces, title case)
    x_label = f"Mean {quant_col.replace('_', ' ').title()}"
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel('Coefficient of Variation (%)', fontsize=12)
    ax.set_title(f'Precision Profile: {variant}', fontsize=14)
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

    # Load data
    df = pd.read_csv(args.input)

    # Step 1: Filter to specified donors
    df = df[df[args.donor_col].isin(args.donors)]
    print(f"After donor filter: {len(df)} rows")

    # Step 2: Get unique variants from mutid column (non-null only)
    df = df[df[args.mutid_col].notna() & (df[args.mutid_col] != '')].copy()
    unique_variants = df[args.mutid_col].unique()
    print(f"Found {len(unique_variants)} unique variants in {args.mutid_col}")
    print(f"Donors: {df[args.donor_col].unique().tolist()}")

    # Compute stats (grouped by donor, variant, and group level)
    stats = compute_precision_stats(df, args.sdmc_col, args.mutid_col,
                                    args.donor_col, args.group_col)
    print(f"Computed stats for {len(stats)} donor-variant-group combinations")

    # Filter to variants with data from 2+ donors
    donors_per_variant = stats.groupby('mutid')['donor'].nunique()
    multi_donor_variants = donors_per_variant[donors_per_variant >= 2].index.tolist()
    print(f"Found {len(multi_donor_variants)} variants with data from 2+ donors")

    # Generate one plot per variant (only multi-donor variants)
    print(f"Generating plots...")
    plot_count = 0
    for variant in multi_donor_variants:
        # Clean variant name for filename
        safe_name = str(variant).replace(':', '_').replace('/', '_').replace(' ', '_')
        plot_path = plots_dir / f"{args.prefix}_{safe_name}.png"

        if plot_variant_precision(stats, variant, plot_path, args.donors, args.sdmc_col):
            plot_count += 1

    print(f"Saved {plot_count} plots to {plots_dir}")

    # Save stats
    stats_path = args.output / f"{args.prefix}_stats.csv"
    stats.to_csv(stats_path, index=False)
    print(f"Saved: {stats_path}")


if __name__ == "__main__":
    main()
