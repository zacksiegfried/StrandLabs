#!/usr/bin/env python3
"""
methyl_precision_profile.py
---------------------------

Computes within-replicate precision (%CV) for every patient × marker and plots
a precision profile for each marker: %CV (y) vs. mean log strands (x).

Input:
  - methyl_data_trim.csv : merged data with up to 3 replicates per patient
        Columns: Study_Subject_ID, Study, mdm, log, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D

Output:
  - methyl_precision_profile_data.csv  : per-patient × marker precision stats
        Columns: Study_Subject_ID, mdm, n_reps, log_mean, log_sd, cv_pct
  - marker_precision_profiles/*.png (optional) : per-marker precision profile plots
"""

import pandas as pd
import numpy as np
from pathlib import Path


def compute_precision_profile(
    input_path: str,
    output_path: str,
    plot: bool = False,
) -> None:

    df = pd.read_csv(input_path)

    required_cols = ["Study_Subject_ID", "Study", "mdm", "log"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    print(f"Loaded {len(df):,} rows | {df['mdm'].nunique()} markers | {df['Study_Subject_ID'].nunique()} patients")

    # --- Within-patient replicate precision per marker ---
    rep_stats = (
        df.groupby(["Study_Subject_ID", "mdm"])["log"]
        .agg(n_reps="count", log_mean="mean", log_sd="std")
        .reset_index()
    )

    # CV is undefined when mean == 0 (below-limit-of-detection); exclude
    n_zero = (rep_stats["log_mean"] == 0).sum()
    precision = rep_stats[rep_stats["log_mean"] > 0].copy()
    precision["cv_pct"] = (precision["log_sd"] / precision["log_mean"]) * 100

    print(f"Excluded {n_zero:,} patient×marker observations with log_mean = 0")
    print(f"Retained {len(precision):,} observations for precision profiling")

    # --- Save per-patient precision data ---
    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_out = out_dir / "methyl_precision_profile_data.csv"
    precision.round(4).to_csv(csv_out, index=False)
    print(f"✅ Saved precision data ({len(precision):,} rows) → {csv_out}")

    if not plot:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    CV_THRESHOLD = 25.0

    def fit_log_linear(x, y):
        # Fit log(%CV) ~ mean_log_strands, back-transform to %CV scale
        # Returns fitted curve (fx, fy) and predicted x where %CV crosses CV_THRESHOLD
        mask = y > 0
        if mask.sum() < 4:
            return None, None, None
        slope, intercept = np.polyfit(x[mask], np.log(y[mask]), 1)
        x_line = np.linspace(x.min(), x.max(), 200)
        y_line = np.exp(intercept + slope * x_line)
        # Solve exp(intercept + slope * x) = CV_THRESHOLD → x = (log(CV_THRESHOLD) - intercept) / slope
        x_pred = (np.log(CV_THRESHOLD) - intercept) / slope if slope != 0 else None
        return x_line, y_line, x_pred

    profiles_dir = out_dir / "marker_precision_profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    markers = sorted(precision["mdm"].unique())

    for marker in markers:
        sub = precision[precision["mdm"] == marker]
        x = sub["log_mean"].values
        y = sub["cv_pct"].values

        fig, ax = plt.subplots(figsize=(6, 4))

        ax.scatter(x, y, s=18, alpha=0.5, color="#2196F3", linewidths=0)

        fx, fy, x_pred = fit_log_linear(x, y)
        if fx is not None:
            ax.plot(fx, fy, color="#E53935", linewidth=1.8, label="log(%CV) ~ mean log")
            ax.legend(fontsize=8, frameon=False)

        # Horizontal reference line at 25% CV
        ax.axhline(CV_THRESHOLD, color="grey", linewidth=1.0, linestyle="--")

        if x_pred is not None:
            ax.text(
                0.97, 0.5, f"25% CV = {x_pred:.2f}",
                transform=ax.transAxes, fontsize=8,
                ha="right", va="center", color="black",
            )

        ax.set_title(marker, fontsize=11)
        ax.set_xlabel("Mean log strands", fontsize=9)
        ax.set_ylabel("%CV", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        plt.tight_layout()

        safe_name = marker.replace("/", "_").replace("(", "").replace(")", "")
        plot_out = profiles_dir / f"{safe_name}_precision_profile.png"
        fig.savefig(plot_out, dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"✅ Saved {len(markers)} precision profile plots → {profiles_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Precision profile: %CV vs. mean log strands for all markers."
    )
    parser.add_argument("--input",  required=True, help="Path to methyl_data_trim.csv")
    parser.add_argument("--output", default="output", help="Output directory")
    parser.add_argument("--plot",   action="store_true", help="Generate precision profile plots")
    args = parser.parse_args()

    compute_precision_profile(args.input, args.output, plot=args.plot)
