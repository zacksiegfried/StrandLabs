#!/usr/bin/env python3
"""
methyl_wide_formatting.py
------------------

Input:
  - data/methyl_data_trim_condensed.csv : condensed data with log_mean per marker per patient
        Columns: Study_Subject_ID, mdm, log_mean, log_sd,
                 CANCER_TYPE, cancer_yn, stage, AGE, SEX_D, log_cv


Output:
  - data/methyl_data_wide.csv : wide format with log_mean as columns per marker
        One row per patient with marker log_mean values as columns
"""

import pandas as pd
from pathlib import Path

def pivot_to_wide_format(input_path: str, output_path: str):

    df = pd.read_csv(input_path)


    # --- Data column checks ---
    required_cols = [
        "Study_Subject_ID", "mdm", "log_mean",
        "CANCER_TYPE", "cancer_yn", "stage", "AGE", "SEX_D"
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Metadata columns to preserve
    metadata_cols = ['CANCER_TYPE', 'cancer_yn', 'stage', 'AGE', 'SEX_D']


    # --- Pivot to wide format for log_mean (methylation markers) ---
    wide_df = df.pivot_table(
        index='Study_Subject_ID',
        columns='mdm',
        values='log_mean',
        aggfunc='first'  # Use first in case of duplicates
    ).reset_index()

    # Rename columns to indicate they are log values
    wide_df.columns.name = None
    log_cols = {col: f"{col}_log" for col in wide_df.columns if col != 'Study_Subject_ID'}
    wide_df = wide_df.rename(columns=log_cols)


    # --- Add metadata back (take first occurrence per patient) ---
    metadata = df[['Study_Subject_ID'] + metadata_cols].drop_duplicates(subset='Study_Subject_ID')
    wide_df = pd.merge(metadata, wide_df, on='Study_Subject_ID', how='right')

    # Reorder columns: metadata first, then marker values
    metadata_first = ['Study_Subject_ID'] + metadata_cols
    marker_cols = [col for col in wide_df.columns if col not in metadata_first]
    wide_df = wide_df[metadata_first + sorted(marker_cols)]

    # Fill missing CANCER_TYPE and stage with "non-cancer"
    wide_df['CANCER_TYPE'] = wide_df['CANCER_TYPE'].fillna('non-cancer')
    wide_df['stage'] = wide_df['stage'].fillna('non-cancer')

    # --- Write output ---
    output_file = Path(output_path) / "methyl_data_wide.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    wide_df.to_csv(output_file, index=False)
    print(f"✅ Saved to: {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Convert condensed methylation data to wide format with log_mean as column per marker.")
    parser.add_argument("--input", required=True, help="Path to methyl_data_trim_condensed.csv file")
    parser.add_argument("--output", default="data", help="Output folder path")
    args = parser.parse_args()

    pivot_to_wide_format(args.input, args.output)
