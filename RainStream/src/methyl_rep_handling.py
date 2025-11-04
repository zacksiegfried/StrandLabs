#!/usr/bin/env python3
"""
methyl_rep_handling.py
------------------

Input:
  - data/methyl_data_trim.csv : merged data set from methyl_data_processing.py output
        Columns: Study_Subject_ID, Study, mdm, log, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D


Output:
  - data/methyl_data_trim_condensed.csv : replicates averaged to return mean, sd, %CV for each marker for each patient
"""

import pandas as pd
from pathlib import Path

def summarize_replicates(input_path: str, output_path: str):

    df = pd.read_csv(input_path)


    # --- Data column checks ---
    required_cols = [
        "Study_Subject_ID", "Study", "mdm", "log",
        "CANCER_TYPE", "cancer_yn", "stage", "AGE", "SEX_D"
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

     # Metadata columns to keep
    metadata_cols = ['CANCER_TYPE', 'cancer_yn', 'stage', 'AGE', 'SEX_D']
    

    # --- Summary stats calculation ---
    condensed = df.groupby(['Study_Subject_ID', 'mdm']).agg(
        log_mean=('log', 'mean'),
        log_sd=('log', 'std'),
        **{col: (col, 'first') for col in metadata_cols}
    ).reset_index()

    # Calculate %CV
    condensed['log_cv'] = (condensed['log_sd'] / condensed['log_mean']) * 100


    # --- Write output ---
    output_file = Path(output_path) / "methyl_data_trim_condensed.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    condensed.to_csv(output_file, index=False)
    print(f"✅ Saved to: {output_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Condense replicate methylation data to mean, SD, and %CV per patient.")
    parser.add_argument("--input", required=True, help="Path to methyl_data_trim.csv file")
    parser.add_argument("--output", default="data", help="Output CSV path")
    args = parser.parse_args()

    summarize_replicates(args.input, args.output)
