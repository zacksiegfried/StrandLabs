#!/usr/bin/env python3
"""
methyl_data_processing.py
------------------

Input:
    - data/clinical.csv : patient clinical metadata 
        Columns: subjid, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D
    - data/methylation.csv : patient methylation data, long format 3 reps per patient
        Columns: Study_Subject_ID, Study, mdm, log


Output:
    - data/methyl_data.csv or data/methyl_data_trim.csv : cleaned, merged table with all replicates retained
"""

import pandas as pd
from pathlib import Path

def methylation_merge(clinical_data_path: str, methylation_data_path: str, output_path: str, trimmed: bool = True):
    
    # Read data
    clinical_df = pd.read_csv(clinical_data_path)
    methylation_df = pd.read_csv(methylation_data_path)

    # Joining variable matching
    clinical_df['subjid'] = clinical_df['subjid'].astype(str).str.upper().str.strip()
    methylation_df['Study_Subject_ID'] = methylation_df['Study_Subject_ID'].astype(str).str.upper().str.strip()

    # Drop Study from clinical if it exists (avoid duplicate columns)
    clinical_df = clinical_df.drop(columns=['Study'], errors='ignore')

    # Merge
    merged_df = pd.merge(
        methylation_df,
        clinical_df,
        left_on='Study_Subject_ID',
        right_on='subjid',
        how='left'
    )

    # Drop unmatched rows
    before = len(merged_df)
    merged_df = merged_df.dropna(subset=['subjid'])
    print(f"Dropped {before - len(merged_df)} methylation records with no clinical match.")

    # Optional: trim to core columns
    if trimmed:
        core_cols = ['Study_Subject_ID', 'Study', 'mdm', 'log',
                     'CANCER_TYPE', 'cancer_yn', 'stage', 'AGE', 'SEX_D']
        merged_df = merged_df[[c for c in core_cols if c in merged_df.columns]]


    # --- Save output ---
    output_file = Path(output_path) / ("methyl_data_trim.csv" if trimmed else "methyl_data.csv")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_file, index=False)
    print(f"✅ Saved to {output_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Merge methylation and clinical datasets")
    parser.add_argument("--clinical", required=True, help="Path of clincial data")
    parser.add_argument("--methylation", required=True, help="Path of methylation data")
    parser.add_argument("--output", default="data", required=True, help="Path of output folder")
    parser.add_argument("--no-trim", dest="trimmed", action="store_false", default=True,
                        help="Disable trimming to output all columns (default: trimmed)")
    args = parser.parse_args()

    methylation_merge(args.clinical, args.methylation, args.output, args.trimmed)
