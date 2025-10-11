#!/usr/bin/env python3
"""
data_processing.py
------------------

Input:
  - data/clinical.csv : patient clinical metadata 
        Columns: subjid, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D
  - data/methylation.csv : patient methylation data, long format 3 reps per patient
        Columns: Study_Subject_ID, Study, mdm, log, strands

Output:
  - data/merged_data.csv : cleaned, merged table with all replicates retained
"""

import pandas as pd
from pathlib import Path

def methylation_merge(clinical_data_path: str, methylation_data_path: str, output_path: str, trimmed: bool = True):

    # --- Convert output path to Path object ---
    out_folder = Path(output_path)
    out_folder.mkdir(parents=True, exist_ok=True)  # ensure folder exists

    clinical_df = pd.read_csv(clinical_data_path)
    methylation_df = pd.read_csv(methylation_data_path)


    # --- Data column checks ---
    if 'subjid' not in clinical_df.columns:
        raise ValueError("incorrect clinical data format")

    if 'Study_Subject_ID' not in methylation_df.columns:
        raise ValueError("incorrect methylation data format")
    
    # Joining variable matching
    clinical_df['subjid'] = clinical_df['subjid'].astype(str).str.upper().str.strip()
    methylation_df['Study_Subject_ID'] = methylation_df['Study_Subject_ID'].astype(str).str.upper().str.strip()

    # Drop Study from clinical_df if it exists
    if "Study" in clinical_df.columns:
        clinical_df = clinical_df.drop(columns=["Study"])

    merged_df = pd.merge(
        methylation_df,
        clinical_df,
        left_on='Study_Subject_ID',
        right_on='subjid',
        how='left',
        validate='many_to_one'  # ensures methylation→clinical is many-to-one
    )

    # Drop methylation rows with no matching clinical data
    before = len(merged_df)
    merged_df = merged_df.dropna(subset=['subjid'])
    after = len(merged_df)
    print(f"Dropped {before - after} methylation records with no clinical match.")


    # --- Save output ---
    if trimmed:
        filename = "merged_data_trim.csv"
        # Trim columns to only the core set
        core_columns = [
            "Study_Subject_ID", "Study", "mdm", "log", "strands",
            "CANCER_TYPE", "cancer_yn", "stage", "AGE", "SEX_D"
        ]
        merged_df = merged_df[[c for c in core_columns if c in merged_df.columns]]
    else:
        filename = "merged_data.csv"

    # --- Save output ---
    out_file = out_folder / filename
    merged_df.to_csv(out_file, index=False)
    print(f"✅ Output saved to {out_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Merge methylation and clinical datasets")
    parser.add_argument("--clinical", required=True, help="Path to clinical metadata CSV")
    parser.add_argument("--methylation", required=True, help="Path to methylation data CSV")
    parser.add_argument("--output", required=True, help="Path to save merged output CSV")
    parser.add_argument("--trimmed", action="store_true", help="Output trimmed dataset")
    args = parser.parse_args()

    methylation_merge(args.clinical, args.methylation, args.output, args.trimmed)
