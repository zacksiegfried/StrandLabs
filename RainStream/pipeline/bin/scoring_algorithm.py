#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from scipy.stats import binom

def pscore_algorithm(input_file, output_dir):
    # Read the input file in .csv format and get file name
    df = pd.read_csv(input_file)
    input_basename = os.path.basename(input_file)

    # --- Your processing logic here ---
    df['new_column'] = 'processed'
    #df["p_score"] = binom.pmf()


    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    # Create output file path
    output_file_path = os.path.join(output_dir, f'{input_basename}_pscore.csv')
    
    # Save the result
    df.to_csv(output_file, index=False)
    print(f"Saved processed file to: {output_file_path}")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input', '-i', required=True, help="")
    parser.add_argument('--output', '-o', required=True, help="")

    args = parser.parse_args()

    pscore_algorithm(args.input, args.output)

if __name__ == '__main__':
    main()