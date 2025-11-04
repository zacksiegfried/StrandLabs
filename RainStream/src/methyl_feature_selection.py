#!/usr/bin/env python3
"""
feature_selection.py
------------------

Input:
  - data/methyl_data_wide.csv : wide format methylation data with cancer labels
        Columns: Study_Subject_ID, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D, [marker]_log columns

Output:
  - data/feature_importance_scores.csv : ranked features with importance scores
  - data/top_markers_plot.png : visualization of top markers (optional)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.impute import SimpleImputer


def analyze_features(
    input_path: str,
    output_path: str,
    n_top_features: int = 20,
    save_plot: bool = False
):
    """
    Perform comprehensive feature importance analysis on methylation markers.

    Uses three methods:
    1. ANOVA F-test (univariate statistical test)
    2. Random Forest feature importance (model-based)
    3. Logistic Regression coefficients (linear model)

    Args:
        input_path: Path to methyl_data_wide.csv
        output_path: Directory to save output files
        n_top_features: Number of top features to report
        save_plot: Whether to generate and save visualization
    """

    # --- Load data ---
    df = pd.read_csv(input_path)

    # Extract methylation marker columns
    marker_cols = [col for col in df.columns if col.endswith('_log')]

    X = df[marker_cols]
    y = df['cancer_yn']
    print(f"Loaded {len(df)} samples with {len(marker_cols)} methylation markers")

    # --- Handle missing values --- (OFF)
    #imputer = SimpleImputer(strategy='median')
    #X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=marker_cols)
    #missing_count = X.isnull().sum().sum()
    #if missing_count > 0:
    #    print(f"Imputed {missing_count} missing values using median strategy")


    # --- Method 1: ANOVA F-test ---
    selector = SelectKBest(f_classif, k='all')
    selector.fit(X, y)

    f_scores = selector.scores_
    p_values = selector.pvalues_


    # --- Method 2: Random Forest importance ---
    rf = RandomForestClassifier(
        n_estimators=200,
        random_state=999,
        max_depth=10,
        n_jobs=-1  # Use all cores
    )
    rf.fit(X, y)
    rf_importance = rf.feature_importances_


    # --- Method 3: Logistic Regression coefficients ---
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X_scaled, y)
    lr_coef = np.abs(lr.coef_[0])  # Absolute value for importance ranking


    # --- Combine results ---
    results = pd.DataFrame({
        'marker': [col.replace('_log', '') for col in marker_cols],
        'marker_full': marker_cols,
        'f_score': f_scores,
        'p_value': p_values,
        'rf_importance': rf_importance,
        'lr_coefficient': lr_coef
    })

    # Normalize scores to 0-1 range for fair comparison
    norm_scaler = MinMaxScaler()
    results[['f_score_norm', 'rf_norm', 'lr_norm']] = norm_scaler.fit_transform(
        results[['f_score', 'rf_importance', 'lr_coefficient']]
    )

    # Consensus score: average of normalized scores
    results['consensus_score'] = results[['f_score_norm', 'rf_norm', 'lr_norm']].mean(axis=1)

    # Sort by consensus score
    results = results.sort_values('consensus_score', ascending=False)

    # Add rank column
    results.insert(0, 'rank', range(1, len(results) + 1))


    # --- Print top features ---
    print(f"\n{'='*80}")
    print(f"Top {n_top_features} Methylation Markers (by consensus score):")
    print('='*80)
    print(results.head(n_top_features).to_string(index=False))


    # --- Save results ---
    output_dir = Path(output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "feature_importance_scores.csv"
    results.to_csv(output_file, index=False)
    print(f"\n✅ Feature importance scores saved to: {output_file}")


    # --- Optional: Generate plot ---
    if save_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend for server use
            import matplotlib.pyplot as plt

            # Create visualization
            top_features = results.head(n_top_features)

            fig, axes = plt.subplots(1, 3, figsize=(18, 6))

            # F-score plot
            axes[0].barh(top_features['marker'], top_features['f_score'], color='steelblue')
            axes[0].set_xlabel('F-Score')
            axes[0].set_ylabel('Methylation Marker')
            axes[0].set_title('ANOVA F-Test')
            axes[0].invert_yaxis()

            # Random Forest importance plot
            axes[1].barh(top_features['marker'], top_features['rf_importance'], color='forestgreen')
            axes[1].set_xlabel('Importance')
            axes[1].set_title('Random Forest')
            axes[1].invert_yaxis()

            # Logistic Regression coefficient plot
            axes[2].barh(top_features['marker'], top_features['lr_coefficient'], color='coral')
            axes[2].set_xlabel('|Coefficient|')
            axes[2].set_title('Logistic Regression')
            axes[2].invert_yaxis()

            plt.suptitle(f'Top {n_top_features} Methylation Markers by Different Methods',
                        fontsize=14, fontweight='bold')
            plt.tight_layout()

            plot_file = output_dir / "top_markers_plot.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"✅ Visualization saved to: {plot_file}")

        except ImportError:
            print("⚠️  matplotlib not available, skipping plot generation")


    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Perform feature importance analysis on methylation markers for cancer detection"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to methyl_data_wide.csv file"
    )
    parser.add_argument(
        "--output",
        default="output",
        help="Output directory path (default: output)"
    )
    parser.add_argument(
        "--n-features",
        type=int,
        default=20,
        help="Number of top features to display (default: 20)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate and save visualization plot"
    )

    args = parser.parse_args()

    analyze_features(
        input_path=args.input,
        output_path=args.output,
        n_top_features=args.n_features,
        save_plot=args.plot
    )
