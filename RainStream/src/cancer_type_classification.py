#!/usr/bin/env python3
"""
cancer_type_classification.py
------------------

Input:
  - data/methyl_data_wide.csv : wide format methylation data with cancer types
        Columns: Study_Subject_ID, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D, [marker]_log columns

Output:
  - data/cancer_type_predictions.csv : classification results and predictions
  - data/cancer_type_model_summary.csv : model performance metrics
  - data/cancer_type_confusion_matrix.png : visualization (optional)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, List
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import warnings
warnings.filterwarnings('ignore')


def classify_cancer_types(
    input_path: str,
    output_path: str,
    cancer_types: Optional[List[str]] = None,
    test_size: float = 0.25,
    top_n_features: Optional[int] = None,
    save_plot: bool = False
):
    """
    Train multi-class classifier to predict cancer type from methylation markers.
    Allows specification of which cancer types to include.

    Args:
        input_path: Path to methyl_data_wide.csv
        output_path: Directory to save output files
        cancer_types: List of specific cancer types to include (default: ['Lung', 'Pancreatic'])
        test_size: Fraction of data for testing (default: 0.25)
        top_n_features: If specified, use only top N features (feature selection)
        save_plot: Whether to generate confusion matrix plot
    """

    # Default to Lung and Pancreatic (known to have distinct profiles)
    if cancer_types is None:
        cancer_types = ['Lung', 'Pancreatic']

    print("="*80)
    print("CANCER TYPE CLASSIFICATION")
    print("="*80)

    # --- Load data ---
    df = pd.read_csv(input_path)

    # Filter to only cancer patients
    cancer_df = df[df['cancer_yn'] == 1].copy()
    print(f"\nLoaded {len(cancer_df)} cancer samples")

    # Extract features
    marker_cols = [col for col in df.columns if col.endswith('_log')]

    print(f"Features: {len(marker_cols)} methylation markers")


    # --- Show all available cancer types ---
    cancer_counts = cancer_df['CANCER_TYPE'].value_counts()
    print(f"\nAvailable cancer types in dataset:")
    for cancer_type, count in cancer_counts.items():
        print(f"  - {cancer_type:25} ({count} samples)")


    # --- Filter to specified cancer types ---
    # Check if requested cancer types exist in data
    available_cancers = set(cancer_df['CANCER_TYPE'].unique())
    requested_cancers = set(cancer_types)
    missing_cancers = requested_cancers - available_cancers

    if missing_cancers:
        print(f"\n⚠️  Warning: The following cancer types were not found in the data:")
        for cancer in missing_cancers:
            print(f"     - {cancer}")
        # Use only the ones that exist
        cancer_types = list(requested_cancers & available_cancers)

    if not cancer_types:
        raise ValueError("No valid cancer types found in the dataset!")

    print(f"\n{'='*80}")
    print(f"Selected Cancer Types for Classification:")
    print("="*80)
    for i, cancer in enumerate(cancer_types, 1):
        count = cancer_counts.get(cancer, 0)
        print(f"  {i}. {cancer:25} ({count} samples)")

    # Filter dataset
    cancer_df_filtered = cancer_df[cancer_df['CANCER_TYPE'].isin(cancer_types)].copy()
    print(f"\nFiltered dataset: {len(cancer_df_filtered)} samples")

    if len(cancer_df_filtered) < 10:
        print(f"\n⚠️  Warning: Very small dataset ({len(cancer_df_filtered)} samples)")
        print("   Model performance may be unreliable with so few samples.")

    X = cancer_df_filtered[marker_cols]
    y = cancer_df_filtered['CANCER_TYPE']


    # --- Handle missing values ---
    imputer = SimpleImputer(strategy='median')
    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=marker_cols, index=X.index)

    missing_count = X.isnull().sum().sum()
    if missing_count > 0:
        print(f"Imputed {missing_count} missing values using median strategy")


    # --- Optional: Feature selection ---
    if top_n_features:
        from sklearn.feature_selection import SelectKBest, f_classif

        print(f"\nPerforming feature selection (top {top_n_features} markers)...")
        selector = SelectKBest(f_classif, k=min(top_n_features, len(marker_cols)))
        X_selected = selector.fit_transform(X_imputed, y)
        selected_features = X_imputed.columns[selector.get_support()].tolist()

        print(f"Selected top {len(selected_features)} features:")
        print(f"  {', '.join([f.replace('_log', '') for f in selected_features[:10]])}...")
        X_imputed = pd.DataFrame(X_selected, columns=selected_features, index=X.index)


    # --- Encode labels ---
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    class_names = le.classes_

    print(f"\n{'='*80}")
    print("Class Encoding:")
    print("="*80)
    for i, name in enumerate(class_names):
        count = (y == name).sum()
        print(f"  {i}: {name:20} ({count} samples)")


    # --- Train/test split ---
    # Use stratified split to maintain class distribution
    X_train, X_test, y_train, y_test = train_test_split(
        X_imputed, y_encoded,
        test_size=test_size,
        random_state=42,
        stratify=y_encoded
    )

    print(f"\nTrain/Test split:")
    print(f"  Training: {len(X_train)} samples")
    print(f"  Testing: {len(X_test)} samples")


    # --- Scale features ---
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)


    # --- Train models ---
    print("\n" + "="*80)
    print("TRAINING MODELS")
    print("="*80)

    models = {}

    # Model 1: Logistic Regression (One-vs-Rest)
    print("\n1. Logistic Regression (One-vs-Rest)")
    lr = LogisticRegression(
        multi_class='ovr',
        max_iter=2000,
        random_state=42,
        class_weight='balanced'
    )
    lr.fit(X_train_scaled, y_train)
    lr_pred = lr.predict(X_test_scaled)
    lr_acc = accuracy_score(y_test, lr_pred)
    print(f"   Test Accuracy: {lr_acc:.3f} ({int(lr_acc * len(y_test))}/{len(y_test)} correct)")
    models['Logistic Regression'] = (lr, lr_pred, lr_acc)

    # Model 2: Multinomial Logistic Regression
    print("\n2. Logistic Regression (Multinomial)")
    lr_multi = LogisticRegression(
        multi_class='multinomial',
        solver='lbfgs',
        max_iter=2000,
        random_state=42,
        class_weight='balanced'
    )
    lr_multi.fit(X_train_scaled, y_train)
    lr_multi_pred = lr_multi.predict(X_test_scaled)
    lr_multi_acc = accuracy_score(y_test, lr_multi_pred)
    print(f"   Test Accuracy: {lr_multi_acc:.3f} ({int(lr_multi_acc * len(y_test))}/{len(y_test)} correct)")
    models['Logistic Regression (Multinomial)'] = (lr_multi, lr_multi_pred, lr_multi_acc)

    # Model 3: Random Forest
    print("\n3. Random Forest")
    rf = RandomForestClassifier(
        n_estimators=200,
        max_depth=10,
        min_samples_split=5,
        random_state=42,
        class_weight='balanced',
        n_jobs=-1
    )
    rf.fit(X_train_scaled, y_train)
    rf_pred = rf.predict(X_test_scaled)
    rf_acc = accuracy_score(y_test, rf_pred)
    print(f"   Test Accuracy: {rf_acc:.3f} ({int(rf_acc * len(y_test))}/{len(y_test)} correct)")
    models['Random Forest'] = (rf, rf_pred, rf_acc)


    # --- Cross-validation ---
    print("\n" + "="*80)
    print("CROSS-VALIDATION SCORES (5-fold Stratified)")
    print("="*80)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    lr_cv_scores = cross_val_score(lr, X_train_scaled, y_train, cv=cv, scoring='accuracy')
    print(f"\nLogistic Regression (OvR):          {lr_cv_scores.mean():.3f} (+/- {lr_cv_scores.std():.3f})")

    lr_multi_cv_scores = cross_val_score(lr_multi, X_train_scaled, y_train, cv=cv, scoring='accuracy')
    print(f"Logistic Regression (Multinomial):  {lr_multi_cv_scores.mean():.3f} (+/- {lr_multi_cv_scores.std():.3f})")

    rf_cv_scores = cross_val_score(rf, X_train_scaled, y_train, cv=cv, scoring='accuracy')
    print(f"Random Forest:                       {rf_cv_scores.mean():.3f} (+/- {rf_cv_scores.std():.3f})")


    # --- Select best model ---
    best_model_name = max(models.items(), key=lambda x: x[1][2])[0]
    best_model, best_pred, best_acc = models[best_model_name]

    print(f"\n{'='*80}")
    print(f"BEST MODEL: {best_model_name}")
    print(f"Test Accuracy: {best_acc:.3f}")
    print("="*80)


    # --- Detailed classification report ---
    print("\nDetailed Classification Report:")
    print(classification_report(y_test, best_pred, target_names=class_names, zero_division=0))


    # --- Confusion Matrix ---
    cm = confusion_matrix(y_test, best_pred)
    print("\nConfusion Matrix:")
    print("(Rows = True, Columns = Predicted)")
    header = " " * 20 + " ".join([f"{name[:8]:>10}" for name in class_names])
    print(header)
    print("-" * len(header))
    for i, name in enumerate(class_names):
        row_label = f"{name[:18]:20}"
        row_values = " ".join([f"{cm[i, j]:>10}" for j in range(len(class_names))])
        print(row_label + row_values)


    # --- Feature importance (for Random Forest) ---
    if best_model_name == 'Random Forest' and hasattr(best_model, 'feature_importances_'):
        print(f"\n{'='*80}")
        print("Top 15 Most Important Markers (Random Forest):")
        print("="*80)
        importance_df = pd.DataFrame({
            'marker': [col.replace('_log', '') for col in X_imputed.columns],
            'importance': best_model.feature_importances_
        }).sort_values('importance', ascending=False)

        for idx, (_, row) in enumerate(importance_df.head(15).iterrows(), 1):
            print(f"  {idx:2}. {row['marker']:25} {row['importance']:.4f}")


    # --- Save results ---
    output_dir = Path(output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Predictions with probabilities
    results_df = pd.DataFrame({
        'sample_id': cancer_df_filtered.loc[X_test.index, 'Study_Subject_ID'],
        'true_type': le.inverse_transform(y_test),
        'predicted_type': le.inverse_transform(best_pred),
        'correct': y_test == best_pred
    })

    # Add prediction probabilities
    if hasattr(best_model, 'predict_proba'):
        proba = best_model.predict_proba(X_test_scaled)
        for i, class_name in enumerate(class_names):
            results_df[f'prob_{class_name}'] = proba[:, i]

    output_file = output_dir / "cancer_type_predictions.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n✅ Predictions saved to: {output_file}")


    # --- Save model summary ---
    summary = {
        'Model': best_model_name,
        'Test_Accuracy': best_acc,
        'CV_Mean': max(lr_cv_scores.mean(), lr_multi_cv_scores.mean(), rf_cv_scores.mean()),
        'Num_Classes': len(class_names),
        'Classes': ', '.join(class_names),
        'Num_Features': X_imputed.shape[1],
        'Train_Samples': len(X_train),
        'Test_Samples': len(X_test),
        'Total_Samples': len(cancer_df_filtered)
    }

    summary_df = pd.DataFrame([summary])
    summary_file = output_dir / "cancer_type_model_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"✅ Model summary saved to: {summary_file}")


    # --- Optional: Plot confusion matrix ---
    if save_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import seaborn as sns  # type: ignore

            plt.figure(figsize=(10, 8))
            sns.heatmap(
                cm,
                annot=True,
                fmt='d',
                cmap='Blues',
                xticklabels=class_names,
                yticklabels=class_names,
                cbar_kws={'label': 'Number of Samples'}
            )
            plt.xlabel('Predicted Cancer Type', fontsize=12)
            plt.ylabel('True Cancer Type', fontsize=12)
            plt.title(f'Confusion Matrix - {best_model_name}\nTest Accuracy: {best_acc:.1%}', fontsize=14)
            plt.tight_layout()

            plot_file = output_dir / "cancer_type_confusion_matrix.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"✅ Confusion matrix plot saved to: {plot_file}")

        except ImportError as e:
            print(f"⚠️  matplotlib/seaborn not available ({e}), skipping plot generation")


    return results_df, summary_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Train multi-class classifier to predict cancer type from methylation profiles"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to methyl_data_wide.csv file"
    )
    parser.add_argument(
        "--output",
        default="data",
        help="Output directory path (default: data)"
    )
    parser.add_argument(
        "--cancer-types",
        nargs="+",
        default=None,
        help="Specific cancer types to classify (e.g., --cancer-types Lung Pancreatic). Default: Lung Pancreatic"
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.25,
        help="Fraction of data for testing (default: 0.25)"
    )
    parser.add_argument(
        "--top-features",
        type=int,
        default=None,
        help="Use only top N features (optional feature selection)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate confusion matrix plot"
    )

    args = parser.parse_args()

    classify_cancer_types(
        input_path=args.input,
        output_path=args.output,
        cancer_types=args.cancer_types,
        test_size=args.test_size,
        top_n_features=args.top_features,
        save_plot=args.plot
    )
