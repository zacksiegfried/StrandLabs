#!/usr/bin/env python3
"""
methyl_cancer_detection.py
--------------------------

Stage 1: Binary cancer / non-cancer classifier.

Input:
  - methyl_data_wide.csv  : wide-format methylation data
        Columns: Study_Subject_ID, CANCER_TYPE, cancer_yn, stage, AGE, SEX_D, [marker]_log
  - feature_importance_scores.csv (optional) : pre-ranked markers from methyl_feature_selection.py
        Used to subset to top N markers before training.

Output:
  - cancer_detection_predictions.csv  : per-sample predicted probability and classification
  - cancer_detection_model_summary.csv : AUC, sensitivity, specificity, accuracy at optimal threshold
  - cancer_detection_roc_curve.png (optional)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, roc_curve, accuracy_score,
    classification_report, confusion_matrix
)
import warnings
warnings.filterwarnings('ignore')


def _optimal_threshold(fpr, tpr, thresholds):
    """Youden's J: maximise sensitivity + specificity - 1."""
    j = tpr - fpr
    idx = np.argmax(j)
    return thresholds[idx], tpr[idx], 1 - fpr[idx]


def detect_cancer(
    input_path: str,
    output_path: str,
    features_path: str = None,
    top_n_features: int = None,
    test_size: float = 0.25,
    plot: bool = False,
) -> None:

    print("=" * 80)
    print("CANCER DETECTION — STAGE 1 BINARY CLASSIFIER")
    print("=" * 80)

    df = pd.read_csv(input_path)
    marker_cols = [col for col in df.columns if col.endswith('_log')]

    # --- Optional: subset to top-ranked markers ---
    if features_path and top_n_features:
        scores = pd.read_csv(features_path)
        top_markers = scores.head(top_n_features)['marker_full'].tolist()
        marker_cols = [m for m in top_markers if m in marker_cols]
        print(f"Using top {len(marker_cols)} markers from feature importance scores")
    else:
        print(f"Using all {len(marker_cols)} markers")

    X = df[marker_cols]
    y = df['cancer_yn']

    print(f"Loaded {len(df)} samples | {y.sum()} cancer | {(y == 0).sum()} non-cancer")

    # --- Impute missing values ---
    imputer = SimpleImputer(strategy='median')
    X_imputed = pd.DataFrame(imputer.fit_transform(X), columns=marker_cols, index=X.index)

    missing = X.isnull().sum().sum()
    if missing > 0:
        print(f"Imputed {missing} missing values (median strategy)")

    # --- Train / test split ---
    X_train, X_test, y_train, y_test = train_test_split(
        X_imputed, y, test_size=test_size, random_state=42, stratify=y
    )
    print(f"\nTrain: {len(X_train)} samples | Test: {len(X_test)} samples")

    # --- Scale ---
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    # --- Train models ---
    print("\n" + "=" * 80)
    print("TRAINING")
    print("=" * 80)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    models = {}

    lr = LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced')
    lr.fit(X_train_sc, y_train)
    lr_cv = cross_val_score(lr, X_train_sc, y_train, cv=cv, scoring='roc_auc')
    print(f"\nLogistic Regression  CV AUC: {lr_cv.mean():.3f} (+/- {lr_cv.std():.3f})")
    models['Logistic Regression'] = (lr, lr_cv.mean())

    rf = RandomForestClassifier(
        n_estimators=200, max_depth=10, min_samples_split=5,
        random_state=42, class_weight='balanced', n_jobs=-1
    )
    rf.fit(X_train_sc, y_train)
    rf_cv = cross_val_score(rf, X_train_sc, y_train, cv=cv, scoring='roc_auc')
    print(f"Random Forest        CV AUC: {rf_cv.mean():.3f} (+/- {rf_cv.std():.3f})")
    models['Random Forest'] = (rf, rf_cv.mean())

    # --- Select best model by CV AUC ---
    best_name = max(models, key=lambda k: models[k][1])
    best_model = models[best_name][0]
    print(f"\nBest model: {best_name}")

    # --- Evaluate on test set ---
    y_prob = best_model.predict_proba(X_test_sc)[:, 1]
    test_auc = roc_auc_score(y_test, y_prob)

    fpr, tpr, thresholds = roc_curve(y_test, y_prob)
    opt_thresh, sensitivity, specificity = _optimal_threshold(fpr, tpr, thresholds)

    y_pred = (y_prob >= opt_thresh).astype(int)
    accuracy = accuracy_score(y_test, y_pred)

    print(f"\n{'=' * 80}")
    print("TEST SET PERFORMANCE (Youden optimal threshold)")
    print("=" * 80)
    print(f"  AUC:         {test_auc:.3f}")
    print(f"  Sensitivity: {sensitivity:.3f}")
    print(f"  Specificity: {specificity:.3f}")
    print(f"  Accuracy:    {accuracy:.3f}")
    print(f"  Threshold:   {opt_thresh:.3f}")
    print(f"\n{classification_report(y_test, y_pred, target_names=['Non-Cancer', 'Cancer'], zero_division=0)}")

    # --- Save outputs ---
    out_dir = Path(output_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    predictions = pd.DataFrame({
        'Study_Subject_ID': df.loc[X_test.index, 'Study_Subject_ID'],
        'true_label': y_test.values,
        'cancer_probability': y_prob.round(4),
        'predicted_label': y_pred,
        'correct': (y_test.values == y_pred),
    })
    pred_file = out_dir / "cancer_detection_predictions.csv"
    predictions.to_csv(pred_file, index=False)
    print(f"✅ Predictions saved → {pred_file}")

    summary = pd.DataFrame([{
        'best_model': best_name,
        'n_markers': len(marker_cols),
        'train_samples': len(X_train),
        'test_samples': len(X_test),
        'cv_auc': round(max(lr_cv.mean(), rf_cv.mean()), 4),
        'test_auc': round(test_auc, 4),
        'sensitivity': round(sensitivity, 4),
        'specificity': round(specificity, 4),
        'accuracy': round(accuracy, 4),
        'optimal_threshold': round(opt_thresh, 4),
    }])
    summary_file = out_dir / "cancer_detection_model_summary.csv"
    summary.to_csv(summary_file, index=False)
    print(f"✅ Model summary saved → {summary_file}")

    # --- Optional ROC curve ---
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(fpr, tpr, color="#2196F3", linewidth=2, label=f"AUC = {test_auc:.3f}")
        ax.plot([0, 1], [0, 1], color="grey", linewidth=1, linestyle="--")
        ax.scatter(1 - specificity, sensitivity, color="#E53935", zorder=5,
                   label=f"Optimal (Se={sensitivity:.2f}, Sp={specificity:.2f})")
        ax.set_xlabel("1 - Specificity", fontsize=10)
        ax.set_ylabel("Sensitivity", fontsize=10)
        ax.set_title(f"ROC Curve — {best_name}", fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()

        plot_file = out_dir / "cancer_detection_roc_curve.png"
        fig.savefig(plot_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"✅ ROC curve saved → {plot_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Stage 1 binary classifier: cancer vs non-cancer."
    )
    parser.add_argument("--input",    required=True, help="Path to methyl_data_wide.csv")
    parser.add_argument("--output",   default="output", help="Output directory")
    parser.add_argument("--features", default=None, help="Path to feature_importance_scores.csv")
    parser.add_argument("--top-features", type=int, default=None,
                        help="Use top N features from importance scores")
    parser.add_argument("--test-size", type=float, default=0.25,
                        help="Test set fraction (default: 0.25)")
    parser.add_argument("--plot", action="store_true", help="Generate ROC curve plot")
    args = parser.parse_args()

    detect_cancer(
        input_path=args.input,
        output_path=args.output,
        features_path=args.features,
        top_n_features=args.top_features,
        test_size=args.test_size,
        plot=args.plot,
    )
