# RainStream is a tool to predict tissue of origin of cancer using methylation signature

## Running the Pipeline

```
source ../.venv/bin/activate
nextflow run main.nf

bash cleanup.sh
```

### Running individual python scripts outside of Nextflow pipeline must be done from project root

## Data Processing Pipeline

The pipeline transforms raw methylation data through four sequential stages:

### 1. Data Merge (`methyl_data_processing.py`)
Joins methylation measurements with clinical metadata. Filters B3GALT6 out.
- **Input**: `clinical.csv` (patient metadata), `methylation.csv` (long-format, 3 reps/patient)
- **Output**: `methyl_data_trim.csv`

### 2. Replicate Handling (`methyl_rep_handling.py`)
Condenses technical replicates into summary statistics per marker per patient.
- **Output**: `methyl_data_trim_condensed.csv` with `log_mean`, `log_sd`, `log_cv`

### 3. Wide Formatting (`methyl_wide_formatting.py`)
Pivots data so each methylation marker becomes a column.
- **Output**: `methyl_data_wide.csv` (one row per patient)

### 4. Feature Selection (`methyl_feature_selection.py`)
Ranks markers by importance using ANOVA F-test, Random Forest, and Logistic Regression.
- **Output**: `feature_importance_scores.csv`

## Analysis Output

### Marker Precision Profiles (`methyl_precision_profile.py`)
Runs in parallel from Step 1. Models the precision of each marker across the quantitative measuring range. Summary statistics computed on replicates of same sample.
- **Output**: `marker_precision_profiles/*.png`

## Predicition Pipeline

### Stage 1 — Cancer Detection (`methyl_cancer_detection.py`)
Binary classifier (cancer vs non-cancer). Uses top-ranked markers from feature selection. Reports AUC, sensitivity, and specificity at the Youden-optimal threshold.
- **Output**: `cancer_detection_predictions.csv`, `cancer_detection_model_summary.csv`, `cancer_detection_roc_curve.png`

### Stage 2 — Tissue of Origin (`cancer_type_classification.py`)
Multi-class classifier trained on cancer samples only to predict `CANCER_TYPE`. Runs LR (OvR), LR (multinomial), and Random Forest with 5-fold CV.
- **Output**: `cancer_type_predictions.csv`, `cancer_type_model_summary.csv`, `cancer_type_confusion_matrix.png`
