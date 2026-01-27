# RainStream is a tool to predict tissue of origin of cancer using methylation signature

## Running indivdual python scripts outside of Nextflow pipeline must be done from project root

## Data Processing Pipeline

The pipeline transforms raw methylation data through four sequential stages:

### 1. Data Merge (`methyl_data_processing.py`)
Joins methylation measurements with clinical metadata.
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
