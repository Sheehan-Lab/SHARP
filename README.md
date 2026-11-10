# SHARP: Sickle Cell Hemoglobinopathy Automated Resolution Program

SHARP is an explainable machine learning clinical decision support framework for sickle cell disease genotyping. It classifies hemoglobin variants (HbSS, HbSC, HbS/beta-plus, HbS/beta-zero, and sickle trait/other variants) using quantitative hemoglobin electrophoresis fractions and optional Complete Blood Count (CBC) parameters.

SHARP uses a two-tier architecture:
1. **Tier-1 Clinical Decision Rules**: Applies established hematological boundary conditions to resolve definitive classic phenotypes and detect recent transfusions.
2. **Tier-2 Transfusion-Aware XGBoost Ensemble**: Resolves complex, borderline, and post-transfusion profiles using a 5-class gradient boosted model integrated with multi-class TreeSHAP explainability.

## Repository Structure

- `1_Model_Training/`
  Scripts for data preparation, rule definitions, and 5-class XGBoost model training (`enhanced_5class_allcases_pipeline.R`, `data_preparation_v2.R`, `utils/clinical_rules.R`), along with trained workflow artifacts (`tidymodels_genotype_workflow_with_transfusion.rds`).
- `2_Shiny_App/`
  Interactive R Shiny web application (`app.R`) for single-patient assessment and batch spreadsheet analysis with interactive SHAP waterfall plots.
- `3_Standalone_Executable/`
  Self-contained Windows executable (`SHARP.exe`), batch template (`SHARP_Batch_Template.csv`), and instructions for offline execution without R or external dependencies.

## Quick Start

### 1. Standalone Windows Executable
To run SHARP locally on Windows without installing R or Python:
1. Navigate to `3_Standalone_Executable/` and open `SHARP.exe`.
2. Use the Single Patient tab for individual case review, or the Batch tab to process a cohort spreadsheet.
3. For command-line batch runs:
   ```cmd
   SHARP.exe --batch input_cases.csv output_results.csv
   ```

### 2. R Shiny Application
To run the interactive Shiny application locally:
```R
install.packages(c("shiny", "bslib", "DT", "xgboost", "tidymodels", "ggplot2", "plotly", "readr", "readxl"))
shiny::runApp("2_Shiny_App")
```

### 3. Model Training Pipeline
Training code and methodology documentation are located in `1_Model_Training/`. See `S.H.A.R.P_Clinical_Constraints_Methods.md` for details on rule thresholds and constraint logic.

## Input Features

### Required Electrophoresis Fractions
- `hgb_s_pct`: Hemoglobin S (%)
- `hgb_a_pct`: Hemoglobin A (%)
- `hgb_c_pct`: Hemoglobin C (%)
- `hgb_f_pct`: Hemoglobin F (%)
- `hgb_a2_pct`: Hemoglobin A2 (%)
- `is_transfused`: Recent blood transfusion within past 90-120 days (1 = Yes, 0 = No)

### Optional CBC Parameters
- `total_hgb`: Total Hemoglobin (g/dL)
- `mcv`: Mean Corpuscular Volume (fL)
- `mchc`: Mean Corpuscular Hemoglobin Concentration (g/dL)
- `rdw`: Red Cell Distribution Width (%)
- `hematocrit`: Hematocrit (%)
- `rbc_count`: Red Blood Cell Count (x10^6/uL)
- `clarity_score`: Hemoglobin fraction sum sanity score (0-10)

A sample batch file is provided in `3_Standalone_Executable/SHARP_Batch_Template.csv`.

## Output and Interpretation

- **Predicted Genotype**: Primary classification across major sickle variants (HbSS, HbSC, HbS/beta-0, HbS/beta-+, Trait/Other).
- **Class Probabilities**: Calibrated probability distribution across all five target genotype categories.
- **Classification Method**: Identifies whether the diagnosis was determined via Tier-1 clinical heuristics or Tier-2 machine learning inference.
- **TreeSHAP Attributions**: Quantitative impact of each electrophoresis fraction and CBC parameter on the predicted class.

## Citation and License

This software is provided under the terms of the repository license. If you use SHARP in your research, please cite the corresponding publication from the Sheehan Lab.
