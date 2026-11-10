================================================================================
SHARP: Sickle-cell Hemoglobin Analysis via Rules and Prediction
Standalone Offline Application for Windows
================================================================================

OVERVIEW:
SHARP is an explainable artificial intelligence (XAI) clinical decision support
framework that predicts hemoglobin genotypes (HbSS, HbSC, HbSβ⁰, HbSβ⁺, HbAS, and
HbAA) from quantitative hemoglobin electrophoresis and optional CBC parameters.
This standalone executable (SHARP.exe) runs 100% offline with zero external runtime
dependencies:
  - No R or Python installation required.
  - No internet connection required.
  - Self-contained Tier-1 Clinical Rules + Tier-2 5-Class Transfusion-Aware XGBoost
    Ensemble (740 decision trees) + exact multi-class TreeSHAP explainability engine
    embedded directly inside the binary.

--------------------------------------------------------------------------------
HOW TO RUN:
--------------------------------------------------------------------------------
1. Interactive Graphical User Interface (GUI):
   - Simply double-click "SHARP.exe".
   - A native desktop control window will launch along with the clean, modern
     interactive SHARP user interface.
   - You can choose between:
       * "Single Patient Prediction": Enter electrophoresis percentages and
         transfusion status to get immediate genotype predictions, confidence scores,
         clinical rationales, and an interactive multi-class TreeSHAP waterfall
         plot explaining the exact feature contributions. Download an HTML report.
       * "Batch Spreadsheet Analysis": Upload any CSV or Excel (.xlsx) file, or
         click "Load Demo Batch", to classify an entire cohort and download an
         exhaustive 75-column results spreadsheet including constrained 5-class
         probabilities, expert-review flags, grouped clinical SHAP, and all 30
         individual feature SHAP values.

2. Direct Batch File Processing (GUI Dialog):
   - Click "Quick Direct Batch Analysis" in the SHARP control window.
   - Select your input spreadsheet (CSV or Excel) and choose where to save the
     resulting CSV file.

3. Command-Line / Automated Batch Processing:
   - Run SHARP headlessly from PowerShell or Command Prompt:
       SHARP.exe --batch input_patients.csv output_results_with_shap.csv
   - To generate a blank template CSV from the command line:
       SHARP.exe --template template.csv

--------------------------------------------------------------------------------
SPREADSHEET FORMAT SPECIFICATIONS:
--------------------------------------------------------------------------------
Required Electrophoresis Columns (case-insensitive headers):
  1. hgb_s_pct      - Hemoglobin S percentage (%)
  2. hgb_a_pct      - Hemoglobin A percentage (%)
  3. hgb_c_pct      - Hemoglobin C percentage (%) (default: 0 if absent)
  4. hgb_f_pct      - Hemoglobin F percentage (%)
  5. hgb_a2_pct     - Hemoglobin A2 percentage (%)
  6. is_transfused  - Recent transfusion within 3 months (0/1 or No/Yes)

Optional CBC & Quality Columns:
  - patient_id      - Patient identifier or sample code
  - total_hgb       - Total Hemoglobin (g/dL)
  - mcv             - Mean Corpuscular Volume (fL)
  - mchc            - Mean Corpuscular Hemoglobin Concentration (g/dL)
  - rdw             - Red Cell Distribution Width (%)
  - hematocrit      - Hematocrit (%)
  - rbc_count       - Red Blood Cell Count (M/uL)
  - clarity_score   - Profile clarity score (1-10)

Output Spreadsheet Columns (75 Total Columns):
  All original columns plus:
  - predicted_genotype (SS, SC, S_beta_0, S_beta_plus, AS, AA)
  - predicted_genotype_label
  - prediction_confidence & prediction_confidence_pct
  - prediction_source (Clinical Rules vs. Machine Learning)
  - clinical_rule_applied
  - needs_expert_review (TRUE / FALSE)
  - clinical_rationale
  - prob_SS, prob_SC, prob_S_beta_0, prob_S_beta_plus, prob_AS
  - shap_target_genotype, shap_baseline_log_odds, shap_prediction_log_odds
  - Grouped SHAP columns: shap_hgb_s_group, shap_hgb_a_group, shap_hgb_c_group,
    shap_hgb_f_group, shap_hgb_a2_group, shap_transfusion_group, shap_mcv_group,
    shap_total_hgb, shap_mchc, shap_rdw, shap_rbc_count, shap_trait_indicator
  - Top positive and negative SHAP feature strings and clinical summary
  - All 30 individual baked model feature SHAP columns (shap_feature_*)

--------------------------------------------------------------------------------
CLINICAL NOTICE:
--------------------------------------------------------------------------------
SHARP is designed for research, educational, and clinical decision support purposes.
Cases with prediction confidence <= 90% or ambiguous features are automatically
flagged as "needs_expert_review = TRUE" for specialist review.
================================================================================
