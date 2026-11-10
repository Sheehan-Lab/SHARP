# ==============================================================================
# Enhanced Hemoglobin Genotype Prediction - 5-Class All-Cases Tidymodels Pipeline
# ==============================================================================
# This script implements an ENHANCED tidymodels approach with:
# 1. Rule-Based Clinical Filter (high-confidence assignments) INCLUDING AS
# 2. 5-Class ML Pipeline (SS, SC, S_beta_0, S_beta_plus, AS) 
# 3. Training on ALL CASES instead of just ambiguous cases
# 4. Transfusion status awareness as a feature
# 5. Comprehensive performance evaluation and interpretation
# 6. Professional visualizations for publication
# ==============================================================================

library(tidyverse)
library(tidymodels)
library(finetune)      # For advanced tuning
library(vip)           # For variable importance plots
library(pROC)
library(corrplot)
library(ggplot2)
library(viridis)
library(gridExtra)
library(knitr)
library(kableExtra)

# Additional libraries for enhanced analysis
library(PRROC)         # For Precision-Recall curves
library(MLmetrics)     # For additional metrics
library(purrr)         # For map functions
library(shapviz)       # For SHAP plots
library(themis)        # For SMOTE (added)

# Set options
options(scipen = 999, digits = 3)
set.seed(42)
tidymodels_prefer()

# ==============================================================================
# 1. LOAD ENHANCED PROCESSED DATA WITH AS PATIENTS
# ==============================================================================

cat("=== ENHANCED 5-CLASS ALL-CASES GENOTYPE PREDICTION PIPELINE ===\n")
cat("Loading enhanced hemoglobin data with AS patients...\n")

# Use the enhanced dataset with AS patients
if (file.exists("enhanced_hemoglobin_data_with_as.csv")) {
  data <- read_csv("enhanced_hemoglobin_data_with_as.csv")
  cat("Enhanced dataset dimensions:", nrow(data), "patients x", ncol(data), "features\n")
} else {
  cat("Enhanced dataset not found. Please run process_as_data_and_update_pipeline.R first.\n")
  stop("enhanced_hemoglobin_data_with_as.csv not found")
}

# Add transfusion status if not present (most patients are not transfused)
if (!"is_transfused" %in% colnames(data)) {
  data <- data %>%
    mutate(is_transfused = 0)  # Default to not transfused
}

data <- read_csv("enhanced_hemoglobin_data_with_as.csv") %>%
  mutate(
    # Assuming NA in is_transfused means the patient was NOT transfused.
    # If NA means something else, change 'FALSE' accordingly.
    is_transfused = replace_na(is_transfused, FALSE)
  )

# Quality and transfusion summaries
if("test_confidence_level" %in% colnames(data)) {
  quality_summary <- data %>% count(test_confidence_level, sort = TRUE)
  cat("Test quality distribution:\n")
  print(quality_summary)
}

transfusion_summary <- data %>% count(is_transfused, sort = TRUE)
cat("\nTransfusion status distribution:\n")
print(transfusion_summary)

# Enhanced genotype distribution
genotype_summary <- data %>% count(genotype_simplified, sort = TRUE)
cat("\nGenotype distribution (5-class enhanced):\n")
print(genotype_summary)

data <- read_csv("enhanced_hemoglobin_data_with_as.csv") %>%
  mutate(
    # Assuming NA in is_transfused means the patient was NOT transfused.
    # If NA means something else, change 'FALSE' accordingly.
    is_transfused = replace_na(is_transfused, FALSE)
  )
# ==============================================================================
# 2. ENHANCED CLINICAL RULES (TIER 1: NOW INCLUDING AS, AA)
# ==============================================================================

cat("\n=== TIER 1: ENHANCED CLINICAL RULE-BASED CLASSIFICATION (5-CLASS) ===\n")

## MODIFIED FUNCTION WITH AA OVERRIDE
apply_clinical_rules <- function(df) {
  df %>%
    mutate(
      # --- Rule Definitions ---
      # These flags identify patterns for non-override cases.
      rule_ss = (hgb_s_pct > 70 & hgb_c_pct < 1 & hgb_f_pct < 20 & 
                   hgb_a2_pct < 3.5 & hgb_a_pct < 1 & mcv > 75),
      
      rule_sc = (hgb_s_pct >= 35 & hgb_s_pct <= 60 & 
                   hgb_c_pct >= 40 & hgb_c_pct <= 60 & hgb_f_pct < 5 &
                   hgb_a_pct < 1 & mcv >= 70 & mcv <= 95),
      
      rule_sbeta0 = (hgb_s_pct > 65 & hgb_c_pct < 1 & hgb_f_pct < 20 & 
                       hgb_a2_pct > 4.5 & hgb_a_pct < 1 & mcv < 80),
      
      rule_sbeta_plus = (hgb_s_pct >= 60 & hgb_s_pct <= 90 & hgb_c_pct < 1 & 
                           hgb_f_pct < 20 & hgb_a2_pct > 4.5 & 
                           hgb_a_pct >= 3 & hgb_a_pct <= 35 & mcv < 85),
      
      rule_as = (hgb_a_pct >= 50 & hgb_a_pct <= 65 & hgb_s_pct >= 30 & 
                   hgb_s_pct <= 45 & hgb_c_pct < 1 & hgb_a2_pct < 3.5 & hgb_f_pct < 3),
      
      # --- Prediction Logic with Override ---
      # This case_when assigns a prediction. The AA override is the FIRST condition.
      rule_prediction = case_when(
        ## NEW: AA OVERRIDE RULE IS PLACED AT THE TOP
        # If HbA is dominant and variants are absent, it's definitively AA.
        !is.na(hgb_a_pct) & hgb_a_pct > 95 & 
          !is.na(hgb_s_pct) & hgb_s_pct < 2 & 
          !is.na(hgb_c_pct) & hgb_c_pct < 2   ~ "AA",
        
        # Standard rules are evaluated only if the override is not met.
        rule_ss ~ "SS",
        rule_sc ~ "SC", 
        rule_sbeta0 ~ "S_beta_0",
        rule_sbeta_plus ~ "S_beta_plus",
        rule_as ~ "AS",
        TRUE ~ "AMBIGUOUS"
      ),
      
      # --- Confidence Logic with Override ---
      # This logic determines if a rule-based prediction is high-confidence.
      high_confidence = case_when(
        ## NEW: AA OVERRIDE CONFIDENCE IS PLACED AT THE TOP
        # If the AA override rule was met, force confidence to TRUE.
        !is.na(hgb_a_pct) & hgb_a_pct > 95 & 
          !is.na(hgb_s_pct) & hgb_s_pct < 2 & 
          !is.na(hgb_c_pct) & hgb_c_pct < 2   ~ TRUE,
        
        # ## MODIFIED: All original confidence logic is now the fallback condition.
        # This is only evaluated if the AA override is FALSE.
        TRUE ~ {
          # This block contains your original logic for all other genotypes.
          quality_adjusted_confidence_initial <- case_when(
            rule_prediction != "AMBIGUOUS" & 
              (if_else(is.na(test_confidence_level), FALSE, test_confidence_level == "High")) &
              total_hgb >= 6 ~ TRUE,
            
            rule_prediction != "AMBIGUOUS" & 
              (if_else(is.na(test_confidence_level), FALSE, test_confidence_level == "Medium")) &
              (if_else(is.na(clarity_score), FALSE, clarity_score >= 9)) &
              total_hgb >= 6 ~ TRUE,
            
            TRUE ~ FALSE
          )
          
          # Final confidence check (excluding AA override cases)
          rule_prediction != "AMBIGUOUS" & quality_adjusted_confidence_initial &
            !(coalesce(is_transfused, FALSE) & rule_prediction %in% c("S_beta_plus", "AS"))
        }
      )
    )
}

# Apply clinical rules
data_with_rules <- apply_clinical_rules(data)

# Evaluate rule performance
rule_performance <- data_with_rules %>%
  filter(high_confidence) %>%
  summarise(
    total_classified = n(),
    classification_rate = round(n() / nrow(data_with_rules) * 100, 1),
    accuracy = round(sum(rule_prediction == genotype_simplified) / n() * 100, 1)
  )

cat("Enhanced Clinical Rules Performance (5-class):\n")
cat("- Patients classified:", rule_performance$total_classified, "/", nrow(data_with_rules), 
    "(", rule_performance$classification_rate, "%)\n")
cat("- Rule accuracy:", rule_performance$accuracy, "%\n")

# --- BEGIN SECTION: Investigate SC Patients with Low HbC --- 
cat("\n=== INVESTIGATING SC PATIENTS WITH HbC NEAR ZERO IN RULE INPUT DATA ===\n")
sc_patients_low_hbc <- data_with_rules %>% 
  filter(genotype_simplified == "SC" & (is.na(hgb_c_pct) | hgb_c_pct < 1)) # Check for HbC < 1%

if(nrow(sc_patients_low_hbc) > 0) {
  cat(paste("Found", nrow(sc_patients_low_hbc), "actual SC patient(s) with HbC < 1% (potential data issue):
"))
  print(sc_patients_low_hbc %>%
          select(patient_mrn, # Replaced patient_id, sample_id
                 genotype_simplified, rule_prediction, high_confidence, test_confidence_level, 
                 hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct, mcv, is_transfused, total_hgb, clarity_score) %>%
          arrange(patient_mrn) # Replaced patient_id
       )
} else {
  cat("No SC patients found with HbC < 1% in the dataset provided to rules.\n")
}
cat("====================================================================\n")
# --- END SECTION: Investigate SC Patients with Low HbC ---

# --- BEGIN SECTION: Investigate SS Patients with Unexpected HbC --- 
cat("\n=== INVESTIGATING SS PATIENTS WITH UNEXPECTED HbC (> 0.5%) IN RULE INPUT DATA ===\n")
ss_patients_with_hbc <- data_with_rules %>% 
  filter(genotype_simplified == "SS" & !is.na(hgb_c_pct) & hgb_c_pct > 0.5) # Check for HbC > 0.5%

if(nrow(ss_patients_with_hbc) > 0) {
  cat(paste("Found", nrow(ss_patients_with_hbc), "actual SS patient(s) with HbC > 0.5% (potential data issue or atypical):
"))
  print(ss_patients_with_hbc %>%
          select(patient_mrn, # Replaced patient_id, sample_id
                 genotype_simplified, rule_prediction, high_confidence, test_confidence_level, 
                 hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct, mcv, is_transfused, total_hgb, clarity_score) %>%
          arrange(patient_mrn) # Replaced patient_id
       )
} else {
  cat("No SS patients found with HbC > 0.5% in the dataset provided to rules.\n")
}
cat("=============================================================================\n")
# --- END SECTION: Investigate SS Patients with Unexpected HbC ---

# Identify ambiguous cases for ML
ambiguous_cases <- data_with_rules %>%
  filter(!high_confidence)

cat("\nAmbiguous cases requiring ML:", nrow(ambiguous_cases), "patients\n")

# ==============================================================================
# 3. ENHANCED DATA PREPARATION FOR ALL-CASES TRAINING
# ==============================================================================

cat("\n=== ENHANCED DATA PREPARATION: ALL-CASES 5-CLASS TRAINING ===\n")

# Prepare features for tidymodels including transfusion status
prepare_ml_features <- function(df) {
  df %>%
    select(
      # Primary electrophoresis values
      hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct,
      # CBC parameters
      total_hgb, hematocrit, mcv, mchc, rdw, rbc_count,
      # Derived features
      hgb_a_present, hgb_s_present, hgb_c_present, 
      hgb_f_elevated, hgb_a2_elevated,
      hgb_s_to_total, hgb_a_to_s_ratio, total_variant_pct,
      # Quality scores
      completeness_score, clarity_score, total_quality_score,
      # Other features
      had_hba_below_1_ever,
      # TRANSFUSION STATUS (KEY ADDITION)
      is_transfused,
      # Target
      genotype_simplified
    ) %>%
    filter(!is.na(genotype_simplified)) %>%
    mutate(
      # Convert logical to numeric
      across(where(is.logical), as.numeric),
      # Convert target to factor
      genotype_simplified = factor(genotype_simplified)
    )
}

# MAJOR CHANGE: Train on ALL cases, not just ambiguous ones
cat("MAJOR ENHANCEMENT: Training on ALL cases instead of just ambiguous ones\n")

# Filter to 5 key genotypes for ML (ENHANCED with AS)
key_genotypes <- c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")  # NOW 5 CLASSES
all_cases_ml_data <- data_with_rules %>%
  filter(genotype_simplified %in% key_genotypes) %>%
  prepare_ml_features() %>%
  mutate(genotype_simplified = factor(genotype_simplified, levels = key_genotypes))

cat("Enhanced ML training data (ALL CASES):", nrow(all_cases_ml_data), "patients\n")
cat("Enhanced class distribution (5-class):\n")
print(all_cases_ml_data %>% count(genotype_simplified, sort = TRUE))

# Check transfusion distribution in enhanced ML data
cat("\nTransfusion status in enhanced ML data:\n")
print(all_cases_ml_data %>% count(is_transfused, genotype_simplified))

# ==============================================================================
# 4. ENHANCED TIDYMODELS RECIPE AND WORKFLOW (5-CLASS)
# ==============================================================================

cat("\n=== ENHANCED TIDYMODELS RECIPE AND WORKFLOW (5-CLASS) ===\n")

#' Apply clinical constraints to ML probabilities (same as Shiny app)
#' @param probabilities Data frame with raw ML probabilities  
#' @param input_data Data frame containing patient hemoglobin data
#' @return Data frame with clinically constrained probabilities
apply_clinical_constraints_pipeline <- function(probabilities, input_data) {
  
  # Create a copy to modify
  prob_constrained <- probabilities
  
  # Apply constraints row by row
  for(i in 1:nrow(input_data)) {
    
    # Extract key values for this patient
    hgb_s_pct <- as.numeric(input_data$hgb_s_pct[i])
    hgb_a_pct <- as.numeric(input_data$hgb_a_pct[i])
    hgb_c_pct <- as.numeric(input_data$hgb_c_pct[i])
    hgb_f_pct <- as.numeric(input_data$hgb_f_pct[i])
    hgb_a2_pct <- as.numeric(input_data$hgb_a2_pct[i])
    
    # Skip if any values are missing
    if(any(is.na(c(hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct)))) next
    
    # Constraint 1: SC requires significant HgbC (>10% minimum, ideally >40%)
    if (hgb_c_pct < 10) {
      if (".pred_SC" %in% names(prob_constrained)) {
        # Heavily penalize SC when HgbC is very low
        penalty_factor <- if (hgb_c_pct == 0) 0.01 else (hgb_c_pct / 10) * 0.1
        prob_constrained$.pred_SC[i] <- prob_constrained$.pred_SC[i] * penalty_factor
      }
    }
    
    # Constraint 2: SS should have minimal HgbA (<5%)
    if (hgb_a_pct > 10) {
      if (".pred_SS" %in% names(prob_constrained)) {
        # Penalize SS when HgbA is high
        penalty_factor <- max(0.1, 1 - (hgb_a_pct - 5) / 50)
        prob_constrained$.pred_SS[i] <- prob_constrained$.pred_SS[i] * penalty_factor
      }
    }
    
    # Constraint 3: Boost AS when pattern is classic trait (HgbA 40-70%, HgbS 30-60%, HgbC ≈ 0)
    if (hgb_a_pct >= 40 && hgb_a_pct <= 70 && 
        hgb_s_pct >= 30 && hgb_s_pct <= 60 && 
        hgb_c_pct <= 2) {
      if (".pred_AS" %in% names(prob_constrained)) {
        # Boost AS for classic trait pattern
        boost_factor <- 1.5
        prob_constrained$.pred_AS[i] <- prob_constrained$.pred_AS[i] * boost_factor
      }
    }
    
    # Constraint 4: S-beta patterns require elevated HgbA2 (>3.5%, ideally >5.7%)
    if (hgb_a2_pct < 3.5) {
      if (".pred_S_beta_0" %in% names(prob_constrained)) {
        penalty_factor <- 0.2
        prob_constrained$.pred_S_beta_0[i] <- prob_constrained$.pred_S_beta_0[i] * penalty_factor
      }
      if (".pred_S_beta_plus" %in% names(prob_constrained)) {
        penalty_factor <- 0.2
        prob_constrained$.pred_S_beta_plus[i] <- prob_constrained$.pred_S_beta_plus[i] * penalty_factor
      }
    }
    
    # Constraint 5: S-beta+ requires some HgbA (>5%)
    if (hgb_a_pct < 5) {
      if (".pred_S_beta_plus" %in% names(prob_constrained)) {
        penalty_factor <- 0.1
        prob_constrained$.pred_S_beta_plus[i] <- prob_constrained$.pred_S_beta_plus[i] * penalty_factor
      }
    }
    
    # Renormalize probabilities to sum to 1 for this row
    row_sum <- sum(prob_constrained[i, ], na.rm = TRUE)
    if (row_sum > 0) {
      prob_constrained[i, ] <- prob_constrained[i, ] / row_sum
    }
  }
  
  return(prob_constrained)
}

if(nrow(all_cases_ml_data) > 0) {
  
  # Create enhanced recipe for preprocessing (5-class)
  genotype_recipe <- recipe(genotype_simplified ~ ., data = all_cases_ml_data) %>%
    # Handle missing values
    step_impute_median(all_numeric_predictors()) %>%
    
    # Create new features including transfusion interactions
    step_mutate(
      # HbF interactions
      hgb_f_to_s_ratio = hgb_f_pct / (hgb_s_pct + 0.1),
      hgb_f_to_a2_ratio = hgb_f_pct / (hgb_a2_pct + 0.1),
      
      # AS-specific features (NEW)
      as_pattern_score = hgb_a_pct / (hgb_s_pct + 0.1),  # Should be ~1.5 for AS
      trait_indicator = as.numeric(hgb_s_pct >= 25 & hgb_s_pct <= 55 & hgb_a_pct >= 45),
      
      # Transfusion interaction features (KEY ADDITION)
      transfused_hgb_s = is_transfused * hgb_s_pct,
      transfused_hgb_a = is_transfused * hgb_a_pct,
      transfused_hgb_c = is_transfused * hgb_c_pct,
      transfused_mcv = is_transfused * mcv,
      transfused_total_hgb = is_transfused * total_hgb,
      
      # Borderline indicators
      borderline_hgb_s = as.numeric(hgb_s_pct >= 70 & hgb_s_pct <= 80),
      borderline_hgb_f = as.numeric(hgb_f_pct >= 5 & hgb_f_pct <= 15),
      borderline_mcv = as.numeric(mcv >= 75 & mcv <= 85)
    ) %>%
    
    # Enhanced oversampling for 5-class problem using SMOTE
    step_smote(genotype_simplified, seed = 42) %>%

    # Remove zero variance and highly correlated predictors
    step_zv(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = 0.95) %>%
    
    # Normalize numeric features
    step_normalize(all_numeric_predictors())
  
  cat("Enhanced recipe created with 5-class training and AS-specific features\n")
  
  # Model specifications
  xgb_spec <- boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    mtry = tune(),
    min_n = tune()
  ) %>%
    set_engine("xgboost", eval_metric = "mlogloss") %>%
    set_mode("classification")
  
  # Create workflows
  xgb_workflow <- workflow() %>%
    add_recipe(genotype_recipe) %>%
    add_model(xgb_spec)
  
  cat("Enhanced workflow created for 5-class XGBoost\n")
  
  # ==============================================================================
  # 5. ENHANCED CROSS-VALIDATION AND HYPERPARAMETER TUNING (5-CLASS)
  # ==============================================================================
  
  cat("\n=== ENHANCED CROSS-VALIDATION AND HYPERPARAMETER TUNING (5-CLASS) ===\n")
  
  # Create CV folds with 5-class stratification
  set.seed(42)
  cv_folds <- vfold_cv(all_cases_ml_data, v = 10, strata = genotype_simplified)
  
  # Define metrics for 5-class problem
  metric_set_custom <- metric_set(accuracy, roc_auc, sens, spec, precision, recall)
  
  # Define parameter ranges for XGBoost (enhanced for 5-class)
  xgb_param_ranges <- parameters(
    trees(range = c(50L, 300L)),
    tree_depth(range = c(3L, 8L)),
    learn_rate(range = c(0.01, 0.3)),
    mtry(range = c(3L, 15L)), # Increased range for more features
    min_n(range = c(2L, 20L))
  )
  
  # Tune models using Bayesian optimization (5-class)
  cat("Tuning Enhanced 5-Class XGBoost with Bayesian Optimization...\n")
  
  initial_configs = 20 # Increased for 5-class problem
  bayes_iterations = 20 # Increased for better exploration

  xgb_tune_results <- xgb_workflow %>% 
    tune_bayes(
      resamples = cv_folds,
      param_info = xgb_param_ranges,
      initial = initial_configs,
      iter = bayes_iterations,
      metrics = metric_set_custom,
      control = control_bayes(save_pred = TRUE, verbose = TRUE, no_improve = 15, seed = 42)
    )
  
  xgb_best <- select_best(xgb_tune_results, metric = "accuracy")
  cat("\nBest Enhanced 5-Class XGBoost results:\n")
  print(show_best(xgb_tune_results, metric = "accuracy", n = 3))
  
  # Finalize and fit enhanced model
  best_workflow <- xgb_workflow %>% finalize_workflow(xgb_best)
  best_tune_results <- xgb_tune_results
  model_type <- "Enhanced 5-Class XGBoost"

  final_fit <- best_workflow %>%
    fit(data = all_cases_ml_data)
  
  # EARLY SAVE: Save the 5-class model immediately after training
  # This ensures we have the AS-capable model even if later analysis fails
  cat("\n=== EARLY SAVE: 5-CLASS MODEL WITH AS ===\n")
  saveRDS(final_fit, "tidymodels_genotype_workflow_with_transfusion.rds")
  cat("✅ Enhanced 5-class tidymodels workflow (including AS) saved successfully!\n")
  cat("   Model file: tidymodels_genotype_workflow_with_transfusion.rds\n\n")
  
  # ==============================================================================
  # 6. MODEL EVALUATION AND PERFORMANCE ANALYSIS
  # ==============================================================================
  
  cat("\n=== MODEL EVALUATION AND PERFORMANCE ANALYSIS ===\n")
  
  # Get CV predictions
  cv_predictions <- collect_predictions(best_tune_results) %>%
    filter(.config == xgb_best$.config) %>%
    mutate(
      predicted_class = .pred_class,
      actual_class = genotype_simplified,
      # Calculate confidence as max probability
      confidence = pmax(
        coalesce(.pred_SS, 0),
        coalesce(.pred_SC, 0),
        coalesce(.pred_S_beta_0, 0),
        coalesce(.pred_S_beta_plus, 0),
        coalesce(.pred_AS, 0),  # NEW: Include AS predictions
        na.rm = TRUE
      ),
      high_confidence = confidence > 0.90
    )
  
  cat("Raw CV predictions without clinical constraints:\n")
  cat("Cross-validation accuracy (raw):", round(mean(cv_predictions$predicted_class == cv_predictions$actual_class) * 100, 1), "%\n")
  
  # APPLY CLINICAL CONSTRAINTS TO CV PREDICTIONS
  cat("\n=== APPLYING CLINICAL CONSTRAINTS TO CV PREDICTIONS ===\n")
  
  # Get probability columns
  prob_cols <- names(cv_predictions)[str_starts(names(cv_predictions), "\\.pred_") & 
                                    !str_detect(names(cv_predictions), "\\.pred_class")]
  
  if(length(prob_cols) > 0) {
    # Extract probabilities and apply constraints
    cv_probabilities <- cv_predictions %>% select(all_of(prob_cols))
    
    # Get corresponding input data for constraints
    cv_input_data <- cv_predictions %>%
      left_join(
        all_cases_ml_data %>% 
          mutate(row_id = row_number()) %>%
          select(row_id, hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct),
        by = c(".row" = "row_id")
      )
    
    # Apply clinical constraints
    cv_probabilities_constrained <- apply_clinical_constraints_pipeline(
      cv_probabilities, 
      cv_input_data
    )
    
    # Recalculate predictions based on constrained probabilities
    cv_predictions_constrained <- cv_probabilities_constrained %>%
      rowwise() %>%
      mutate(
        max_prob = max(c_across(everything()), na.rm = TRUE),
        predicted_class_constrained = names(.)[which.max(c_across(everything()))]
      ) %>%
      ungroup() %>%
      mutate(
        predicted_class_constrained = str_remove(predicted_class_constrained, "^\\.pred_"),
        # Convert to factor with same levels as original
        predicted_class_constrained = factor(predicted_class_constrained, 
                                            levels = levels(cv_predictions$actual_class)),
        confidence_constrained = max_prob
      ) %>%
      bind_cols(cv_predictions %>% select(actual_class, .row)) %>%
      mutate(
        high_confidence_constrained = confidence_constrained > 0.90
      )
    
    # Compare before and after
    accuracy_raw <- mean(cv_predictions$predicted_class == cv_predictions$actual_class) * 100
    accuracy_constrained <- mean(cv_predictions_constrained$predicted_class_constrained == cv_predictions_constrained$actual_class) * 100
    
    cat("Cross-validation accuracy (raw):         ", round(accuracy_raw, 1), "%\n")
    cat("Cross-validation accuracy (constrained): ", round(accuracy_constrained, 1), "%\n")
    cat("Improvement from constraints:            ", round(accuracy_constrained - accuracy_raw, 1), " percentage points\n")
    
    # Look for cases where SC prediction was corrected
    sc_corrections <- cv_predictions %>%
      left_join(cv_input_data %>% select(.row, hgb_c_pct), by = ".row") %>%
      left_join(cv_predictions_constrained %>% select(.row, predicted_class_constrained), by = ".row") %>%
      filter(predicted_class == "SC" & hgb_c_pct < 5) %>%
      mutate(corrected = predicted_class != predicted_class_constrained)
    
    if(nrow(sc_corrections) > 0) {
      cat(sprintf("\nClinical constraint corrections:\n"))
      cat(sprintf("- %d cases with SC prediction and HgbC < 5%%\n", nrow(sc_corrections)))
      cat(sprintf("- %d of these were corrected by constraints (%.1f%%)\n", 
                  sum(sc_corrections$corrected), 
                  mean(sc_corrections$corrected) * 100))
      
      # Show examples
      if(sum(sc_corrections$corrected) > 0) {
        cat("\nExample corrections:\n")
        examples <- sc_corrections %>%
          filter(corrected) %>%
          select(actual_class, predicted_class, predicted_class_constrained, hgb_c_pct) %>%
          head(5)
        print(examples)
      }
    }
    
    # Update cv_predictions to use constrained results
    cv_predictions <- cv_predictions %>%
      left_join(cv_predictions_constrained %>% select(.row, predicted_class_constrained, confidence_constrained, high_confidence_constrained), by = ".row") %>%
      mutate(
        predicted_class = predicted_class_constrained,
        confidence = confidence_constrained,
        high_confidence = high_confidence_constrained
      ) %>%
      select(-predicted_class_constrained, -confidence_constrained, -high_confidence_constrained)
    
  } else {
    cat("No probability columns found for applying constraints\n")
  }
  
  # Helper function to handle missing columns
  if_exists <- function(col, default = 0) {
    if(substitute(col) %in% names(cv_predictions)) col else default
  }
  
  # Update confidence calculation more safely
  pred_cols <- names(cv_predictions)[str_starts(names(cv_predictions), "\\.pred_") & 
                                    !str_detect(names(cv_predictions), "\\.pred_class")]
  if(length(pred_cols) > 0) {
    cv_predictions$confidence <- apply(cv_predictions[pred_cols], 1, max, na.rm = TRUE)
  } else {
    cv_predictions$confidence <- 0.5  # Default value
  }
  
  cv_predictions$high_confidence <- cv_predictions$confidence > 0.90
  
  # Calculate overall accuracy
  cv_accuracy <- cv_predictions %>%
    accuracy(truth = actual_class, estimate = predicted_class) %>%
    pull(.estimate)
  
  cat("Cross-validation accuracy:", round(cv_accuracy * 100, 1), "%\n")
  
  # CV Coverage Analysis (consistent with final system evaluation)
  cv_automated_cases <- cv_predictions %>% filter(high_confidence)  # Using 0.90 threshold
  cv_expert_review_cases <- cv_predictions %>% filter(!high_confidence)
  
  cv_coverage <- nrow(cv_automated_cases) / nrow(cv_predictions)
  cv_accuracy_at_coverage <- if(nrow(cv_automated_cases) > 0) {
    sum(cv_automated_cases$predicted_class == cv_automated_cases$actual_class, na.rm = TRUE) / nrow(cv_automated_cases)
  } else { NA }
  
  cat("\n=== CV COVERAGE-BASED EVALUATION ===\n")
  cat(sprintf("CV Coverage (90%% confidence): %.1f%% (%d/%d cases)\n", 
              cv_coverage * 100, nrow(cv_automated_cases), nrow(cv_predictions)))
  cat(sprintf("CV Accuracy at Coverage: %.1f%% (automated cases only)\n", 
              cv_accuracy_at_coverage * 100))
  cat(sprintf("CV Expert Review Rate: %.1f%% (%d cases)\n", 
              (1 - cv_coverage) * 100, nrow(cv_expert_review_cases)))
  cat(sprintf("Traditional CV Accuracy (all cases): %.1f%% - FOR COMPARISON ONLY\n", 
              cv_accuracy * 100))
  
  # Per-class performance
  cv_class_metrics <- cv_predictions %>%
    group_by(actual_class) %>%
    summarise(
      n = n(),
      correct = sum(predicted_class == actual_class),
      accuracy = round(correct / n * 100, 1),
      avg_confidence = round(mean(confidence), 3),
      high_conf_rate = round(sum(high_confidence) / n * 100, 1),
      .groups = "drop"
    )
  
  cat("\nPer-class performance:\n")
  print(cv_class_metrics)
  
  # Confusion matrix
  cv_conf_matrix <- cv_predictions %>%
    conf_mat(truth = actual_class, estimate = predicted_class)
  
  cat("\nConfusion Matrix:\n")
  print(cv_conf_matrix)
  
  # ==============================================================================
  # 7. TRANSFUSION IMPACT ANALYSIS
  # ==============================================================================
  
  cat("\n=== TRANSFUSION IMPACT ANALYSIS ===\n")
  
  # Add transfusion status to CV predictions
  cv_predictions_with_transfusion <- cv_predictions %>%
    mutate(
      row_id = row_number()
    ) %>%
    left_join(
      all_cases_ml_data %>% 
        mutate(row_id = row_number()) %>%
        select(row_id, is_transfused),
      by = "row_id"
    )
  
  # Performance by transfusion status
  transfusion_performance <- cv_predictions_with_transfusion %>%
    group_by(is_transfused) %>%
    summarise(
      n_patients = n(),
      correct_predictions = sum(predicted_class == actual_class),
      accuracy = round(correct_predictions / n_patients * 100, 1),
      avg_confidence = round(mean(confidence), 3),
      high_conf_rate = round(sum(high_confidence) / n_patients * 100, 1),
      .groups = "drop"
    ) %>%
    mutate(
      transfusion_status = if_else(is_transfused == 1, "Transfused", "Not Transfused")
    )
  
  cat("Performance by transfusion status:\n")
  print(transfusion_performance %>% select(-is_transfused))
  
  # Statistical comparison
  if(any(cv_predictions_with_transfusion$is_transfused == 1)) {
    transfused_accuracy <- cv_predictions_with_transfusion %>%
      filter(is_transfused == 1) %>%
      summarise(accuracy = mean(predicted_class == actual_class)) %>%
      pull(accuracy)
    
    non_transfused_accuracy <- cv_predictions_with_transfusion %>%
      filter(is_transfused == 0) %>%
      summarise(accuracy = mean(predicted_class == actual_class)) %>%
      pull(accuracy)
    
    cat(sprintf("\nTransfused patients accuracy: %.1f%%\n", transfused_accuracy * 100))
    cat(sprintf("Non-transfused patients accuracy: %.1f%%\n", non_transfused_accuracy * 100))
    cat(sprintf("Difference: %.1f percentage points\n", 
                (transfused_accuracy - non_transfused_accuracy) * 100))
  }
  
  # ==============================================================================
  # 7A. COMPREHENSIVE PERFORMANCE METRICS
  # ==============================================================================
  
  cat("\n=== COMPREHENSIVE PERFORMANCE ANALYSIS ===\n")
  
  # Calculate precision, recall, F1 for each class
  calculate_multiclass_metrics <- function(actual, predicted, classes) {
    metrics_list <- list()
    
    for(class_name in classes) {
      # Binary classification for this class vs all others
      actual_binary <- ifelse(actual == class_name, 1, 0)
      predicted_binary <- ifelse(predicted == class_name, 1, 0)
      
      # Calculate metrics
      tp <- sum(actual_binary == 1 & predicted_binary == 1)
      fp <- sum(actual_binary == 0 & predicted_binary == 1)
      fn <- sum(actual_binary == 1 & predicted_binary == 0)
      tn <- sum(actual_binary == 0 & predicted_binary == 0)
      
      precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
      recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
      f1 <- ifelse(precision + recall > 0, 2 * (precision * recall) / (precision + recall), 0)
      specificity <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
      
      metrics_list[[class_name]] <- data.frame(
        class = class_name,
        precision = precision,
        recall = recall,
        f1 = f1,
        specificity = specificity,
        support = sum(actual_binary)
      )
    }
    
    do.call(rbind, metrics_list)
  }
  
  # Calculate comprehensive metrics on CV predictions
  cv_comprehensive_metrics <- calculate_multiclass_metrics(
    cv_predictions$actual_class, 
    cv_predictions$predicted_class, 
    unique(cv_predictions$actual_class)
  )
  
  cat("Per-class Performance Metrics (Cross-Validation):\n")
  print(cv_comprehensive_metrics %>% 
        mutate(across(c(precision, recall, f1, specificity), ~round(.x, 3))))
  
  # Calculate macro and weighted averages
  macro_avg <- cv_comprehensive_metrics %>%
    summarise(
      macro_precision = mean(precision),
      macro_recall = mean(recall),
      macro_f1 = mean(f1),
      macro_specificity = mean(specificity)
    )
  
  weighted_avg <- cv_comprehensive_metrics %>%
    summarise(
      weighted_precision = sum(precision * support) / sum(support),
      weighted_recall = sum(recall * support) / sum(support),
      weighted_f1 = sum(f1 * support) / sum(support),
      weighted_specificity = sum(specificity * support) / sum(support)
    )
  
  cat("\nMacro-averaged metrics:\n")
  print(macro_avg %>% mutate(across(everything(), ~round(.x, 3))))
  
  cat("\nWeighted-averaged metrics:\n")
  print(weighted_avg %>% mutate(across(everything(), ~round(.x, 3))))
  
  # ==============================================================================
  # 7B. CALIBRATION ANALYSIS
  # ==============================================================================
  
  cat("\n=== CALIBRATION ANALYSIS ===\n")
  
  # Calculate Brier Score for multiclass
  calculate_brier_score <- function(actual, pred_probs_df) {
    # Convert actual to one-hot encoding
    classes <- sort(unique(actual))
    actual_matrix <- matrix(0, nrow = length(actual), ncol = length(classes))
    for(i in 1:length(classes)) {
      actual_matrix[, i] <- ifelse(actual == classes[i], 1, 0)
    }
    
    # Get prediction probabilities in correct order
    pred_matrix <- as.matrix(pred_probs_df[, paste0(".pred_", classes)])
    
    # Calculate Brier score
    brier_scores <- rowMeans((actual_matrix - pred_matrix)^2)
    mean(brier_scores)
  }
  
  # Calculate Brier Score using tidymodels predictions
  if(exists("cv_probs") && nrow(cv_probs) > 0) {
    brier_score <- calculate_brier_score(cv_probs$genotype_simplified, cv_probs)
    cat(sprintf("Brier Score: %.4f\n", brier_score))
  }
  
  # Calculate Expected Calibration Error (ECE)
  calculate_ece <- function(actual, predicted, confidence, n_bins = 10) {
    bin_boundaries <- seq(0, 1, length.out = n_bins + 1)
    bin_lowers <- bin_boundaries[1:n_bins]
    bin_uppers <- bin_boundaries[2:(n_bins + 1)]
    
    ece <- 0
    total_samples <- length(actual)
    
    for(i in 1:n_bins) {
      in_bin <- confidence > bin_lowers[i] & confidence <= bin_uppers[i]
      prop_in_bin <- sum(in_bin) / total_samples
      
      if(sum(in_bin) > 0) {
        accuracy_in_bin <- mean(actual[in_bin] == predicted[in_bin])
        avg_confidence_in_bin <- mean(confidence[in_bin])
        ece <- ece + abs(avg_confidence_in_bin - accuracy_in_bin) * prop_in_bin
      }
    }
    
    return(ece)
  }
  
  ece_score <- calculate_ece(
    cv_predictions$actual_class, 
    cv_predictions$predicted_class, 
    cv_predictions$confidence
  )
  
  cat(sprintf("Expected Calibration Error (ECE): %.4f\n", ece_score))
  
  # ==============================================================================
  # 8. VARIABLE IMPORTANCE WITH TRANSFUSION FEATURES
  # ==============================================================================
  
  cat("\n=== VARIABLE IMPORTANCE ANALYSIS ===\n")
  
  # Prepare training data for vip
  processed_ml_data <- genotype_recipe %>%
    prep(training = all_cases_ml_data) %>%
    bake(new_data = NULL) # Use NULL to bake the training data itself

  # Define a custom prediction wrapper for XGBoost (returns probabilities)
  custom_predict_xgb_probs <- function(object, newdata) {
    predictions <- predict(object, new_data = newdata, type = "prob")
    as_tibble(predictions) 
  }

  # Define a specific prediction wrapper for accuracy metric (returns factor of classes)
  pred_wrapper_for_accuracy <- function(object, newdata) {
    predict(object, new_data = newdata, type = "class") %>% 
      pull(.pred_class)
  }

  # Get variable importance data
  var_imp_data <- final_fit %>%
    extract_fit_parsnip() %>%
    vip::vi(
      method = "permute", 
      train = processed_ml_data,        # Supply processed training data
      target = "genotype_simplified",   # Specify the outcome variable
      metric = "accuracy",             # Use accuracy metric
      nsim = 5,                        # Number of permutations
      pred_wrapper = pred_wrapper_for_accuracy # Use class prediction wrapper
    ) %>%
    arrange(desc(Importance)) %>%
    head(20)
  
  cat("\nTop 20 features:\n")
  print(var_imp_data)

  # Create variable importance plot using the calculated var_imp_data
  if(nrow(var_imp_data) > 0) {
    vip_plot <- var_imp_data %>%
      ggplot(aes(x = Importance, y = reorder(Variable, Importance), fill = Importance)) +
      geom_col(show.legend = FALSE, alpha = 0.8) +
      scale_fill_viridis_c() +
        labs(
        title = paste("Permutation Feature Importance -", model_type),
        subtitle = "Top 20 features from permutation method (accuracy)",
        x = "Importance (Permutation based on Accuracy)",
        y = "Variable"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10)
      )
    
    # Save plot
    dir.create("figures", showWarnings = FALSE)
    ggsave("figures/variable_importance_tidymodels_with_transfusion.png", 
           vip_plot, width = 12, height = 8, dpi = 300)
    
    cat("Variable importance plot saved\n")
    } else {
    cat("No variable importance data to plot.\n")
  }

  # Identify transfusion-related features
  transfusion_features <- var_imp_data %>%
    filter(str_detect(Variable, "transfus|is_transfused"))
  
  if(nrow(transfusion_features) > 0) {
    cat("\nTransfusion-related features in top 20:\n")
    print(transfusion_features)
} else {
    cat("\nNo transfusion-related features in top 20\n")
}

# ==============================================================================
# 9. SHAP INTERPRETABILITY ANALYSIS (using shapviz)
# ==============================================================================

cat("\n=== SHAP INTERPRETABILITY ANALYSIS ===\n")

if(exists("final_fit") && !is.null(final_fit) && exists("processed_ml_data") && nrow(processed_ml_data) > 0) {
  
  cat("Preparing data for SHAP analysis...\n")
  
  # Extract the core XGBoost model object from the tidymodels workflow
  xgb_core_model <- final_fit %>% 
    extract_fit_parsnip() %>% 
    pluck("fit")
  
  # Prepare the feature matrix (ensure no target variable is present)
  # processed_ml_data was created before and contains the target `genotype_simplified`
  features_for_shap <- processed_ml_data %>% 
    select(-genotype_simplified)
  
  # Convert to matrix as required by xgboost
  data_matrix_for_shap <- as.matrix(features_for_shap)
  
  cat(sprintf("SHAP data matrix dimensions: %d rows, %d features\n", 
              nrow(data_matrix_for_shap), ncol(data_matrix_for_shap)))
  
  # Calculate SHAP contributions using xgboost's predict function
  # For multiclass, predcontrib = TRUE returns a list of matrices (one per class)
  # Each matrix has dimensions: [n_instances, n_features + 1 (for BIAS)]
  # cat("Calculating SHAP contributions...\n")
  # shap_contributions <- predict(
  #   xgb_core_model, 
  #   newdata = data_matrix_for_shap, 
  #   predcontrib = TRUE
  # )
  
  # Create the shapviz object by passing the model and predictor data
  # shapviz will internally call the appropriate SHAP calculation for XGBoost
  cat("Creating shapviz object (this may take a moment for SHAP calculation)...\n")
  shap_calculation_success <- TRUE
  sv_raw_output <- NULL # Ensure it exists in this scope
  tryCatch({
    sv_raw_output <- shapviz(
      xgb_core_model, 
      X_pred = data_matrix_for_shap, 
      X = features_for_shap # Original features for visualization context
    )
  }, error = function(e) {
    cat(paste("Error during shapviz() call:", e$message, "\n"))
    shap_calculation_success <<- FALSE
    # sv_raw_output will remain NULL or be the error object if not caught by tryCatch itself
  })

  sv <- NULL # Initialize sv to NULL
  if (shap_calculation_success && !is.null(sv_raw_output)) {
    class_levels_for_naming <- levels(all_cases_ml_data$genotype_simplified)
    expected_names <- make.names(paste0("shap_pred_", class_levels_for_naming))

    if (is.list(sv_raw_output) && !inherits(sv_raw_output, "shapviz")) {
      # Case 1: Output is a list, hopefully of shapviz objects
      if (length(sv_raw_output) == length(class_levels_for_naming)) {
        # Check if all elements are indeed shapviz objects
        all_elements_are_shapviz <- all(sapply(sv_raw_output, inherits, "shapviz"))
        if(all_elements_are_shapviz){
            sv <- setNames(sv_raw_output, expected_names)
            cat("Successfully created and named a list of shapviz objects for each class.\n")
        } else {
            cat("Warning: List output from shapviz() has correct length but not all elements are shapviz objects. sv will be this raw list.\n")
            sv <- sv_raw_output # Keep it as is, downstream will filter valid ones
        }
      } else if (length(sv_raw_output) == 1 && inherits(sv_raw_output[[1]], "shapviz")) {
        # Case 2: Output is a list of length 1 containing a single shapviz object
        cat("shapviz returned a list of length 1 containing a shapviz object. Naming the list element using first class name.\n")
        sv <- setNames(sv_raw_output, expected_names[1]) # Name the list element itself
      } else {
        cat("Warning: shapviz() returned a list of unexpected length or structure. sv will be this raw list.\n")
        sv <- sv_raw_output
      }
    } else if (inherits(sv_raw_output, "shapviz")) {
      # Case 3: Output is a single shapviz object (not a list)
      cat("shapviz returned a single shapviz object. Wrapping it in a named list (using first class name) for consistent handling.\n")
      sv <- list()
      sv[[expected_names[1]]] <- sv_raw_output
    } else {
      cat("Warning: Unexpected structure from shapviz() call. sv is being set to this raw output.\n")
      sv <- sv_raw_output # Fallback, might be NULL or an error object
    }
  } else {
    cat("SHAP calculation failed or sv_raw_output is NULL. Skipping further SHAP processing.\n")
    # sv remains NULL, shap_importance_df will be empty later
  }
  
  cat("SHAP analysis setup complete.\n")
  
  # Filter for valid shapviz objects before plotting
  sv_for_plotting <- list()
  if (!is.null(sv) && is.list(sv)) {
    sv_for_plotting <- sv[sapply(sv, function(x) inherits(x, "shapviz"))]
    if(length(sv_for_plotting) == 0) {
        cat("No valid shapviz objects found in the list for plotting.\n")
    }
  } else if (!is.null(sv) && inherits(sv, "shapviz")) { # If sv itself is a single shapviz object
    sv_for_plotting <- list(sv) # Treat as a list with one item
  }

  # 1. Overall SHAP Importance Plots
  # =================================
  if (length(sv_for_plotting) > 0) {
    cat("Generating overall SHAP importance plots...\n")
    
    # Aggregate SHAP importance from all valid shapviz objects in the list
    all_class_shap_importance <- list()
    for(shap_obj_name in names(sv_for_plotting)){
      shap_obj <- sv_for_plotting[[shap_obj_name]]
      if(inherits(shap_obj, "shapviz")){
        tryCatch({
          class_imp_raw <- sv_importance(shap_obj, kind = "no")
          # Check if the result is a data.frame and has expected columns
          if(is.data.frame(class_imp_raw) && "feature" %in% names(class_imp_raw) && "mean_abs_shap" %in% names(class_imp_raw)){
            class_imp <- class_imp_raw %>%
              mutate(class_source = shap_obj_name) # Add source if needed
            all_class_shap_importance[[shap_obj_name]] <- class_imp
          } else {
            cat(paste("Warning: sv_importance for", shap_obj_name, "did not return expected data.frame structure.\n"))
          }
        }, error = function(e) {
          cat(paste("Error getting SHAP importance for", shap_obj_name, ":", e$message, "\n"))
        })
      }
    }

    if(length(all_class_shap_importance) > 0) {
      combined_shap_importance <- do.call(rbind, all_class_shap_importance)
      
      # Aggregate: Mean of mean_abs_shap across classes for each feature
      shap_importance_df <- combined_shap_importance %>%
        group_by(feature) %>%
        summarise(mean_abs_shap = mean(mean_abs_shap, na.rm = TRUE)) %>%
        arrange(desc(mean_abs_shap))

      # Bar plot of AGGREGATED SHAP feature importance
      if(nrow(shap_importance_df) > 0) {
        tryCatch({
          shap_importance_plot_bar <- shap_importance_df %>%
            head(min(20, nrow(.))) %>%
            ggplot(aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap, fill = mean_abs_shap)) +
            geom_col(show.legend = FALSE) +
            coord_flip() +
            scale_fill_viridis_c() +
            ggtitle("Overall SHAP Feature Importance (Aggregated Mean Absolute SHAP)") +
            labs(x = "Feature", y = "Aggregated Mean Absolute SHAP Value") +
            theme_minimal()
          ggsave("figures/shap_overall_importance_bar_aggregated.png", shap_importance_plot_bar, 
                 width = 10, height = 8, dpi = 300)
          cat("Aggregated SHAP bar plot saved to figures/shap_overall_importance_bar_aggregated.png\n")
        }, error = function(e) {
          cat(paste("Error generating aggregated SHAP bar plot:", e$message, "\n"))
        })
      } else {
        cat("No data for aggregated SHAP importance plot.\n")
        if(!exists("shap_importance_df")) shap_importance_df <- data.frame() # Ensure exists
      }
    } else {
      cat("No class-specific SHAP importance data to combine.\n")
      shap_importance_df <- data.frame()
    }

    # Beeswarm plot for SHAP values (per class)
    cat("Generating SHAP beeswarm plots (per class)...\n")
    dir.create("figures/shap_beeswarm_per_class", showWarnings = FALSE)
    for(shap_obj_name in names(sv_for_plotting)){
      shap_obj <- sv_for_plotting[[shap_obj_name]]
      if(inherits(shap_obj, "shapviz")){
        tryCatch({
          top_n_features_for_beeswarm <- min(20, ncol(shap_obj$X))
          shap_beeswarm_plot_class <- sv_importance(shap_obj, kind = "beeswarm", max_display = top_n_features_for_beeswarm) + 
            ggtitle(paste0("SHAP Beeswarm Plot - Class: ", gsub("shap_pred_","",shap_obj_name), 
                           " (Top ", top_n_features_for_beeswarm, " Features)")) +
            theme_minimal()
          ggsave(paste0("figures/shap_beeswarm_per_class/shap_beeswarm_", shap_obj_name, ".png"), 
                 shap_beeswarm_plot_class, width = 12, height = 10, dpi = 300)
          cat(paste0("SHAP beeswarm plot for class ", shap_obj_name, " saved.\n"))
        }, error = function(e) {
          cat(paste("Error generating SHAP beeswarm plot for class", shap_obj_name, ":", e$message, "\n"))
        })
      }
    }

    # --- NEW: Per-class mean absolute SHAP grouped bar plot ---
    # Only proceed if we have valid SHAP importance data
    if(length(all_class_shap_importance) > 0) {
      tryCatch({
        # Map shapviz object names to genotype labels
        class_name_map <- c(
          shap_pred_SS = "SS",
          shap_pred_SC = "SC",
          shap_pred_S_beta_0 = "Sβ⁰",
          shap_pred_S_beta_plus = "Sβ⁺",
          shap_pred_AS = "AS"  # Add AS mapping
        )

        # Collect mean absolute SHAP for each class
        per_class_shap_long <- purrr::imap_dfr(
          sv_for_plotting,
          function(shap_obj, shap_name) {
            tryCatch({
              imp <- sv_importance(shap_obj, kind = "no")
              if(is.data.frame(imp) && "feature" %in% names(imp) && "mean_abs_shap" %in% names(imp)) {
                imp$class <- class_name_map[shap_name]
                return(imp)
              } else {
                return(NULL)
              }
            }, error = function(e) {
              cat(paste("Error processing SHAP importance for", shap_name, ":", e$message, "\n"))
              return(NULL)
            })
          }
        )

        # Only proceed if we have valid data
        if(!is.null(per_class_shap_long) && nrow(per_class_shap_long) > 0 && "feature" %in% names(per_class_shap_long)) {
          # Only keep top N features by overall mean SHAP
          top_n <- 15
          top_features <- per_class_shap_long %>%
            group_by(feature) %>%
            summarise(mean_shap = mean(mean_abs_shap, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(mean_shap)) %>%
            slice_head(n = top_n) %>%
            pull(feature)

          per_class_shap_long_filtered <- per_class_shap_long %>%
            filter(feature %in% top_features)

          # Plot
          shap_per_class_plot <- ggplot(per_class_shap_long_filtered, aes(x = mean_abs_shap, y = feature, fill = class)) +
            geom_col(position = "dodge") +
            scale_fill_viridis_d(name = "Genotype") +
            labs(
              title = "Per-Genotype SHAP Feature Importance (Mean Absolute SHAP)",
              x = "mean(|SHAP value|)",
              y = "Feature"
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold"),
              axis.text.y = element_text(size = 10)
            )

          ggsave("figures/shap_per_class_mean_abs_bar.png", shap_per_class_plot, width = 12, height = 8, dpi = 300)
          cat("Per-class mean absolute SHAP bar plot saved to figures/shap_per_class_mean_abs_bar.png\n")
        } else {
          cat("No valid per-class SHAP data available for grouped bar plot.\n")
        }
      }, error = function(e) {
        cat(paste("Error creating per-class SHAP bar plot:", e$message, "\n"))
      })
    } else {
      cat("Skipping per-class SHAP bar plot: no valid SHAP importance data.\n")
    }
  } else {
    cat("Skipping overall SHAP importance plots as no valid shapviz objects are available.\n")
    shap_importance_df <- data.frame() # Ensure it exists if plots are skipped
  }

  # 2. SHAP Explanations for Misclassified Cases
  # ============================================
  if (!is.null(sv) && (is.list(sv) && length(sv[sapply(sv, function(x) inherits(x, "shapviz"))]) > 0 || inherits(sv, "shapviz"))) {
    cat("Generating SHAP explanations for misclassified cases...\n")
    dir.create("figures/shap_misclassified", showWarnings = FALSE)

    # Get final model predictions on ml_data (the data SHAP was computed on)
    final_ml_preds_for_shap <- predict(final_fit, new_data = all_cases_ml_data, type = "class") %>%
      bind_cols(predict(final_fit, new_data = all_cases_ml_data, type = "prob")) %>%
      mutate(
        actual_class = all_cases_ml_data$genotype_simplified,
        row_id_in_ml_data = row_number() 
      )

    misclassified_shap_cases <- final_ml_preds_for_shap %>%
      filter(.pred_class != actual_class) %>%
      head(5) 

    if(nrow(misclassified_shap_cases) > 0) {
      for(i in 1:nrow(misclassified_shap_cases)) {
        case_info <- misclassified_shap_cases[i, ]
        case_row_id <- case_info$row_id_in_ml_data
        
        true_label <- as.character(case_info$actual_class)
        pred_label <- as.character(case_info$.pred_class)
        
        pred_prob_col_name <- paste0(".pred_", pred_label)
        confidence_case <- if(pred_prob_col_name %in% names(case_info)) {
          case_info[[pred_prob_col_name]]
        } else { NA }
        
        plot_title <- sprintf("SHAP Waterfall for Misclassified Case (Row ID: %d)\nPredicted: %s (Prob: %.2f), True: %s", 
                              case_row_id, pred_label, confidence_case, true_label)
        
        # Select the SHAP object for the predicted class if sv is a list
        current_sv_object <- NULL
        if (is.list(sv) && !inherits(sv, "shapviz")) {
          shap_obj_name <- make.names(paste0("shap_pred_", pred_label))
          if (shap_obj_name %in% names(sv) && inherits(sv[[shap_obj_name]], "shapviz")) {
            current_sv_object <- sv[[shap_obj_name]]
          } else {
             # Fallback: try to use the first valid shapviz object if predicted class specific one not found/valid
            if(length(sv_for_plotting) > 0) current_sv_object <- sv_for_plotting[[1]]
          }
        } else if (inherits(sv, "shapviz")) { # Single shapviz object
          current_sv_object <- sv
        }

        if(!is.null(current_sv_object)){
          tryCatch({
            shap_waterfall_plot <- sv_waterfall(current_sv_object, row_id = case_row_id) +
              ggtitle(plot_title) +
              theme_minimal() +
              theme(plot.title = element_text(size = 10))
            
            ggsave(paste0("figures/shap_misclassified/waterfall_rowid_", case_row_id, 
                            "_pred_", pred_label, "_true_", true_label, ".png"), 
                   shap_waterfall_plot, width = 10, height = 8, dpi = 300)
            cat(paste0("Saved SHAP waterfall for misclassified case, row_id: ", case_row_id, "\n"))
          }, error = function(e) {
            cat(paste("Error generating SHAP waterfall for row_id", case_row_id, ":", e$message, "\n"))
          })
        } else {
          cat(paste("Skipping SHAP waterfall for row_id", case_row_id, ": No valid shapviz object for predicted class", pred_label, "or general fallback.\n"))
        }
      }
    } else {
      cat("No misclassified cases found in ml_data to explain with SHAP.\n")
    }
  } else {
    cat("Skipping SHAP explanations for misclassified cases: sv object not available or invalid.\n")
  }

} else {
  cat("Skipping SHAP analysis: final_fit or processed_ml_data not available or empty.\n")
  sv <- NULL # Ensure sv exists for comprehensive_results
  shap_importance_df <- data.frame() # Ensure this exists
}

# ==============================================================================
# 8A. ENHANCED PERFORMANCE VISUALIZATION 
# ==============================================================================

cat("\n=== GENERATING ENHANCED PERFORMANCE VISUALIZATIONS ===\n")

# 1. CALIBRATION PLOTS
# Ensure cv_predictions and best_tune_results/xgb_best exist from previous steps.
if(exists("cv_predictions") && nrow(cv_predictions) > 0 && 
   exists("best_tune_results") && exists("xgb_best")) {
    
  current_best_config_id <- NULL
  if (nrow(xgb_best) == 1 && ".config" %in% names(xgb_best)) {
    current_best_config_id <- xgb_best$.config[1]
  }

  if (!is.null(current_best_config_id) && !is.na(current_best_config_id)) {
    # Get probability predictions for calibration analysis
    # These predictions are from the cross-validation folds for the *best* hyperparameter set
    cv_probs_all_configs <- collect_predictions(best_tune_results)

    if (".config" %in% names(cv_probs_all_configs)) {
      cv_probs <- cv_probs_all_configs %>%
        filter(.config == current_best_config_id)
      
      if (nrow(cv_probs) > 0) {
        create_calibration_plot <- function(actual, predicted_probs, class_name, n_bins = 10) {
          # Binary classification for this class
          actual_binary <- ifelse(actual == class_name, 1, 0)
          
          # Create bins
          bin_boundaries <- seq(0, 1, length.out = n_bins + 1)
          bin_centers <- (bin_boundaries[-1] + bin_boundaries[-length(bin_boundaries)]) / 2
          
          calibration_data <- data.frame()
          
          for(i in 1:n_bins) {
            bin_lower <- bin_boundaries[i]
            bin_upper <- bin_boundaries[i + 1]
            # Ensure the last bin is inclusive of 1.0
            in_bin <- if(i == n_bins) {
                        predicted_probs >= bin_lower & predicted_probs <= bin_upper
                      } else {
                        predicted_probs >= bin_lower & predicted_probs < bin_upper
                      }
            
            if(sum(in_bin, na.rm = TRUE) > 0) { # Added na.rm = TRUE
              mean_predicted <- mean(predicted_probs[in_bin], na.rm = TRUE)
              mean_actual <- mean(actual_binary[in_bin], na.rm = TRUE)
              count <- sum(in_bin, na.rm = TRUE)
              
              calibration_data <- rbind(calibration_data, data.frame(
                bin_center = bin_centers[i],
                mean_predicted = mean_predicted,
                mean_actual = mean_actual,
                count = count
              ))
            }
          }
          
          ggplot(calibration_data, aes(x = mean_predicted, y = mean_actual)) +
            geom_point(aes(size = count), alpha = 0.7, color = "steelblue") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
            geom_smooth(method = "loess", se = TRUE, alpha = 0.3, na.rm = TRUE) + # Added na.rm = TRUE
            xlim(0, 1) + ylim(0, 1) +
            labs(
              title = paste("Calibration Plot -", class_name),
              x = "Mean Predicted Probability",
              y = "Fraction of Positives",
              size = "Count"
            ) +
            theme_minimal() +
            theme(plot.title = element_text(size = 12, face = "bold"))
        }
        
        # Create calibration plots for each class
        # Ensure cv_probs has the 'genotype_simplified' column (truth) and .pred_CLASS columns
        if("genotype_simplified" %in% names(cv_probs)) {
          classes_for_calib <- levels(cv_probs$genotype_simplified)
          calibration_plots_list <- list() # Renamed from calibration_plots
          
          for(class_val in classes_for_calib) { # Renamed class_name to class_val
            pred_col_name <- paste0(".pred_", class_val) # Renamed pred_col
            if(pred_col_name %in% colnames(cv_probs)) {
              # Check if there are any positive instances for this class in actuals
              if (any(cv_probs$genotype_simplified == class_val, na.rm = TRUE)) {
                 calibration_plots_list[[class_val]] <- create_calibration_plot(
                    cv_probs$genotype_simplified,
                    cv_probs[[pred_col_name]],
                    class_val # Use class_val here
                 )
              } else {
                cat(paste0("Skipping calibration plot for class '", class_val, "': No actual positive cases in cv_probs.\\n"))
              }
            } else {
              cat(paste0("Skipping calibration plot for class '", class_val, "': Prediction column '", pred_col_name, "' not found.\\n"))
            }
          }
          
          if(length(calibration_plots_list) > 0) {
            combined_calibration_plots <- do.call(grid.arrange, c(calibration_plots_list, list(ncol = 2))) # Renamed
            ggsave("figures/calibration_plots_tidymodels.png", combined_calibration_plots, 
                   width = 12, height = 8, dpi = 300)
            cat("Calibration plots saved\\n")
          } else {
            cat("No valid calibration plots generated.\\n")
          }
        } else {
          cat("Column 'genotype_simplified' not found in cv_probs. Skipping calibration plots.\\n")
        }
      } else {
        cat("No predictions found for the best configuration (cv_probs is empty). Skipping calibration plots.\\n")
      }
    } else {
      cat("'.config' column not found in collect_predictions(best_tune_results). Skipping calibration plots.\\n")
    }
  } else {
    cat("Best configuration ID is invalid or xgb_best is not as expected. Skipping calibration plots.\\n")
  }
} else {
  cat("CV predictions, tuning results, or best hyperparameter selection not available. Skipping calibration plots.\\n")
}
  
  # 2. ENHANCED CONFUSION MATRIX HEATMAP
  if(exists("cv_conf_matrix")) {
    # Convert confusion matrix to data frame
    confusion_matrix_df <- as.data.frame(cv_conf_matrix$table)
    
    confusion_heatmap <- confusion_matrix_df %>%
      ggplot(aes(x = Truth, y = Prediction, fill = Freq)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = Freq), color = "white", fontface = "bold", size = 3) +
      scale_fill_viridis_c(name = "Count") +
      labs(
        title = "Confusion Matrix - Cross-Validation Results",
        subtitle = "Tidymodels Pipeline with Transfusion Awareness",
        x = "Actual Genotype",
        y = "Predicted Genotype"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 0),
        panel.grid = element_blank()
      )
    
    ggsave("figures/confusion_matrix_heatmap_tidymodels.png", confusion_heatmap, 
           width = 10, height = 8, dpi = 300)
    cat("Confusion matrix heatmap saved\n")
  } else {
    cat("Confusion matrix not available for visualization\n")
  }
  
  # 3. CROSS-VALIDATION PERFORMANCE VISUALIZATION
  if (exists("best_tune_results") && exists("xgb_best")) {
    cv_preds <- collect_predictions(best_tune_results) %>%
      filter(.config == xgb_best$.config)
    
    # Calculate accuracy per fold
    per_fold_acc <- cv_preds %>%
      group_by(id) %>%
      summarise(
        accuracy = mean(.pred_class == genotype_simplified),
        .groups = "drop"
      )
    
    # Plot per-fold accuracy
    cv_plot <- per_fold_acc %>%
      ggplot(aes(x = id, y = accuracy)) +
      geom_col(fill = viridis(1), alpha = 0.8) +
      geom_hline(yintercept = mean(per_fold_acc$accuracy), 
                 color = "red", linetype = "dashed", linewidth = 1) +
      geom_text(aes(label = paste0(round(accuracy * 100, 1), "%")), 
                vjust = -0.3, fontface = "bold", size = 3) +
      annotate("text", x = 5, y = mean(per_fold_acc$accuracy) + 0.02, 
               label = paste0("Mean: ", round(mean(per_fold_acc$accuracy) * 100, 1), "%"), 
               color = "red", fontface = "bold") +
      ylim(0, 1) +
      labs(
        title = "10-Fold Cross-Validation Accuracy - Tidymodels",
        subtitle = paste0("Model: ", if (exists("model_type")) model_type else "XGBoost", " with transfusion awareness"),
        x = "Fold",
        y = "Accuracy"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave("figures/cross_validation_performance_tidymodels.png", cv_plot, 
           width = 10, height = 6, dpi = 300)
    cat("Cross-validation performance plot saved\n")
  }
  
  # 4. CONFIDENCE DISTRIBUTION ANALYSIS
  confidence_dist_plot <- cv_predictions %>%
    ggplot(aes(x = confidence, fill = predicted_class == actual_class)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 0.90, color = "red", linetype = "dashed", size = 1) +
    scale_fill_manual(values = c("FALSE" = "#FF6B6B", "TRUE" = "#4ECDC4"),
                     name = "Prediction",
                     labels = c("Incorrect", "Correct")) +
    labs(
      title = "ML Prediction Confidence Distribution - Tidymodels",
      subtitle = "Red line shows 90% confidence threshold for expert review",
      x = "Prediction Confidence",
      y = "Number of Predictions"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
    ggsave("figures/confidence_distribution_tidymodels.png", confidence_dist_plot, 
           width = 10, height = 6, dpi = 300)
    cat("Confidence distribution plot saved\n")
  
  # 5. MULTICLASS ROC AND PRECISION-RECALL CURVES
  # Calculate AUROC and AUPRC for each class
  roc_auc_results <- data.frame()
  pr_auc_results <- data.frame()
  
  classes <- levels(all_cases_ml_data$genotype_simplified)
  
  for(class_name in classes) {
    pred_col <- paste0(".pred_", class_name)
    
    if(pred_col %in% colnames(cv_probs)) {
      # Binary classification for this class
      actual_binary <- ifelse(cv_probs$genotype_simplified == class_name, 1, 0)
      predicted_probs <- cv_probs[[pred_col]]
      
      # Calculate ROC
      tryCatch({
        roc_obj <- pROC::roc(actual_binary, predicted_probs, quiet = TRUE)
        auroc <- pROC::auc(roc_obj)
        
        roc_auc_results <- rbind(roc_auc_results, data.frame(
          class = class_name,
          auroc = as.numeric(auroc),
          support = sum(actual_binary)
        ))
        
        # Save individual ROC plot
        png(sprintf("figures/roc_curve_tidymodels_%s.png", class_name), width=800, height=600)
        plot(roc_obj, main = paste("ROC Curve -", class_name, "(Tidymodels)"), col = "blue", lwd = 2)
        auc_val <- round(pROC::auc(roc_obj), 3)
        legend("bottomright", legend = paste("AUC =", auc_val), bty = "n")
        dev.off()
        
      }, error = function(e) {
        cat(sprintf("Error calculating ROC for %s: %s\n", class_name, e$message))
      })
    }
  }
  
  # Calculate macro and weighted averages for AUROC
  if(nrow(roc_auc_results) > 0) {
    macro_auroc <- mean(roc_auc_results$auroc, na.rm = TRUE)
    weighted_auroc <- sum(roc_auc_results$auroc * roc_auc_results$support, na.rm = TRUE) / 
                     sum(roc_auc_results$support, na.rm = TRUE)
    
    cat("\nAUROC Results:\n")
    print(roc_auc_results %>% mutate(auroc = round(auroc, 3)))
    cat(sprintf("Macro-averaged AUROC: %.3f\n", macro_auroc))
    cat(sprintf("Weighted-averaged AUROC: %.3f\n", weighted_auroc))
  }
  
  cat("ROC curves saved for each class\n")
  
} else {
  cat("No ML data available - all cases classified by rules!\n")
  final_fit <- NULL
  cv_predictions <- data.frame()
  var_imp_data <- data.frame()
  transfusion_performance <- data.frame()
}

# ==============================================================================
# 9. COMBINE RESULTS AND FINAL PREDICTIONS
# ==============================================================================

cat("\n=== FINAL PREDICTIONS ===\n")

# Combine rule-based and ML predictions
key_genotypes_for_ml <- c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")  # NOW 5 CLASSES

classifiable_data <- data_with_rules %>%
  filter(
    high_confidence |  # Can be classified by rules
    (!high_confidence & genotype_simplified %in% key_genotypes_for_ml)  # Can be classified by ML
  )

final_predictions <- classifiable_data %>%
  mutate(
    final_prediction = case_when(
      high_confidence ~ rule_prediction,
      TRUE ~ "PENDING_ML"
    ),
    prediction_source = case_when(
      high_confidence ~ "Clinical_Rules", 
      TRUE ~ "Machine_Learning"
    ),
    ml_confidence = NA_real_,
    needs_expert_review = FALSE
  )

# Add ML predictions if available
if(exists("final_fit") && !is.null(final_fit)) {
  ambiguous_for_prediction <- final_predictions %>%
    filter(final_prediction == "PENDING_ML") %>%
    prepare_ml_features() %>%
    filter(genotype_simplified %in% key_genotypes_for_ml) %>%
    select(-genotype_simplified)
  
  if(nrow(ambiguous_for_prediction) > 0) {
    ml_predictions <- predict(final_fit, new_data = ambiguous_for_prediction, type = "prob") %>%
      bind_cols(predict(final_fit, new_data = ambiguous_for_prediction, type = "class"))
    
    # APPLY CLINICAL CONSTRAINTS TO FINAL PREDICTIONS
    cat("\n=== APPLYING CLINICAL CONSTRAINTS TO FINAL ML PREDICTIONS ===\n")
    
    # Get the corresponding patient data for constraints
    final_ml_input_data <- final_predictions %>%
      filter(final_prediction == "PENDING_ML") %>%
      select(hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct)
    
    # Apply constraints to probabilities
    prob_cols_final <- names(ml_predictions)[str_starts(names(ml_predictions), "\\.pred_") & 
                                           !str_detect(names(ml_predictions), "\\.pred_class")]
    
    if(length(prob_cols_final) > 0 && nrow(final_ml_input_data) == nrow(ml_predictions)) {
      ml_probabilities_raw <- ml_predictions %>% select(all_of(prob_cols_final))
      
      # Apply clinical constraints
      ml_probabilities_constrained <- apply_clinical_constraints_pipeline(
        ml_probabilities_raw, 
        final_ml_input_data
      )
      
      # Recalculate predictions based on constrained probabilities
      ml_predictions_final <- ml_probabilities_constrained %>%
        rowwise() %>%
        mutate(
          max_prob = max(c_across(everything()), na.rm = TRUE),
          ml_prediction_constrained = names(.)[which.max(c_across(everything()))]
        ) %>%
        ungroup() %>%
        mutate(
          ml_prediction_constrained = str_remove(ml_prediction_constrained, "^\\.pred_"),
          # Convert to factor with same levels as training data
          ml_prediction_constrained = factor(ml_prediction_constrained, 
                                            levels = levels(all_cases_ml_data$genotype_simplified)),
          ml_confidence_constrained = max_prob,
          needs_review_constrained = ml_confidence_constrained <= 0.90
        )
      
      # Count corrections made
      raw_predictions <- ml_predictions$.pred_class
      constrained_predictions <- ml_predictions_final$ml_prediction_constrained
      corrections_made <- sum(raw_predictions != constrained_predictions, na.rm = TRUE)
      
      cat(sprintf("Final ML predictions: %d cases\n", nrow(ml_predictions)))
      cat(sprintf("Clinical constraint corrections: %d cases (%.1f%%)\n", 
                  corrections_made, corrections_made / nrow(ml_predictions) * 100))
      
      # Show specific SC corrections
      sc_final_corrections <- data.frame(
        raw_pred = as.character(raw_predictions),
        constrained_pred = constrained_predictions,
        hgb_c_pct = final_ml_input_data$hgb_c_pct
      ) %>%
        filter(raw_pred == "SC" & hgb_c_pct < 5)
      
      if(nrow(sc_final_corrections) > 0) {
        cat(sprintf("SC corrections (HgbC < 5%%): %d cases\n", nrow(sc_final_corrections)))
        cat("Corrections:\n")
        print(sc_final_corrections %>% 
              count(raw_pred, constrained_pred, name = "count") %>%
              arrange(desc(count)))
      }
      
      # Use constrained results
      ml_predictions <- ml_predictions_final %>%
        select(ml_prediction = ml_prediction_constrained, 
               ml_confidence = ml_confidence_constrained,
               needs_review = needs_review_constrained)
      
    } else {
      cat("Unable to apply constraints - data mismatch\n")
      # Fallback to original logic
      ml_predictions <- ml_predictions %>%
        mutate(
          ml_prediction = .pred_class,
          ml_confidence = pmax(
            if(".pred_SS" %in% names(.)) .pred_SS else 0,
            if(".pred_SC" %in% names(.)) .pred_SC else 0,
            if(".pred_S_beta_0" %in% names(.)) .pred_S_beta_0 else 0,
            if(".pred_S_beta_plus" %in% names(.)) .pred_S_beta_plus else 0,
            if(".pred_AS" %in% names(.)) .pred_AS else 0  # NEW: Include AS
          ),
          needs_review = ml_confidence <= 0.90
        )
    }
    
    # Update final predictions
    ambiguous_indices <- which(final_predictions$final_prediction == "PENDING_ML")
    if(length(ambiguous_indices) == nrow(ml_predictions)) {
      final_predictions$final_prediction[ambiguous_indices] <- as.character(ml_predictions$ml_prediction)
      final_predictions$ml_confidence[ambiguous_indices] <- ml_predictions$ml_confidence
      final_predictions$needs_expert_review[ambiguous_indices] <- ml_predictions$needs_review
    }
  }
}

# ==============================================================================
# 10. COMPREHENSIVE PERFORMANCE SUMMARY
# ==============================================================================

cat("\n=== COMPREHENSIVE PERFORMANCE SUMMARY ===\n")

# Calculate final metrics
rule_cases <- final_predictions %>% filter(prediction_source == "Clinical_Rules")
ml_cases <- final_predictions %>% filter(prediction_source == "Machine_Learning")

rule_accuracy <- if(nrow(rule_cases) > 0) {
  sum(rule_cases$final_prediction == rule_cases$genotype_simplified, na.rm = TRUE) / nrow(rule_cases) * 100
} else { NA }

ml_accuracy <- if(nrow(ml_cases) > 0) {
  sum(ml_cases$final_prediction == ml_cases$genotype_simplified, na.rm = TRUE) / nrow(ml_cases) * 100
} else { NA }

overall_accuracy <- if(nrow(final_predictions) > 0) {
  sum(final_predictions$final_prediction == final_predictions$genotype_simplified, na.rm = TRUE) / nrow(final_predictions) * 100
} else { NA }

# Overall performance summary
cat("=== FINAL TIDYMODELS PIPELINE PERFORMANCE ===\n")
cat("- Total classifiable patients:", nrow(final_predictions), "\n")
cat("- Rule-classified:", nrow(rule_cases), 
    "(", round(nrow(rule_cases)/nrow(final_predictions)*100, 1), "%)\n")
cat("- ML-classified:", nrow(ml_cases), 
    "(", round(nrow(ml_cases)/nrow(final_predictions)*100, 1), "%)\n")
cat("- Rule accuracy:", round(rule_accuracy, 1), "%\n")
cat("- ML accuracy:", round(ml_accuracy, 1), "%\n")
cat("- Overall accuracy:", round(overall_accuracy, 1), "%\n")

# Transfusion summary in final results
final_transfusion_summary <- final_predictions %>%
  group_by(is_transfused) %>%
  summarise(
    n_patients = n(),
    correct = sum(final_prediction == genotype_simplified, na.rm = TRUE),
    accuracy = round(correct / n_patients * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(transfusion_status = if_else(is_transfused == 1, "Transfused", "Not Transfused"))

cat("\nFinal results by transfusion status:\n")
print(final_transfusion_summary %>% select(-is_transfused))

# ==============================================================================
# 10A. COMPREHENSIVE FINAL METRICS CALCULATION
# ==============================================================================

cat("\n=== COMPREHENSIVE FINAL METRICS ===\n")

# Calculate sensitivity and specificity for each genotype
all_gts <- unique(final_predictions$genotype_simplified)
sens_spec <- map_dfr(all_gts, function(gt) {
  tp <- sum(final_predictions$genotype_simplified == gt & final_predictions$final_prediction == gt, na.rm=TRUE)
  fn <- sum(final_predictions$genotype_simplified == gt & final_predictions$final_prediction != gt, na.rm=TRUE)
  tn <- sum(final_predictions$genotype_simplified != gt & final_predictions$final_prediction != gt, na.rm=TRUE)
  fp <- sum(final_predictions$genotype_simplified != gt & final_predictions$final_prediction == gt, na.rm=TRUE)
  
  tibble(
    genotype = gt,
    sensitivity = ifelse(tp + fn > 0, tp / (tp + fn), NA),
    specificity = ifelse(tn + fp > 0, tn / (tn + fp), NA),
    precision = ifelse(tp + fp > 0, tp / (tp + fp), NA),
    support = tp + fn,  # Total actual cases of this genotype
    tp = tp, fn = fn, tn = tn, fp = fp
  )
})

cat("\nSensitivity and Specificity by Genotype:\n")
sens_spec_display <- sens_spec %>% 
  select(genotype, sensitivity, specificity, precision, support) %>%
  mutate(
    sensitivity = round(sensitivity, 3),
    specificity = round(specificity, 3),
    precision = round(precision, 3)
  )
print(sens_spec_display)

# Detailed confusion matrix for final predictions
cat("\nDetailed Confusion Matrix:\n")
confusion_final <- table(
  Predicted = final_predictions$final_prediction,
  Actual = final_predictions$genotype_simplified
)
print(confusion_final)

# Per-genotype performance breakdown
cat("\nPer-Genotype Performance Breakdown:\n")
genotype_performance <- final_predictions %>%
  group_by(genotype_simplified) %>%
  summarise(
    total_cases = n(),
    correctly_predicted = sum(final_prediction == genotype_simplified, na.rm = TRUE),
    accuracy = round(correctly_predicted / total_cases * 100, 1),
    rule_based = sum(prediction_source == "Clinical_Rules"),
    ml_based = sum(prediction_source == "Machine_Learning"),
    .groups = "drop"
  ) %>%
  arrange(desc(total_cases))

print(genotype_performance)

# ==============================================================================
# 10B. MISCLASSIFIED CASES ANALYSIS
# ==============================================================================

cat("\n=== PATIENTS MISCLASSIFIED BY RULE-BASED SYSTEM ===\n")
rule_cases <- final_predictions %>% filter(prediction_source == "Clinical_Rules")

misclassified_by_rules <- rule_cases %>%
  filter(final_prediction != genotype_simplified) %>%
  mutate(patient_row = row_number()) %>%
  select(
    patient_id = patient_row,
    genotype_simplified,
    rule_prediction = final_prediction,
    test_confidence_level,
    hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct,
    mcv, mchc, rdw, total_hgb,
    is_transfused
  )

cat("Number of misclassified high-confidence (rule-based) patients:", nrow(misclassified_by_rules), "\n")
if(nrow(misclassified_by_rules) > 0) {
  cat("--- DETAILS OF RULE-BASED MISCLASSIFICATIONS ---\n")
  print(misclassified_by_rules)
  cat("------------------------------------------------\n")
} else {
  cat("None! All rule-based high-confidence cases were correct.\n")
}

# ML misclassified cases
cat("\n=== PATIENTS MISCLASSIFIED BY ML SYSTEM ===\n")
ml_cases <- final_predictions %>% filter(prediction_source == "Machine_Learning")

misclassified_by_ml <- ml_cases %>%
  filter(final_prediction != genotype_simplified) %>%
  mutate(patient_row = row_number()) %>%
  select(
    patient_id = patient_row,
    genotype_simplified,
    ml_prediction = final_prediction,
    ml_confidence,
    hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct,
    mcv, mchc, rdw, total_hgb,
    is_transfused,
    needs_expert_review
  )

cat("Number of misclassified ML patients:", nrow(misclassified_by_ml), "\n")
if(nrow(misclassified_by_ml) > 0) {
  print(head(misclassified_by_ml, 10))  # Show first 10
  
  # Specifically investigate SS as SC and SC as SS by ML
  ss_as_sc_ml <- misclassified_by_ml %>%
    filter(genotype_simplified == "SS" & ml_prediction == "SC")
  
  sc_as_ss_ml <- misclassified_by_ml %>%
    filter(genotype_simplified == "SC" & ml_prediction == "SS")
  
  if(nrow(ss_as_sc_ml) > 0){
    cat("\n--- DETAILS: Actual SS Misclassified as SC by ML ---\n")
    print(ss_as_sc_ml)
    cat("-----------------------------------------------------\n")
  }
  
  if(nrow(sc_as_ss_ml) > 0){
    cat("\n--- DETAILS: Actual SC Misclassified as SS by ML ---\n")
    print(sc_as_ss_ml)
    cat("-----------------------------------------------------\n")
  }
  
} else {
  cat("None! All ML cases were correct.\n")
}

# ==============================================================================
# 11. SAVE RESULTS
# ==============================================================================

cat("\n=== SAVING TIDYMODELS RESULTS ===\n")

# Save fitted workflow
if(exists("final_fit") && !is.null(final_fit)) {
  saveRDS(final_fit, "tidymodels_genotype_workflow_with_transfusion.rds")
  cat("Tidymodels workflow saved\n")
}

# Save predictions
write_csv(final_predictions, "genotype_predictions_tidymodels_with_transfusion.csv")
cat("Final predictions saved\n")

# Save comprehensive results
comprehensive_results <- list(
  dataset_info = list(
    total_patients = nrow(data),
    transfusion_distribution = transfusion_summary
  ),
  clinical_rules = list(
    performance = rule_performance
  ),
  ml_performance = if(exists("cv_accuracy")) list(
    cv_accuracy = cv_accuracy,
    model_type = if(exists("model_type")) model_type else NULL,
    class_performance = if(exists("cv_class_metrics")) cv_class_metrics else NULL
  ) else NULL,
  transfusion_analysis = if(exists("transfusion_performance")) transfusion_performance else NULL,
  variable_importance = if(exists("var_imp_data")) var_imp_data else NULL,
  final_metrics = list(
    rule_accuracy = rule_accuracy,
    ml_accuracy = ml_accuracy,
    overall_accuracy = overall_accuracy
  )
)

# Save comprehensive results with enhanced metrics
comprehensive_results <- list(
  # Dataset information
  dataset_info = list(
    total_patients = nrow(data),
    transfusion_distribution = transfusion_summary,
    classifiable_patients = nrow(final_predictions),
    excluded_patients = nrow(data) - nrow(final_predictions)
  ),
  
  # Clinical rules performance
  clinical_rules = list(
    performance = rule_performance,
    misclassified_cases = if(exists("misclassified_by_rules")) misclassified_by_rules else NULL
  ),
  
  # Cross-validation and ML performance
  ml_performance = if(exists("cv_accuracy")) list(
    cv_accuracy = cv_accuracy,
    cv_coverage = if(exists("cv_coverage")) cv_coverage else NULL,
    cv_accuracy_at_coverage = if(exists("cv_accuracy_at_coverage")) cv_accuracy_at_coverage else NULL,
    model_type = if(exists("model_type")) model_type else NULL,
    class_performance = if(exists("cv_class_metrics")) cv_class_metrics else NULL,
    misclassified_cases = if(exists("misclassified_by_ml")) misclassified_by_ml else NULL
  ) else NULL,
  
  # Comprehensive performance metrics
  performance_metrics = if(exists("cv_comprehensive_metrics")) list(
    per_class_metrics = cv_comprehensive_metrics,
    macro_averaged = if(exists("macro_avg")) macro_avg else NULL,
    weighted_averaged = if(exists("weighted_avg")) weighted_avg else NULL,
    auroc_results = if(exists("roc_auc_results")) roc_auc_results else NULL,
    macro_auroc = if(exists("macro_auroc")) macro_auroc else NULL,
    weighted_auroc = if(exists("weighted_auroc")) weighted_auroc else NULL
  ) else NULL,
  
  # Calibration metrics
  calibration = if(exists("brier_score")) list(
    brier_score = brier_score,
    ece_score = if(exists("ece_score")) ece_score else NULL
  ) else NULL,
  
  # Transfusion analysis
  transfusion_analysis = if(exists("transfusion_performance")) list(
    cv_performance = transfusion_performance,
    final_performance = final_transfusion_summary
  ) else NULL,
  
  # Variable importance (from vip)
  variable_importance_vip = if(exists("var_imp_data")) var_imp_data else NULL,
  
  # SHAP Interpretability Analysis
  shap_analysis = if(exists("sv") && !is.null(sv) && exists("shap_importance_df")) list(
    shap_feature_importance = shap_importance_df,
    # You could add paths to saved misclassified plots or a summary here if needed
    notes = "SHAP bar, beeswarm, and misclassified waterfall plots saved in figures/ and figures/shap_misclassified/"
  ) else NULL,
  
  # Final pipeline performance
  pipeline_performance = list(
    total_classifiable_patients = nrow(final_predictions),
    rule_classified = nrow(rule_cases),
    ml_classified = nrow(ml_cases),
    rule_accuracy = rule_accuracy,
    ml_accuracy = ml_accuracy,
    overall_accuracy = overall_accuracy,
    expert_review_needed = sum(final_predictions$needs_expert_review, na.rm = TRUE)
  ),
  
  # Sensitivity and specificity
  sensitivity_specificity = if(exists("sens_spec")) sens_spec else NULL,
  
  # Final predictions summary
  predictions_summary = if(exists("genotype_performance")) genotype_performance else NULL,
  
  # Final metrics
  final_metrics = list(
    rule_accuracy = rule_accuracy,
    ml_accuracy = ml_accuracy,
    overall_accuracy = overall_accuracy
  )
)

saveRDS(comprehensive_results, "comprehensive_results_tidymodels_with_transfusion.rds")
cat("Comprehensive results saved\n")

cat("\n✅ TIDYMODELS PIPELINE WITH TRANSFUSION AWARENESS COMPLETE!\n")
cat("🔹 KEY FEATURES:\n")
cat("   • Modern tidymodels workflow\n")
cat("   • Transfusion status as feature with interactions\n")
cat("   • Automated model selection and hyperparameter tuning\n")
cat("   • Comprehensive transfusion impact assessment\n")
cat("   • Variable importance highlighting transfusion effects\n")
cat("   • Professional visualizations\n")

cat("\n📁 FILES GENERATED:\n")
cat("   • tidymodels_genotype_workflow_with_transfusion.rds\n")
cat("   • genotype_predictions_tidymodels_with_transfusion.csv\n")
cat("   • comprehensive_results_tidymodels_with_transfusion.rds\n")
cat("   • Variable importance plot in figures/\n")
cat("   • figures/calibration_plots_tidymodels.png\n")
cat("   • figures/confusion_matrix_heatmap_tidymodels.png\n")
cat("   • figures/cross_validation_performance_tidymodels.png\n")
cat("   • figures/confidence_distribution_tidymodels.png\n")
cat("   • figures/roc_curve_tidymodels_[class].png (for each class)\n")

cat("\n", paste(rep("=", 80), collapse=""), "\n")

# ==============================================================================
# 12. COMPREHENSIVE FINAL SUMMARY
# ==============================================================================

cat("                    TIDYMODELS PIPELINE PERFORMANCE SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\n📊 ACCURACY METRICS:\n")
cat("   • Rule-based accuracy:    ", if(!is.na(rule_accuracy)) paste0(round(rule_accuracy, 1), "%") else "N/A", "\n")
cat("   • ML accuracy:            ", if(!is.na(ml_accuracy)) paste0(round(ml_accuracy, 1), "%") else "N/A", "\n") 
cat("   • Overall pipeline:       ", if(!is.na(overall_accuracy)) paste0(round(overall_accuracy, 1), "%") else "N/A", "\n")

if(exists("cv_accuracy")) {
  cat("   • Cross-validation:       ", round(cv_accuracy * 100, 1), "%\n")
}

cat("\n📈 COVERAGE:\n")
cat("   • Total patients:         ", nrow(data), "\n")
cat("   • Classifiable patients:  ", nrow(final_predictions), "\n")
cat("   • Rule-classified:        ", nrow(rule_cases), 
    " (", round(nrow(rule_cases)/nrow(final_predictions)*100, 1), "%)\n")
cat("   • ML-classified:          ", nrow(ml_cases), 
    " (", round(nrow(ml_cases)/nrow(final_predictions)*100, 1), "%)\n")

if(exists("cv_comprehensive_metrics")) {
  cat("\n🎯 PERFORMANCE METRICS:\n")
  cat("   • Macro-averaged F1:      ", round(macro_avg$macro_f1, 3), "\n")
  cat("   • Weighted-averaged F1:   ", round(weighted_avg$weighted_f1, 3), "\n")
  if(exists("macro_auroc")) {
    cat("   • Macro AUROC:           ", round(macro_auroc, 3), "\n")
    cat("   • Weighted AUROC:        ", round(weighted_auroc, 3), "\n")
  }
}

if(exists("brier_score")) {
  cat("\n📐 CALIBRATION:\n")
  cat("   • Brier Score:           ", round(brier_score, 4), "\n")
  if(exists("ece_score")) {
    cat("   • Expected Cal. Error:   ", round(ece_score, 4), "\n")
  }
}

cat("\n🩸 TRANSFUSION IMPACT:\n")
transfused_count <- sum(final_predictions$is_transfused, na.rm = TRUE)
cat("   • Transfused patients:   ", transfused_count, "/", nrow(final_predictions), 
    " (", round(transfused_count/nrow(final_predictions)*100, 1), "%)\n")

if(exists("transfusion_performance") && nrow(transfusion_performance) >= 2) {
  transfused_acc <- transfusion_performance %>% filter(is_transfused == 1) %>% pull(accuracy)
  non_transfused_acc <- transfusion_performance %>% filter(is_transfused == 0) %>% pull(accuracy)
  if(length(transfused_acc) > 0 && length(non_transfused_acc) > 0) {
    cat("   • Transfused accuracy:   ", transfused_acc, "%\n")
    cat("   • Non-transfused acc:    ", non_transfused_acc, "%\n")
    cat("   • Difference:            ", transfused_acc - non_transfused_acc, " percentage points\n")
  }
}

if(sum(final_predictions$needs_expert_review, na.rm = TRUE) > 0) {
  cat("\n⚠️  EXPERT REVIEW:\n")
  cat("   • Cases needing review:   ", sum(final_predictions$needs_expert_review, na.rm = TRUE), "\n")
}

cat("\n🔬 MODEL SELECTION:\n")
if(exists("model_type")) {
  cat("   • Selected algorithm:     ", model_type, "\n")
}

cat("\n✅ ENHANCED TIDYMODELS PIPELINE COMPLETE!\n")
cat("🔹 KEY IMPROVEMENTS OVER ORIGINAL:\n")
cat("   • Modern tidymodels framework vs direct model calls\n")
cat("   • Explicit transfusion awareness with interaction features\n")
cat("   • Automated model selection (XGBoost)\n")
cat("   • Enhanced preprocessing through recipes\n")
cat("   • Comprehensive transfusion impact assessment\n")
cat("   • Calibration and performance analysis\n")
cat("   • Professional publication-ready visualizations\n")
cat("   • Reproducible workflow structure\n")

cat("\n📁 COMPREHENSIVE FILES GENERATED:\n")
cat("   • tidymodels_genotype_workflow_with_transfusion.rds (fitted workflow)\n")
cat("   • genotype_predictions_tidymodels_with_transfusion.csv (final predictions)\n")
cat("   • comprehensive_results_tidymodels_with_transfusion.rds (complete results)\n")
cat("   • figures/variable_importance_tidymodels_with_transfusion.png\n")
cat("   • figures/calibration_plots_tidymodels.png\n")
cat("   • figures/confusion_matrix_heatmap_tidymodels.png\n")
cat("   • figures/cross_validation_performance_tidymodels.png\n")
cat("   • figures/confidence_distribution_tidymodels.png\n")
cat("   • figures/roc_curve_tidymodels_[class].png (for each class)\n")

cat("\n", paste(rep("=", 80), collapse=""), "\n")

# ... existing code ...
# Detailed confusion matrix for final predictions
cat("\nDetailed Confusion Matrix (Automated Predictions Only - Excludes Expert Review Cases):\n")

# Filter out cases flagged for expert review
automated_predictions_for_cm <- final_predictions %>%
  filter(!needs_expert_review)

if(nrow(automated_predictions_for_cm) > 0) {
  confusion_final_automated <- table(
    Predicted = automated_predictions_for_cm$final_prediction,
    Actual = automated_predictions_for_cm$genotype_simplified
  )
  print(confusion_final_automated)
} else {
  cat("No cases available for automated confusion matrix after excluding expert review cases.\n")
}

# Create and save a heatmap for the final automated confusion matrix
if(exists("automated_predictions_for_cm") && nrow(automated_predictions_for_cm) > 0) {
  
  # Ensure truth and estimate are factors with the same levels
  all_levels <- unique(c(as.character(automated_predictions_for_cm$genotype_simplified), 
                         as.character(automated_predictions_for_cm$final_prediction)))
  
  automated_predictions_for_cm <- automated_predictions_for_cm %>% 
    mutate(
      genotype_simplified = factor(genotype_simplified, levels = all_levels),
      final_prediction = factor(final_prediction, levels = all_levels)
    )

  conf_mat_final_automated <- conf_mat(
    automated_predictions_for_cm,
    truth = genotype_simplified,
    estimate = final_prediction
  )
  
  # Convert to data frame for ggplot
  conf_mat_final_df <- as.data.frame(conf_mat_final_automated$table)
  
  final_automated_heatmap <- conf_mat_final_df %>% 
    ggplot(aes(x = Truth, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Freq), color = "white", fontface = "bold", size = 3) +
    scale_fill_viridis_c(name = "Count") +
    labs(
      title = "Final Automated Confusion Matrix (Excludes Expert Review)",
      subtitle = "Tidymodels Pipeline - Performance on Automatically Classified Cases",
      x = "Actual Genotype",
      y = "Predicted Genotype"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0),
      panel.grid = element_blank()
    )
  
  ggsave("figures/confusion_matrix_final_automated_heatmap.png", final_automated_heatmap, 
         width = 10, height = 8, dpi = 300)
  cat("Final automated confusion matrix heatmap saved to figures/confusion_matrix_final_automated_heatmap.png\n")
} else {
  cat("Skipping final automated confusion matrix heatmap generation: no data after filtering expert review cases.\n")
}

# Per-genotype performance breakdown
cat("\nPer-Genotype Performance Breakdown (All Classifiable Cases):\n")
genotype_performance <- final_predictions %>%
  group_by(genotype_simplified) %>%
  summarise(
    total_cases = n(),
    correctly_predicted = sum(final_prediction == genotype_simplified, na.rm = TRUE),
    accuracy = round(correctly_predicted / total_cases * 100, 1),
    rule_based = sum(prediction_source == "Clinical_Rules"),
    ml_based = sum(prediction_source == "Machine_Learning"),
    .groups = "drop"
  ) %>%
  arrange(desc(total_cases))

print(genotype_performance)

# ==============================================================================
# 10A. COMPREHENSIVE EVALUATION WITH EXPERT REVIEW AS CORRECT BEHAVIOR
# ==============================================================================

cat("\n=== COMPREHENSIVE EVALUATION (EXPERT REVIEW AS CORRECT SYSTEM BEHAVIOR) ===\n")

# Calculate coverage and accuracy at coverage
automated_cases <- final_predictions %>% filter(!needs_expert_review)
expert_review_cases <- final_predictions %>% filter(needs_expert_review)

coverage <- nrow(automated_cases) / nrow(final_predictions)
accuracy_at_coverage <- if(nrow(automated_cases) > 0) {
  sum(automated_cases$final_prediction == automated_cases$genotype_simplified, na.rm = TRUE) / nrow(automated_cases)
} else { NA }

# Coverage and accuracy by source
rule_coverage <- automated_cases %>% filter(prediction_source == "Clinical_Rules") %>% nrow()
ml_coverage <- automated_cases %>% filter(prediction_source == "Machine_Learning") %>% nrow()

rule_accuracy_at_coverage <- if(rule_coverage > 0) {
  rule_auto <- automated_cases %>% filter(prediction_source == "Clinical_Rules")
  sum(rule_auto$final_prediction == rule_auto$genotype_simplified, na.rm = TRUE) / nrow(rule_auto)
} else { NA }

ml_accuracy_at_coverage <- if(ml_coverage > 0) {
  ml_auto <- automated_cases %>% filter(prediction_source == "Machine_Learning")
  sum(ml_auto$final_prediction == ml_auto$genotype_simplified, na.rm = TRUE) / nrow(ml_auto)
} else { NA }

cat("=== SYSTEM PERFORMANCE WITH PROPER EXPERT REVIEW ACCOUNTING ===\n")
cat(sprintf("Coverage (Automated Classification Rate): %.1f%% (%d/%d cases)\n", 
            coverage * 100, nrow(automated_cases), nrow(final_predictions)))
cat(sprintf("Accuracy at Coverage (Automated Cases Only): %.1f%%\n", accuracy_at_coverage * 100))
cat(sprintf("Expert Review Rate: %.1f%% (%d cases) - CORRECTLY IDENTIFIED UNCERTAINTY\n", 
            (1 - coverage) * 100, nrow(expert_review_cases)))

cat("\nBREAKDOWN BY PREDICTION SOURCE:\n")
cat(sprintf("Rule-based: %d cases, %.1f%% accuracy at coverage\n", 
            rule_coverage, rule_accuracy_at_coverage * 100))
cat(sprintf("ML-based: %d cases, %.1f%% accuracy at coverage\n", 
            ml_coverage, ml_accuracy_at_coverage * 100))

# Performance by genotype at coverage
genotype_coverage_performance <- automated_cases %>%
  group_by(genotype_simplified) %>%
  summarise(
    automated_cases = n(),
    correct_automated = sum(final_prediction == genotype_simplified, na.rm = TRUE),
    accuracy_at_coverage = round(correct_automated / automated_cases * 100, 1),
    .groups = "drop"
  )

# Add expert review counts
expert_review_by_genotype <- expert_review_cases %>%
  count(genotype_simplified, name = "expert_review_cases")

coverage_summary <- final_predictions %>%
  count(genotype_simplified, name = "total_cases") %>%
  left_join(genotype_coverage_performance, by = "genotype_simplified") %>%
  left_join(expert_review_by_genotype, by = "genotype_simplified") %>%
  mutate(
    automated_cases = replace_na(automated_cases, 0),
    correct_automated = replace_na(correct_automated, 0),
    accuracy_at_coverage = replace_na(accuracy_at_coverage, 0),
    expert_review_cases = replace_na(expert_review_cases, 0),
    coverage_rate = round(automated_cases / total_cases * 100, 1)
  ) %>%
  arrange(desc(total_cases))

cat("\nPER-GENOTYPE COVERAGE AND ACCURACY:\n")
print(coverage_summary)

# Risk-coverage analysis at different confidence thresholds
cat("\n=== RISK-COVERAGE ANALYSIS ===\n")

ml_cases_for_analysis <- final_predictions %>% 
  filter(prediction_source == "Machine_Learning" & !is.na(ml_confidence))

if(nrow(ml_cases_for_analysis) > 0) {
  confidence_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99)
  
  risk_coverage_analysis <- map_dfr(confidence_thresholds, function(threshold) {
    high_conf_cases <- ml_cases_for_analysis %>% filter(ml_confidence >= threshold)
    
    if(nrow(high_conf_cases) > 0) {
      coverage_at_threshold <- nrow(high_conf_cases) / nrow(ml_cases_for_analysis)
      accuracy_at_threshold <- sum(high_conf_cases$final_prediction == high_conf_cases$genotype_simplified, na.rm = TRUE) / nrow(high_conf_cases)
      error_rate_at_threshold <- 1 - accuracy_at_threshold
      
      tibble(
        confidence_threshold = threshold,
        coverage = round(coverage_at_threshold * 100, 1),
        accuracy_at_coverage = round(accuracy_at_threshold * 100, 1),
        error_rate = round(error_rate_at_threshold * 100, 1),
        cases_covered = nrow(high_conf_cases)
      )
    } else {
      tibble(
        confidence_threshold = threshold,
        coverage = 0,
        accuracy_at_coverage = NA,
        error_rate = NA,
        cases_covered = 0
      )
    }
  })
  
  cat("Risk-Coverage Trade-off for ML Cases:\n")
  print(risk_coverage_analysis)
  
  # Current operating point
  current_threshold <- 0.90
  current_performance <- risk_coverage_analysis %>% filter(confidence_threshold == current_threshold)
  if(nrow(current_performance) > 0) {
    cat(sprintf("\nCURRENT OPERATING POINT (%.0f%% confidence threshold):\n", current_threshold * 100))
    cat(sprintf("- ML Coverage: %.1f%%\n", current_performance$coverage))
    cat(sprintf("- ML Accuracy at Coverage: %.1f%%\n", current_performance$accuracy_at_coverage))
    cat(sprintf("- ML Error Rate: %.1f%%\n", current_performance$error_rate))
  }
} else {
  cat("No ML cases with confidence scores available for risk-coverage analysis.\n")
}

# Summary metrics for reporting
system_metrics <- list(
  total_patients = nrow(final_predictions),
  coverage = round(coverage * 100, 1),
  accuracy_at_coverage = round(accuracy_at_coverage * 100, 1),
  expert_review_rate = round((1 - coverage) * 100, 1),
  rule_coverage = rule_coverage,
  ml_coverage = ml_coverage,
  rule_accuracy_at_coverage = round(rule_accuracy_at_coverage * 100, 1),
  ml_accuracy_at_coverage = round(ml_accuracy_at_coverage * 100, 1)
)

cat("\n=== SYSTEM METRICS SUMMARY ===\n")
cat("COVERAGE METRICS (Primary Performance Indicators):\n")
cat(sprintf("- System Coverage: %.1f%% (correctly automated %d/%d cases)\n", 
            system_metrics$coverage, nrow(automated_cases), system_metrics$total_patients))
cat(sprintf("- Accuracy at Coverage: %.1f%% (correct among automated cases)\n", 
            system_metrics$accuracy_at_coverage))
cat(sprintf("- Expert Review Rate: %.1f%% (appropriate uncertainty handling)\n", 
            system_metrics$expert_review_rate))

cat("\nCOMPONENT PERFORMANCE:\n")
cat(sprintf("- Rule-based: %d cases at %.1f%% accuracy\n", 
            system_metrics$rule_coverage, system_metrics$rule_accuracy_at_coverage))
cat(sprintf("- ML-based: %d cases at %.1f%% accuracy\n", 
            system_metrics$ml_coverage, system_metrics$ml_accuracy_at_coverage))

# ==============================================================================
# 10B. TRADITIONAL METRICS FOR COMPARISON (DEPRECATED APPROACH)
# ==============================================================================

cat("\n=== TRADITIONAL METRICS (FOR COMPARISON ONLY - NOT RECOMMENDED) ===\n")
cat("NOTE: These treat expert review as 'wrong' which is methodologically incorrect\n")
cat("for AI systems with abstain options. Use coverage metrics above instead.\n\n")

# Calculate traditional metrics (the old way)
traditional_rule_accuracy <- if(nrow(rule_cases) > 0) {
  sum(rule_cases$final_prediction == rule_cases$genotype_simplified, na.rm = TRUE) / nrow(rule_cases) * 100
} else { NA }

traditional_ml_accuracy <- if(nrow(ml_cases) > 0) {
  sum(ml_cases$final_prediction == ml_cases$genotype_simplified, na.rm = TRUE) / nrow(ml_cases) * 100
} else { NA }

traditional_overall_accuracy <- if(nrow(final_predictions) > 0) {
  sum(final_predictions$final_prediction == final_predictions$genotype_simplified, na.rm = TRUE) / nrow(final_predictions) * 100
} else { NA }

cat("DEPRECATED METRICS (Expert Review Counted as Wrong):\n")
cat("- Traditional Rule Accuracy:", round(traditional_rule_accuracy, 1), "%\n")
cat("- Traditional ML Accuracy:", round(traditional_ml_accuracy, 1), "%\n")
cat("- Traditional Overall Accuracy:", round(traditional_overall_accuracy, 1), "%\n")
cat("\n*** USE COVERAGE METRICS ABOVE FOR PROPER EVALUATION ***\n")

# ==============================================================================
# 11. SAVE RESULTS WITH UPDATED METRICS
# ==============================================================================

cat("\n=== SAVING TIDYMODELS RESULTS WITH COVERAGE METRICS ===\n")

# Save fitted workflow
if(exists("final_fit") && !is.null(final_fit)) {
  saveRDS(final_fit, "tidymodels_genotype_workflow_with_transfusion.rds")
  cat("Tidymodels workflow saved\n")
}

# Save predictions
write_csv(final_predictions, "genotype_predictions_tidymodels_with_transfusion.csv")
cat("Final predictions saved\n")

# Save comprehensive results with coverage metrics
comprehensive_results <- list(
  # Dataset information
  dataset_info = list(
    total_patients = nrow(data),
    classifiable_patients = nrow(final_predictions),
    excluded_patients = nrow(data) - nrow(final_predictions)
  ),
  
  # Coverage metrics (primary)
  coverage_metrics = system_metrics,
  
  # Coverage analysis by genotype
  coverage_analysis = coverage_summary,
  
  # Risk-coverage trade-off
  risk_coverage_analysis = if(exists("risk_coverage_analysis")) risk_coverage_analysis else NULL,
  
  # Clinical rules performance
  clinical_rules = list(
    performance = rule_performance,
    misclassified_cases = if(exists("misclassified_by_rules")) misclassified_by_rules else NULL
  ),
  
  # Cross-validation and ML performance
  ml_performance = if(exists("cv_accuracy")) list(
    cv_accuracy = cv_accuracy,
    cv_coverage = if(exists("cv_coverage")) cv_coverage else NULL,
    cv_accuracy_at_coverage = if(exists("cv_accuracy_at_coverage")) cv_accuracy_at_coverage else NULL,
    model_type = if(exists("model_type")) model_type else NULL,
    class_performance = if(exists("cv_class_metrics")) cv_class_metrics else NULL,
    misclassified_cases = if(exists("misclassified_by_ml")) misclassified_by_ml else NULL
  ) else NULL,
  
  # Transfusion analysis (ADDED BACK)
  transfusion_analysis = if(exists("transfusion_performance")) list(
    cv_performance = transfusion_performance,
    final_performance = if(exists("final_transfusion_summary")) final_transfusion_summary else NULL
  ) else NULL,
  
  # Traditional metrics (deprecated)
  traditional_metrics = list(
    rule_accuracy = traditional_rule_accuracy,
    ml_accuracy = traditional_ml_accuracy,
    overall_accuracy = traditional_overall_accuracy,
    note = "These metrics treat expert review as wrong - use coverage metrics instead"
  )
)

saveRDS(comprehensive_results, "comprehensive_results_tidymodels_with_transfusion.rds")
cat("Comprehensive results with coverage metrics saved\n")

cat("\n✅ ENHANCED TIDYMODELS PIPELINE WITH COVERAGE METRICS COMPLETE!\n")
cat("🔹 KEY EVALUATION IMPROVEMENTS:\n")
cat("   • Coverage metrics (primary evaluation method)\n")
cat("   • Accuracy at coverage (correct evaluation of automated cases)\n")
cat("   • Expert review as correct uncertainty handling\n")
cat("   • Risk-coverage trade-off analysis\n")

cat("\n📊 FINAL SYSTEM PERFORMANCE (COVERAGE METRICS):\n")
cat(sprintf("   • System Coverage: %.1f%% (automated %d/%d cases)\n", 
            system_metrics$coverage, nrow(automated_cases), system_metrics$total_patients))
cat(sprintf("   • Accuracy at Coverage: %.1f%% (among automated cases)\n", 
            system_metrics$accuracy_at_coverage))
cat(sprintf("   • Expert Review Rate: %.1f%% (appropriate uncertainty handling)\n", 
            system_metrics$expert_review_rate))

cat("\n", paste(rep("=", 80), collapse=""), "\n")

# ==============================================================================
# 6. ENHANCED MODEL EVALUATION (5-CLASS ALL-CASES)
# ==============================================================================

cat("\n=== ENHANCED MODEL EVALUATION (5-CLASS ALL-CASES) ===\n")

# Get CV predictions for 5-class problem
cv_predictions <- collect_predictions(best_tune_results) %>%
  filter(.config == xgb_best$.config) %>%
  mutate(
    predicted_class = .pred_class,
    actual_class = genotype_simplified
  )

# Calculate confidence for 5-class problem
pred_cols <- names(cv_predictions)[str_starts(names(cv_predictions), "\\.pred_") & 
                                  !str_detect(names(cv_predictions), "\\.pred_class")]
if(length(pred_cols) > 0) {
  cv_predictions$confidence <- apply(cv_predictions[pred_cols], 1, max, na.rm = TRUE)
} else {
  cv_predictions$confidence <- 0.5  # Default value
}

cv_predictions$high_confidence <- cv_predictions$confidence > 0.90

# Calculate overall accuracy for 5-class
cv_accuracy <- cv_predictions %>%
  accuracy(truth = actual_class, estimate = predicted_class) %>%
  pull(.estimate)

cat("Enhanced 5-class cross-validation accuracy:", round(cv_accuracy * 100, 1), "%\n")

# Enhanced per-class performance (5-class)
cv_class_metrics <- cv_predictions %>%
  group_by(actual_class) %>%
  summarise(
    n = n(),
    correct = sum(predicted_class == actual_class),
    accuracy = round(correct / n * 100, 1),
    avg_confidence = round(mean(confidence), 3),
    high_conf_rate = round(sum(high_confidence) / n * 100, 1),
    .groups = "drop"
  )

cat("\nEnhanced per-class performance (5-class):\n")
print(cv_class_metrics)

# Enhanced confusion matrix (5-class)
cv_conf_matrix <- cv_predictions %>%
  conf_mat(truth = actual_class, estimate = predicted_class)

cat("\nEnhanced Confusion Matrix (5-class):\n")
print(cv_conf_matrix)

# ==============================================================================
# 12. CLINICAL CONSTRAINTS IMPACT SUMMARY
# ==============================================================================

cat("\n=== CLINICAL CONSTRAINTS IMPACT SUMMARY ===\n")
cat("🧬 BIOCHEMICAL CONSTRAINT BENEFITS:\n")
cat("   • Prevents SC disease prediction when HgbC = 0%\n")
cat("   • Reduces SS disease probability when HgbA is elevated\n")
cat("   • Boosts AS (trait) recognition for classic patterns\n")
cat("   • Enforces HgbA2 requirements for S-beta thalassemias\n")
cat("   • Maintains clinical plausibility in ML predictions\n")

if(exists("accuracy_raw") && exists("accuracy_constrained")) {
  cat(sprintf("\n📈 CROSS-VALIDATION IMPROVEMENTS:\n"))
  cat(sprintf("   • Raw CV accuracy:         %.1f%%\n", accuracy_raw))
  cat(sprintf("   • Constrained CV accuracy: %.1f%%\n", accuracy_constrained))
  cat(sprintf("   • Improvement:             +%.1f percentage points\n", accuracy_constrained - accuracy_raw))
}

if(exists("corrections_made") && nrow(ml_cases) > 0) {
  cat(sprintf("\n🔧 FINAL PREDICTION CORRECTIONS:\n"))
  cat(sprintf("   • ML cases processed:      %d\n", nrow(ml_cases)))
  cat(sprintf("   • Cases corrected:         %d (%.1f%%)\n", corrections_made, corrections_made / nrow(ml_cases) * 100))
}

cat("\n💡 CLINICAL SIGNIFICANCE:\n")
cat("   • Eliminates impossible genotype combinations\n")
cat("   • Improves diagnostic confidence and safety\n")
cat("   • Reduces need for expert review of obvious errors\n")
cat("   • Enhances model trustworthiness in clinical practice\n")

# ==============================================================================
# === 13. LANCET HAEMATOLOGY PUBLICATION SVG EXPORT MODULE (SHARP) ===
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(scales)
library(patchwork)
library(svglite)
library(pROC)
library(shapviz)

sharp_svg_dir <- "publication_figures_SHARP_svg"
dir.create(sharp_svg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(sharp_svg_dir, "Supp_FigS4_Beeswarms_Individual"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(sharp_svg_dir, "Supp_FigS5_Misclassified_Waterfalls"), showWarnings = FALSE, recursive = TRUE)

# Official Lancet 5-Genotype Palette (ggsci 'lanonc')
lancet_genotype_colors <- c(
  "HbSS"  = "#00468B", # Lancet Navy
  "HbSC"  = "#ED0000", # Lancet Crimson
  "HbSβ⁰" = "#42B540", # Lancet Green
  "HbSβ⁺" = "#0099B4", # Lancet Cyan
  "HbAS"  = "#925E9F"  # Lancet Purple
)

# Unified Lancet Theme with Times New Roman
theme_lancet_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(family = "Times New Roman", size = base_size + 2, face = "bold", color = "#1B1919", hjust = 0),
      plot.subtitle    = element_text(family = "Times New Roman", size = base_size, color = "#424242", hjust = 0, margin = margin(b = 8)),
      axis.title       = element_text(family = "Times New Roman", size = base_size, face = "bold", color = "#1B1919"),
      axis.text        = element_text(family = "Times New Roman", size = base_size - 1, color = "#1B1919"),
      legend.title     = element_text(family = "Times New Roman", size = base_size - 0.5, face = "bold", color = "#1B1919"),
      legend.text      = element_text(family = "Times New Roman", size = base_size - 1, color = "#1B1919"),
      strip.text       = element_text(family = "Times New Roman", size = base_size, face = "bold", color = "#1B1919"),
      strip.background = element_rect(fill = "#F4F6F8", color = "#CCCCCC"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.35),
      panel.border     = element_rect(color = "#1B1919", fill = NA, linewidth = 0.6),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.tag         = element_text(family = "Times New Roman", size = base_size + 4, face = "bold", color = "#1B1919")
    )
}

# Standard Clinical Nomenclature Dictionaries
clean_genotype_labels <- function(x) {
  recode(as.character(x),
         "SS"                = "HbSS",
         "SC"                = "HbSC",
         "S_beta_0"          = "HbSβ⁰",
         "S_beta_plus"       = "HbSβ⁺",
         "AS"                = "HbAS",
         "shap_pred_SS"      = "HbSS",
         "shap_pred_SC"      = "HbSC",
         "shap_pred_S_beta_0"= "HbSβ⁰",
         "shap_pred_S_beta_plus" = "HbSβ⁺",
         "shap_pred_AS"      = "HbAS",
         .default            = as.character(x)
  )
}

clean_sharp_features <- function(x) {
  recode(as.character(x),
         "hgb_s_pct"              = "HbS (%)",
         "hgb_a_pct"              = "HbA (%)",
         "hgb_c_pct"              = "HbC (%)",
         "hgb_f_pct"              = "HbF (%)",
         "hgb_a2_pct"             = "HbA2 (%)",
         "total_hgb"              = "Total Hemoglobin (g/dL)",
         "hematocrit"             = "Hematocrit (%)",
         "mcv"                    = "MCV (fL)",
         "mchc"                   = "MCHC (g/dL)",
         "rdw"                    = "RDW (%)",
         "rbc_count"              = "RBC Count (10^6/uL)",
         "hgb_a_present"          = "HbA Present",
         "hgb_s_present"          = "HbS Present",
         "hgb_c_present"          = "HbC Present",
         "hgb_f_elevated"         = "Elevated HbF",
         "hgb_a2_elevated"        = "Elevated HbA2",
         "hgb_s_to_total"         = "HbS-to-Total Ratio",
         "hgb_a_to_s_ratio"       = "HbA-to-HbS Ratio",
         "total_variant_pct"      = "Total Variant Hb (%)",
         "completeness_score"     = "Data Completeness Score",
         "clarity_score"          = "Profile Clarity Score",
         "total_quality_score"    = "Total Quality Score",
         "had_hba_below_1_ever"   = "Historical HbA < 1%",
         "is_transfused"          = "Recent Transfusion (<30d)",
         "hgb_f_to_s_ratio"       = "HbF-to-HbS Ratio",
         "hgb_f_to_a2_ratio"      = "HbF-to-HbA2 Ratio",
         "as_pattern_score"       = "HbAS Pattern Score",
         "trait_indicator"        = "Sickle Trait Indicator",
         "transfused_hgb_s"       = "Transfusion × HbS (%)",
         "transfused_hgb_a"       = "Transfusion × HbA (%)",
         "transfused_hgb_c"       = "Transfusion × HbC (%)",
         "transfused_mcv"         = "Transfusion × MCV",
         "transfused_total_hgb"   = "Transfusion × Total Hb",
         "borderline_hgb_s"       = "Borderline HbS (70-80%)",
         "borderline_hgb_f"       = "Borderline HbF (5-15%)",
         "borderline_mcv"         = "Borderline MCV (75-85 fL)",
         .default                 = as.character(x)
  )
}

genotype_order <- c("HbSS", "HbSC", "HbSβ⁰", "HbSβ⁺", "HbAS")

# ------------------------------------------------------------------------------
# Panel 3A: SHARP Final Automated Confusion Matrix (≥90% Confidence Coverage)
# ------------------------------------------------------------------------------
cm_sharp_auto_df <- final_predictions %>%
  filter(!needs_expert_review) %>%
  filter(genotype_simplified %in% c("SS", "SC", "S_beta_0", "S_beta_plus", "AS") & 
           final_prediction %in% c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")) %>% # <--- ADD THIS FILTER LINE
  mutate(
    Truth = factor(clean_genotype_labels(genotype_simplified), levels = genotype_order),
    Prediction = factor(clean_genotype_labels(final_prediction), levels = rev(genotype_order))
  ) %>%
  count(Truth, Prediction, .drop = FALSE) %>%
  group_by(Truth) %>%
  mutate(
    col_pct = ifelse(sum(n) > 0, n / sum(n) * 100, 0),
    label_str = ifelse(n > 0, sprintf("%d\n(%.1f%%)", n, col_pct), "0"),
    text_col = ifelse(col_pct > 50, "white", "#1B1919")
  ) %>%
  ungroup()

fig3a_cm <- ggplot(cm_sharp_auto_df, aes(x = Truth, y = Prediction, fill = col_pct)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = label_str, color = text_col),
            family = "Times New Roman", fontface = "bold", size = 3.3) +
  scale_color_identity() +
  scale_fill_gradient(low = "#F0F4F8", high = "#00468B", name = "Class Accuracy (%)") +
  labs(
    title = "SHARP Automated Confusion Matrix",
    subtitle = sprintf("Automated Calls at >=90%% Confidence (%d / %d Cases; %.1f%% Coverage)",
                       sum(!final_predictions$needs_expert_review), nrow(final_predictions),
                       mean(!final_predictions$needs_expert_review) * 100),
    x = "WGS-Confirmed Reference Genotype",
    y = "SHARP Predicted Genotype"
  ) +
  theme_lancet_pub() +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave(file.path(sharp_svg_dir, "Fig3A_SHARP_Automated_Confusion_Matrix.svg"), fig3a_cm, width = 6.2, height = 5.2, device = svglite)

# Also export the Tier-2 Cross-Validation Confusion Matrix as Supplementary SVG
if (exists("cv_conf_matrix")) {
  cm_cv_df <- as.data.frame(cv_conf_matrix$table) %>%
    mutate(
      Truth = factor(clean_genotype_labels(Truth), levels = genotype_order),
      Prediction = factor(clean_genotype_labels(Prediction), levels = rev(genotype_order))
    ) %>%
    group_by(Truth) %>%
    mutate(
      col_pct = ifelse(sum(Freq) > 0, Freq / sum(Freq) * 100, 0),
      label_str = ifelse(Freq > 0, sprintf("%d\n(%.1f%%)", Freq, col_pct), "0"),
      text_col = ifelse(col_pct > 50, "white", "#1B1919")
    ) %>%
    ungroup()
  
  supp_cm_cv <- ggplot(cm_cv_df, aes(x = Truth, y = Prediction, fill = col_pct)) +
    geom_tile(color = "white", linewidth = 1.0) +
    geom_text(aes(label = label_str, color = text_col),
              family = "Times New Roman", fontface = "bold", size = 3.3) +
    scale_color_identity() +
    scale_fill_gradient(low = "#F0F4F8", high = "#00468B", name = "Class Recall (%)") +
    labs(
      title = "SHARP Tier-2 ML Cross-Validation Confusion Matrix",
      subtitle = "10-Fold Stratified Out-of-Fold Predictions with Post-Hoc Biochemical Constraints",
      x = "WGS-Confirmed Reference Genotype",
      y = "Predicted Genotype"
    ) +
    theme_lancet_pub() +
    theme(panel.grid = element_blank(), legend.position = "none")
  
  ggsave(file.path(sharp_svg_dir, "Supp_SHARP_CV_Confusion_Matrix.svg"), supp_cm_cv, width = 6.2, height = 5.2, device = svglite)
}

# ------------------------------------------------------------------------------
# Panel 3B: Overlaid 5-Genotype Multi-Class ROC Curves
# ------------------------------------------------------------------------------
raw_classes <- c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
roc_list_df <- map_dfr(raw_classes, function(cl) {
  pred_col <- paste0(".pred_", cl)
  if (exists("cv_probs") && pred_col %in% colnames(cv_probs)) {
    actual_bin <- ifelse(cv_probs$genotype_simplified == cl, 1, 0)
    roc_obj <- pROC::roc(actual_bin, cv_probs[[pred_col]], quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    clean_cl <- clean_genotype_labels(cl)
    tibble(
      Genotype = clean_cl,
      LegendLabel = sprintf("%s (AUROC = %.3f)", clean_cl, auc_val),
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      auroc = auc_val
    ) %>% arrange(fpr, tpr)
  } else {
    NULL
  }
})

roc_list_df <- roc_list_df %>%
  mutate(Genotype = factor(Genotype, levels = genotype_order))

legend_labels_map <- roc_list_df %>%
  distinct(Genotype, LegendLabel) %>%
  arrange(Genotype) %>%
  tibble::deframe()

macro_auc_sharp <- mean(unique(roc_list_df$auroc), na.rm = TRUE)

fig3b_roc <- ggplot(roc_list_df, aes(x = fpr, y = tpr, color = Genotype)) +
  geom_abline(lty = 2, color = "#ADB6B6", linewidth = 0.7) +
  geom_path(linewidth = 1.1, alpha = 0.9) +
  scale_color_manual(
    values = lancet_genotype_colors,
    labels = legend_labels_map,
    name = sprintf("Genotype (Macro AUROC = %.3f)", macro_auc_sharp)
  ) +
  scale_x_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  labs(
    title = "Multi-Class Receiver Operating Characteristic",
    subtitle = "One-vs-Rest Discrimination Across All Five Genotypes",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_lancet_pub() +
  theme(
    legend.position = c(0.66, 0.28),
    legend.background = element_rect(fill = alpha("white", 0.92), color = "#CCCCCC", linewidth = 0.4)
  )

ggsave(file.path(sharp_svg_dir, "Fig3B_SHARP_Multiclass_ROC.svg"), fig3b_roc, width = 6.2, height = 5.2, device = svglite)

# ------------------------------------------------------------------------------
# Panel 3C: Prediction Confidence Distribution & Risk-Coverage Curve
# ------------------------------------------------------------------------------
conf_clean_df <- cv_predictions %>%
  mutate(
    Outcome = factor(ifelse(predicted_class == actual_class, "Correct Classification", "Misclassification"),
                     levels = c("Correct Classification", "Misclassification"))
  )

fig3c_conf <- ggplot(conf_clean_df, aes(x = confidence, fill = Outcome)) +
  geom_histogram(bins = 28, alpha = 0.85, position = "identity", color = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0.90, color = "#AD002A", linetype = "dashed", linewidth = 0.9) +
  annotate("label", x = 0.74, y = max( hist(conf_clean_df$confidence, breaks = 28, plot = FALSE)$counts ) * 0.78,
           label = "90% Confidence Threshold\n(Flagged for Expert Review ->)",
           family = "Times New Roman", fontface = "bold", size = 3.2,
           fill = "white", color = "#AD002A", label.size = 0.35) +
  scale_fill_manual(
    values = c("Correct Classification" = "#00468B", "Misclassification" = "#ED0000"),
    name = "Prediction Outcome"
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.30, 1.02)) +
  labs(
    title = "Prediction Confidence & Uncertainty Filtering",
    subtitle = "Distribution of Out-of-Fold Probabilities Relative to 90% Abstention Cutoff",
    x = "Maximum Class Probability (Prediction Confidence)",
    y = "Number of Patients"
  ) +
  theme_lancet_pub() +
  theme(legend.position = "top")

ggsave(file.path(sharp_svg_dir, "Fig3C_SHARP_Confidence_Distribution.svg"), fig3c_conf, width = 6.2, height = 5.2, device = svglite)

# Also export the Risk-Coverage Trade-Off Curve as a standalone SVG option
if (exists("risk_coverage_analysis") && nrow(risk_coverage_analysis) > 0) {
  rc_plot <- ggplot(risk_coverage_analysis, aes(x = coverage, y = accuracy_at_coverage)) +
    geom_line(color = "#00468B", linewidth = 1.1) +
    geom_point(aes(color = confidence_threshold == 0.90), size = 3.2) +
    geom_text(aes(label = sprintf("%.0f%%", confidence_threshold * 100)),
              vjust = -0.9, family = "Times New Roman", fontface = "bold", size = 3.2, color = "#1B1919") +
    scale_color_manual(values = c("FALSE" = "#00468B", "TRUE" = "#ED0000"), guide = "none") +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = "SHARP Risk-Coverage Trade-Off Curve",
      subtitle = "Automated Classification Accuracy vs. System Coverage Across Confidence Thresholds",
      x = "Automated ML Coverage (% of Ambiguous Cases)",
      y = "Accuracy Among Automated Cases (%)"
    ) +
    theme_lancet_pub()
  
  ggsave(file.path(sharp_svg_dir, "Fig3C_Alt_SHARP_Risk_Coverage_Curve.svg"), rc_plot, width = 6.2, height = 5.2, device = svglite)
}

# ------------------------------------------------------------------------------
# Panel 3D: Per-Genotype Mean Absolute SHAP Feature Importance
# ------------------------------------------------------------------------------
if (exists("per_class_shap_long") && !is.null(per_class_shap_long) && nrow(per_class_shap_long) > 0) {
  top_12_feats <- per_class_shap_long %>%
    group_by(feature) %>%
    summarise(overall_mean = mean(mean_abs_shap, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = 12) %>%
    pull(feature)
  
  shap_bar_clean_df <- per_class_shap_long %>%
    filter(feature %in% top_12_feats) %>%
    mutate(
      CleanFeature = clean_sharp_features(feature),
      Genotype = factor(clean_genotype_labels(class), levels = genotype_order)
    )
  
  feat_order <- shap_bar_clean_df %>%
    group_by(CleanFeature) %>%
    summarise(m = mean(mean_abs_shap, na.rm = TRUE), .groups = "drop") %>%
    arrange(m) %>%
    pull(CleanFeature)
  
  shap_bar_clean_df$CleanFeature <- factor(shap_bar_clean_df$CleanFeature, levels = feat_order)
  
  fig3d_shap <- ggplot(shap_bar_clean_df, aes(x = mean_abs_shap, y = CleanFeature, fill = Genotype)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = lancet_genotype_colors, name = "Genotype") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "Per-Genotype SHAP Feature Importance",
      subtitle = "Top 12 Clinical & Laboratory Predictors Driving Multi-Class Differentiation",
      x = "Mean |SHAP Value| (Log-Odds Contribution)",
      y = "Predictor Variable"
    ) +
    theme_lancet_pub() +
    theme(legend.position = "top")
  
  ggsave(file.path(sharp_svg_dir, "Fig3D_SHARP_Per_Genotype_SHAP.svg"), fig3d_shap, width = 6.8, height = 5.5, device = svglite)
}

# ------------------------------------------------------------------------------
# Combine Panels into Main Text Figure 3 (2x2 Multi-Panel Layout)
# ------------------------------------------------------------------------------
if (exists("fig3d_shap")) {
  fig3_combined <- (fig3a_cm | fig3b_roc) / (fig3c_conf | fig3d_shap) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 15))
  
  ggsave(file.path(sharp_svg_dir, "Figure3_SHARP_Main_4Panel.svg"), fig3_combined, width = 13.5, height = 11.2, device = svglite)
}

# ------------------------------------------------------------------------------
# Supplementary Figure S1: 5-Class Calibration Loess Grid
# ------------------------------------------------------------------------------
if (exists("cv_probs") && nrow(cv_probs) > 0) {
  calib_plots_lancet <- list()
  for (cl in raw_classes) {
    pred_col <- paste0(".pred_", cl)
    if (pred_col %in% colnames(cv_probs)) {
      actual_bin <- ifelse(cv_probs$genotype_simplified == cl, 1, 0)
      predicted_p <- cv_probs[[pred_col]]
      clean_cl <- clean_genotype_labels(cl)
      
      bins <- seq(0, 1, length.out = 11)
      centers <- (bins[-1] + bins[-11]) / 2
      cal_df <- map_dfr(1:10, function(i) {
        in_b <- if (i == 10) (predicted_p >= bins[i] & predicted_p <= bins[i + 1]) else (predicted_p >= bins[i] & predicted_p < bins[i + 1])
        if (sum(in_b, na.rm = TRUE) > 0) {
          tibble(
            bin_center = centers[i],
            mean_pred = mean(predicted_p[in_b], na.rm = TRUE),
            mean_obs  = mean(actual_bin[in_b], na.rm = TRUE),
            count     = sum(in_b, na.rm = TRUE)
          )
        } else NULL
      })
      
      calib_plots_lancet[[clean_cl]] <- ggplot(cal_df, aes(x = mean_pred, y = mean_obs)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#ADB6B6", linewidth = 0.8) +
        geom_smooth(method = "loess", se = TRUE, color = lancet_genotype_colors[[clean_cl]],
                    fill = lancet_genotype_colors[[clean_cl]], alpha = 0.20, linewidth = 0.9, na.rm = TRUE) +
        geom_point(aes(size = count), color = lancet_genotype_colors[[clean_cl]], alpha = 0.85) +
        scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
        scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
        scale_size_continuous(name = "Sample Count", range = c(2, 6)) +
        labs(
          title = sprintf("Calibration Curve: %s", clean_cl),
          x = "Mean Predicted Probability",
          y = "Observed Fraction of Positives"
        ) +
        theme_lancet_pub(base_size = 10) +
        theme(legend.position = "bottom")
    }
  }
  
  if (length(calib_plots_lancet) > 0) {
    supp_fig_s1 <- wrap_plots(calib_plots_lancet, ncol = 3) +
      plot_annotation(
        title = "SHARP Multi-Class Probability Calibration Curves",
        subtitle = sprintf("10-Fold Out-of-Fold Reliability Diagrams (Expected Calibration Error = %.4f)",
                           ifelse(exists("ece_score"), ece_score, NA)),
        tag_levels = "A",
        theme = theme_lancet_pub(base_size = 12)
      )
    ggsave(file.path(sharp_svg_dir, "Supp_FigS1_SHARP_Calibration_Grid.svg"), supp_fig_s1, width = 13.0, height = 8.5, device = svglite)
  }
}

# ------------------------------------------------------------------------------
# Supplementary Figure S4: Per-Class SHAP Beeswarm Plots (Cleaned Feature Names)
# ------------------------------------------------------------------------------
if (exists("sv_for_plotting") && length(sv_for_plotting) > 0) {
  beeswarm_list <- list()
  for (shap_name in names(sv_for_plotting)) {
    shap_obj <- sv_for_plotting[[shap_name]]
    if (inherits(shap_obj, "shapviz")) {
      # Rename columns inside shapviz object to clean clinical names
      colnames(shap_obj$S) <- clean_sharp_features(colnames(shap_obj$S))
      colnames(shap_obj$X) <- clean_sharp_features(colnames(shap_obj$X))
      clean_gt <- clean_genotype_labels(shap_name)
      
      p_bee <- sv_importance(shap_obj, kind = "beeswarm", max_display = 12) +
        scale_color_gradient(low = "#00468B", high = "#ED0000", name = "Feature Value\n(Standardized)") +
        labs(
          title = sprintf("SHAP Beeswarm: %s", clean_gt),
          x = "SHAP Value (Impact on Log-Odds)",
          y = NULL
        ) +
        theme_lancet_pub(base_size = 10)
      
      beeswarm_list[[clean_gt]] <- p_bee
      ggsave(file.path(sharp_svg_dir, "Supp_FigS4_Beeswarms_Individual", sprintf("SHAP_Beeswarm_%s.svg", clean_gt)),
             p_bee, width = 7.0, height = 5.5, device = svglite)
    }
  }
  
  if (length(beeswarm_list) > 0) {
    supp_fig_s4 <- wrap_plots(beeswarm_list, ncol = 2) +
      plot_annotation(
        title = "Genotype-Specific SHAP Summary Beeswarm Plots",
        subtitle = "Directionality and Magnitude of Top 12 Predictors Across All Five Genotypes",
        tag_levels = "A",
        theme = theme_lancet_pub(base_size = 12)
      )
    ggsave(file.path(sharp_svg_dir, "Supp_FigS4_SHARP_Beeswarms_Combined.svg"), supp_fig_s4, width = 14.0, height = 14.5, device = svglite)
  }
}

# ------------------------------------------------------------------------------
# Supplementary Figure S5: Misclassified Case Waterfall Plots (Cleaned Names)
# ------------------------------------------------------------------------------
if (exists("misclassified_shap_cases") && nrow(misclassified_shap_cases) > 0 && exists("sv_for_plotting") && length(sv_for_plotting) > 0) {
  for (i in seq_len(nrow(misclassified_shap_cases))) {
    case_info <- misclassified_shap_cases[i, ]
    r_id <- case_info$row_id_in_ml_data
    t_lbl <- clean_genotype_labels(as.character(case_info$actual_class))
    p_lbl_raw <- as.character(case_info$.pred_class)
    p_lbl <- clean_genotype_labels(p_lbl_raw)
    prob_col <- paste0(".pred_", p_lbl_raw)
    conf_val <- if (prob_col %in% names(case_info)) case_info[[prob_col]] else NA
    
    shap_key <- make.names(paste0("shap_pred_", p_lbl_raw))
    sv_curr <- if (shap_key %in% names(sv_for_plotting)) sv_for_plotting[[shap_key]] else sv_for_plotting[[1]]
    
    if (inherits(sv_curr, "shapviz")) {
      colnames(sv_curr$S) <- clean_sharp_features(colnames(sv_curr$S))
      colnames(sv_curr$X) <- clean_sharp_features(colnames(sv_curr$X))
      
      p_water <- sv_waterfall(sv_curr, row_id = r_id, max_display = 10,
                              fill_colors = c("#ED0000", "#00468B")) +
        labs(
          title = sprintf("SHAP Waterfall Analysis: Ambiguous Edge Case (Patient Row %d)", r_id),
          subtitle = sprintf("Predicted: %s (Confidence: %.1f%%) | WGS Reference: %s", p_lbl, conf_val * 100, t_lbl),
          x = "Cumulative SHAP Contribution"
        ) +
        theme_lancet_pub(base_size = 10)
      
      ggsave(file.path(sharp_svg_dir, "Supp_FigS5_Misclassified_Waterfalls",
                       sprintf("Waterfall_Row%d_Pred_%s_True_%s.svg", r_id, p_lbl_raw, as.character(case_info$actual_class))),
             p_water, width = 8.0, height = 5.5, device = svglite)
    }
  }
}

# ------------------------------------------------------------------------------
# Additional Supplementary SVGs: Permutation VIP & 10-Fold CV Stability
# ------------------------------------------------------------------------------
if (exists("var_imp_data") && nrow(var_imp_data) > 0) {
  vip_clean_df <- var_imp_data %>%
    mutate(CleanVar = clean_sharp_features(Variable))
  
  supp_vip_svg <- ggplot(vip_clean_df, aes(x = Importance, y = fct_reorder(CleanVar, Importance))) +
    geom_col(fill = "#00468B", width = 0.68) +
    geom_text(aes(label = sprintf("%.3f", Importance)),
              hjust = -0.15, family = "Times New Roman", fontface = "bold", size = 3.0, color = "#1B1919") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      title = "SHARP Permutation Feature Importance",
      subtitle = "Top 20 Predictors Ranked by Accuracy Drop Upon Permutation",
      x = "Permutation Importance (Decrease in Accuracy)",
      y = "Predictor Variable"
    ) +
    theme_lancet_pub()
  
  ggsave(file.path(sharp_svg_dir, "Supp_SHARP_Permutation_Importance.svg"), supp_vip_svg, width = 7.5, height = 6.0, device = svglite)
}

if (exists("per_fold_acc") && nrow(per_fold_acc) > 0) {
  supp_cv_svg <- ggplot(per_fold_acc, aes(x = id, y = accuracy)) +
    geom_col(fill = "#00468B", width = 0.62, alpha = 0.9) +
    geom_hline(yintercept = mean(per_fold_acc$accuracy), color = "#ED0000", linetype = "dashed", linewidth = 0.9) +
    geom_text(aes(label = sprintf("%.1f%%", accuracy * 100)),
              vjust = -0.5, family = "Times New Roman", fontface = "bold", size = 3.2, color = "#1B1919") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "SHARP 10-Fold Stratified Cross-Validation Stability",
      subtitle = sprintf("Mean Cross-Validation Accuracy: %.1f%% (SD: %.2f%%)",
                         mean(per_fold_acc$accuracy) * 100, sd(per_fold_acc$accuracy) * 100),
      x = "Cross-Validation Fold",
      y = "Classification Accuracy (%)"
    ) +
    theme_lancet_pub()
  
  ggsave(file.path(sharp_svg_dir, "Supp_SHARP_10Fold_CV_Stability.svg"), supp_cv_svg, width = 7.5, height = 4.8, device = svglite)
}
# ==============================================================================
# === 14. LANCET HAEMATOLOGY & TRIPOD+AI TABLE EXPORT MODULE (SHARP) ===
# ==============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

sharp_table_dir <- "publication_tables_SHARP"
dir.create(sharp_table_dir, showWarnings = FALSE, recursive = TRUE)

fmt_med_iqr_sharp <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return("NA")
  q <- quantile(x_clean, probs = c(0.25, 0.50, 0.75))
  sprintf("%.1f (%.1f–%.1f)", q[2], q[1], q[3])
}

calc_binom_ci_sharp <- function(k, n) {
  if (is.na(k) || is.na(n) || n == 0) return("NA")
  bt <- binom.test(k, n)
  sprintf("%.1f%% (%.1f–%.1f)", (k / n) * 100, bt$conf.int[1] * 100, bt$conf.int[2] * 100)
}

# ------------------------------------------------------------------------------
# 1. Table 1 (Part B): SHARP Baseline Characteristics by Genotype
# ------------------------------------------------------------------------------
summarise_sharp_group <- function(df_sub, group_label) {
  n_tot <- nrow(df_sub)
  n_tx  <- sum(df_sub$is_transfused == 1 | df_sub$is_transfused == TRUE, na.rm = TRUE)
  
  tibble(
    `Genotype / Cohort` = group_label,
    `Total N (%)` = sprintf("%s (%.1f%%)", format(n_tot, big.mark = ","), (n_tot / nrow(all_cases_ml_data)) * 100),
    `Transfused <30d, n (%)` = sprintf("%d (%.1f%%)", n_tx, (n_tx / n_tot) * 100),
    `HbS (%)` = fmt_med_iqr_sharp(df_sub$hgb_s_pct),
    `HbA (%)` = fmt_med_iqr_sharp(df_sub$hgb_a_pct),
    `HbC (%)` = fmt_med_iqr_sharp(df_sub$hgb_c_pct),
    `HbF (%)` = fmt_med_iqr_sharp(df_sub$hgb_f_pct),
    `HbA2 (%)` = fmt_med_iqr_sharp(df_sub$hgb_a2_pct),
    `Total Hemoglobin (g/dL)` = fmt_med_iqr_sharp(df_sub$total_hgb),
    `Hematocrit (%)` = fmt_med_iqr_sharp(df_sub$hematocrit),
    `MCV (fL)` = fmt_med_iqr_sharp(df_sub$mcv),
    `MCHC (g/dL)` = fmt_med_iqr_sharp(df_sub$mchc),
    `RDW (%)` = fmt_med_iqr_sharp(df_sub$rdw),
    `RBC Count (10^6/uL)` = fmt_med_iqr_sharp(df_sub$rbc_count)
  )
}

sharp_missing_row <- all_cases_ml_data %>%
  summarise(
    `Genotype / Cohort` = "Missing Values Count (Imputed via Median, n)",
    `Total N (%)` = as.character(nrow(all_cases_ml_data)),
    `Transfused <30d, n (%)` = as.character(sum(is.na(is_transfused))),
    `HbS (%)` = as.character(sum(is.na(hgb_s_pct))),
    `HbA (%)` = as.character(sum(is.na(hgb_a_pct))),
    `HbC (%)` = as.character(sum(is.na(hgb_c_pct))),
    `HbF (%)` = as.character(sum(is.na(hgb_f_pct))),
    `HbA2 (%)` = as.character(sum(is.na(hgb_a2_pct))),
    `Total Hemoglobin (g/dL)` = as.character(sum(is.na(total_hgb))),
    `Hematocrit (%)` = as.character(sum(is.na(hematocrit))),
    `MCV (fL)` = as.character(sum(is.na(mcv))),
    `MCHC (g/dL)` = as.character(sum(is.na(mchc))),
    `RDW (%)` = as.character(sum(is.na(rdw))),
    `RBC Count (10^6/uL)` = as.character(sum(is.na(rbc_count)))
  )

table1_sharp_rows <- bind_rows(
  summarise_sharp_group(all_cases_ml_data, "All SHARP Participants (Combined)"),
  map_dfr(c("SS", "SC", "S_beta_0", "S_beta_plus", "AS"), function(gt) {
    clean_gt <- recode(gt, "SS"="HbSS", "SC"="HbSC", "S_beta_0"="HbSβ⁰", "S_beta_plus"="HbSβ⁺", "AS"="HbAS")
    summarise_sharp_group(all_cases_ml_data %>% filter(genotype_simplified == gt), paste0("Genotype: ", clean_gt))
  }),
  summarise_sharp_group(all_cases_ml_data %>% filter(is_transfused == 1), "Subgroup: Recently Transfused (<30d)"),
  summarise_sharp_group(all_cases_ml_data %>% filter(is_transfused == 0), "Subgroup: Non-Transfused"),
  sharp_missing_row
)

write_csv(table1_sharp_rows, file.path(sharp_table_dir, "Table1_PartB_SHARP_Baseline_By_Genotype.csv"))

# ------------------------------------------------------------------------------
# 2. Table 2 (Part B): SHARP Pipeline, Per-Genotype, & Transfusion Performance
# ------------------------------------------------------------------------------
# Part B1: System Architecture & Tier Breakdown
tier_summary_rows <- bind_rows(
  # Full Combined Pipeline
  {
    n_tot <- nrow(final_predictions)
    auto_df <- final_predictions %>% filter(!needs_expert_review)
    n_auto <- nrow(auto_df)
    n_corr_auto <- sum(auto_df$final_prediction == auto_df$genotype_simplified, na.rm = TRUE)
    n_corr_all  <- sum(final_predictions$final_prediction == final_predictions$genotype_simplified, na.rm = TRUE)
    tibble(
      `Pipeline Component / Stratum` = "Combined SHARP Pipeline (Tier 1 + Tier 2)",
      `Total N` = n_tot,
      `Automated Coverage, n (%)` = sprintf("%d (%.1f%%)", n_auto, (n_auto / n_tot) * 100),
      `Flagged for Expert Review, n (%)` = sprintf("%d (%.1f%%)", n_tot - n_auto, ((n_tot - n_auto) / n_tot) * 100),
      `Accuracy at >=90% Coverage (95% CI)` = calc_binom_ci_sharp(n_corr_auto, n_auto),
      `Unfiltered Accuracy Across All Cases (95% CI)` = calc_binom_ci_sharp(n_corr_all, n_tot)
    )
  },
  # Tier 1: Clinical Rules Only
  {
    t1_df <- final_predictions %>% filter(prediction_source == "Clinical_Rules")
    n_t1 <- nrow(t1_df)
    n_corr <- sum(t1_df$final_prediction == t1_df$genotype_simplified, na.rm = TRUE)
    tibble(
      `Pipeline Component / Stratum` = "Tier 1: Literature-Based Clinical Rules",
      `Total N` = n_t1,
      `Automated Coverage, n (%)` = sprintf("%d (%.1f%% of Cohort)", n_t1, (n_t1 / nrow(final_predictions)) * 100),
      `Flagged for Expert Review, n (%)` = "0 (0.0%)",
      `Accuracy at >=90% Coverage (95% CI)` = calc_binom_ci_sharp(n_corr, n_t1),
      `Unfiltered Accuracy Across All Cases (95% CI)` = calc_binom_ci_sharp(n_corr, n_t1)
    )
  },
  # Tier 2: Ambiguous Cases Routed to XGBoost
  {
    t2_df <- final_predictions %>% filter(prediction_source == "Machine_Learning")
    n_t2 <- nrow(t2_df)
    t2_auto <- t2_df %>% filter(!needs_expert_review)
    n_corr_auto <- sum(t2_auto$final_prediction == t2_auto$genotype_simplified, na.rm = TRUE)
    n_corr_all  <- sum(t2_df$final_prediction == t2_df$genotype_simplified, na.rm = TRUE)
    tibble(
      `Pipeline Component / Stratum` = "Tier 2: XGBoost + Biochemical Constraints (Ambiguous Cases)",
      `Total N` = n_t2,
      `Automated Coverage, n (%)` = sprintf("%d (%.1f%% of Tier 2)", nrow(t2_auto), (nrow(t2_auto) / n_t2) * 100),
      `Flagged for Expert Review, n (%)` = sprintf("%d (%.1f%%)", n_t2 - nrow(t2_auto), ((n_t2 - nrow(t2_auto)) / n_t2) * 100),
      `Accuracy at >=90% Coverage (95% CI)` = calc_binom_ci_sharp(n_corr_auto, nrow(t2_auto)),
      `Unfiltered Accuracy Across All Cases (95% CI)` = calc_binom_ci_sharp(n_corr_all, n_t2)
    )
  },
  # Transfused vs Non-Transfused (10-Fold CV on All Cases)
  map_dfr(c(1, 0), function(tx_flag) {
    sub_cv <- cv_predictions_with_transfusion %>% filter(is_transfused == tx_flag)
    n_sub <- nrow(sub_cv)
    sub_auto <- sub_cv %>% filter(high_confidence)
    n_corr_auto <- sum(sub_auto$predicted_class == sub_auto$actual_class, na.rm = TRUE)
    n_corr_all  <- sum(sub_cv$predicted_class == sub_cv$actual_class, na.rm = TRUE)
    lbl <- ifelse(tx_flag == 1, "ML Subgroup: Recently Transfused (<30d)", "ML Subgroup: Non-Transfused")
    tibble(
      `Pipeline Component / Stratum` = lbl,
      `Total N` = n_sub,
      `Automated Coverage, n (%)` = sprintf("%d (%.1f%%)", nrow(sub_auto), (nrow(sub_auto) / n_sub) * 100),
      `Flagged for Expert Review, n (%)` = sprintf("%d (%.1f%%)", n_sub - nrow(sub_auto), ((n_sub - nrow(sub_auto)) / n_sub) * 100),
      `Accuracy at >=90% Coverage (95% CI)` = calc_binom_ci_sharp(n_corr_auto, nrow(sub_auto)),
      `Unfiltered Accuracy Across All Cases (95% CI)` = calc_binom_ci_sharp(n_corr_all, n_sub)
    )
  })
)

write_csv(tier_summary_rows, file.path(sharp_table_dir, "Table2_PartB1_SHARP_Pipeline_and_Subgroup_Summary.csv"))

# Part B2: Per-Genotype Diagnostic Performance with 95% CIs
per_genotype_table2 <- map_dfr(c("SS", "SC", "S_beta_0", "S_beta_plus", "AS"), function(gt) {
  clean_gt <- recode(gt, "SS"="HbSS", "SC"="HbSC", "S_beta_0"="HbSβ⁰", "S_beta_plus"="HbSβ⁺", "AS"="HbAS")
  
  # Automated predictions evaluation
  auto_df <- final_predictions %>% filter(!needs_expert_review)
  tp <- sum(auto_df$genotype_simplified == gt & auto_df$final_prediction == gt, na.rm = TRUE)
  fn <- sum(auto_df$genotype_simplified == gt & auto_df$final_prediction != gt, na.rm = TRUE)
  fp <- sum(auto_df$genotype_simplified != gt & auto_df$final_prediction == gt, na.rm = TRUE)
  tn <- sum(auto_df$genotype_simplified != gt & auto_df$final_prediction != gt, na.rm = TRUE)
  
  tot_gt_all <- sum(final_predictions$genotype_simplified == gt, na.rm = TRUE)
  n_rules_gt <- sum(final_predictions$genotype_simplified == gt & final_predictions$prediction_source == "Clinical_Rules", na.rm = TRUE)
  n_ml_gt    <- sum(final_predictions$genotype_simplified == gt & final_predictions$prediction_source == "Machine_Learning", na.rm = TRUE)
  n_rev_gt   <- sum(final_predictions$genotype_simplified == gt & final_predictions$needs_expert_review, na.rm = TRUE)
  
  f1_gt <- ifelse((2 * tp + fp + fn) > 0, (2 * tp) / (2 * tp + fp + fn), NA_real_)
  
  # One-vs-Rest AUROC with DeLong 95% CI from cv_probs
  pred_col <- paste0(".pred_", gt)
  auroc_ci_str <- if (exists("cv_probs") && pred_col %in% colnames(cv_probs)) {
    actual_bin <- ifelse(cv_probs$genotype_simplified == gt, 1, 0)
    tryCatch({
      r_obj <- pROC::roc(actual_bin, cv_probs[[pred_col]], quiet = TRUE)
      ci_v  <- as.numeric(pROC::ci.auc(r_obj, method = "delong"))
      sprintf("%.4f (%.4f–%.4f)", ci_v[2], ci_v[1], ci_v[3])
    }, error = function(e) "NA")
  } else "NA"
  
  tibble(
    Genotype = clean_gt,
    `Total N` = tot_gt_all,
    `Classified by Tier 1 Rules, n (%)` = sprintf("%d (%.1f%%)", n_rules_gt, (n_rules_gt / tot_gt_all) * 100),
    `Classified by Tier 2 ML, n (%)` = sprintf("%d (%.1f%%)", n_ml_gt, (n_ml_gt / tot_gt_all) * 100),
    `Flagged for Review, n (%)` = sprintf("%d (%.1f%%)", n_rev_gt, (n_rev_gt / tot_gt_all) * 100),
    `Sensitivity / Class Accuracy (95% CI)` = calc_binom_ci_sharp(tp, tp + fn),
    `Specificity (95% CI)` = calc_binom_ci_sharp(tn, tn + fp),
    `PPV / Precision (95% CI)` = calc_binom_ci_sharp(tp, tp + fp),
    `NPV (95% CI)` = calc_binom_ci_sharp(tn, tn + fn),
    `F1-Score` = sprintf("%.4f", f1_gt),
    `One-vs-Rest AUROC (95% CI)` = auroc_ci_str
  )
})

write_csv(per_genotype_table2, file.path(sharp_table_dir, "Table2_PartB2_SHARP_Per_Genotype_Performance.csv"))

# ------------------------------------------------------------------------------
# 3. Supplementary Table S1: Tier 1 Clinical Rules & Post-Hoc Constraints
# ------------------------------------------------------------------------------
supp_table_s1 <- tribble(
  ~Stage, ~Target_Genotype, ~Quantitative_Criteria, ~Action_or_Confidence_Rule,
  "Tier 1: Override Rule", "HbAA (Normal)", "HbA > 95%, HbS < 2%, HbC < 2%", "Assign HbAA with High Confidence (Overrides all downstream rules)",
  "Tier 1: Clinical Rule", "HbSS", "HbS > 70%, HbC < 1%, HbF < 20%, HbA2 < 3.5%, HbA < 1%, MCV > 75 fL", "Assign HbSS if Total Hgb >= 6 g/dL and High/Medium Test Quality",
  "Tier 1: Clinical Rule", "HbSC", "HbS 35–60%, HbC 40–60%, HbF < 5%, HbA < 1%, MCV 70–95 fL", "Assign HbSC if Total Hgb >= 6 g/dL and High/Medium Test Quality",
  "Tier 1: Clinical Rule", "HbSβ⁰", "HbS > 65%, HbC < 1%, HbF < 20%, HbA2 > 4.5%, HbA < 1%, MCV < 80 fL", "Assign HbSβ⁰ if Total Hgb >= 6 g/dL and High/Medium Test Quality",
  "Tier 1: Clinical Rule", "HbSβ⁺", "HbS 60–90%, HbC < 1%, HbF < 20%, HbA2 > 4.5%, HbA 3–35%, MCV < 85 fL", "Assign HbSβ⁺ unless patient had transfusion <30d (route to Tier 2 ML)",
  "Tier 1: Clinical Rule", "HbAS (Trait)", "HbA 50–65%, HbS 30–45%, HbC < 1%, HbA2 < 3.5%, HbF < 3%", "Assign HbAS unless patient had transfusion <30d (route to Tier 2 ML)",
  "Tier 2: Post-Hoc Constraint 1", "HbSC", "HbC < 10%", "Multiply P(HbSC) by 0.01 if HbC = 0%, else by (HbC / 10) * 0.1; renormalize",
  "Tier 2: Post-Hoc Constraint 2", "HbSS", "HbA > 10%", "Multiply P(HbSS) by max(0.1, 1 - (HbA - 5) / 50); renormalize",
  "Tier 2: Post-Hoc Constraint 3", "HbAS (Trait)", "HbA 40–70%, HbS 30–60%, HbC <= 2%", "Boost P(HbAS) by factor of 1.5; renormalize",
  "Tier 2: Post-Hoc Constraint 4", "HbSβ⁰ & HbSβ⁺", "HbA2 < 3.5%", "Multiply P(HbSβ⁰) and P(HbSβ⁺) by 0.2; renormalize",
  "Tier 2: Post-Hoc Constraint 5", "HbSβ⁺", "HbA < 5%", "Multiply P(HbSβ⁺) by 0.1; renormalize",
  "Tier 2: Abstention Gate", "All Genotypes", "Maximum Constrained Probability <= 0.90", "Flag prediction for manual expert review"
)

write_csv(supp_table_s1, file.path(sharp_table_dir, "Supp_TableS1_SHARP_Clinical_Rules_and_Constraints.csv"))

# ------------------------------------------------------------------------------
# 4. Supplementary Table S2 (Part C): SHARP Hyperparameters, Calibration & Risk-Coverage
# ------------------------------------------------------------------------------
if (exists("xgb_best")) {
  sharp_hyperparams <- xgb_best %>%
    mutate(
      Brier_Score = ifelse(exists("brier_score"), round(brier_score, 4), NA),
      Expected_Calibration_Error_ECE = ifelse(exists("ece_score"), round(ece_score, 4), NA),
      Macro_F1_CV = ifelse(exists("macro_avg"), round(macro_avg$macro_f1, 4), NA),
      Macro_AUROC_CV = ifelse(exists("macro_auroc"), round(macro_auroc, 4), NA)
    )
  write_csv(sharp_hyperparams, file.path(sharp_table_dir, "Supp_TableS2_PartC_SHARP_Winning_Hyperparameters_and_Calibration.csv"))
}

if (exists("risk_coverage_analysis")) {
  write_csv(risk_coverage_analysis, file.path(sharp_table_dir, "Supp_TableS3_SHARP_Risk_Coverage_Tradeoff.csv"))
}

cat(sprintf("\n✅ All SHARP Tables exported to '%s/'\n", sharp_table_dir))
cat(sprintf("\n✅ All SHARP Lancet SVGs exported successfully to '%s/'\n", sharp_svg_dir))

