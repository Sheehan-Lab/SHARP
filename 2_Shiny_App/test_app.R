# Test Script for Hemoglobin Genotype Prediction App
# This script demonstrates the core functionality without the Shiny interface

# Load required libraries
library(tidyverse)
library(tidymodels)

# Source the utility functions
source("utils/data_validation.R")
source("utils/clinical_rules.R")
source("utils/shap_generator.R")
source("utils/helpers.R")

# Load the model
cat("Loading trained model...\n")
model <- readRDS("data/tidymodels_genotype_workflow_with_transfusion.rds")
cat("Model loaded successfully!\n\n")

# Test Case 1: Classic SS pattern
cat("=== Test Case 1: Classic SS Pattern ===\n")
test_case_1 <- list(
  hgb_s_pct = 85.2,
  hgb_a_pct = 0.0,
  hgb_c_pct = 0.0,
  hgb_f_pct = 12.8,
  hgb_a2_pct = 2.0,
  is_transfused = 0,
  total_hgb = 8.5,
  mcv = 95
)

# Test clinical rules
rule_result_1 <- apply_clinical_rules(test_case_1)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_1$prediction, "\n")
cat("Confidence:", round(rule_result_1$confidence, 3), "\n")
cat("High Confidence:", rule_result_1$high_confidence, "\n")
cat("Rationale:", rule_result_1$rationale, "\n\n")

# Test Case 2: SC pattern
cat("=== Test Case 2: SC Pattern ===\n")
test_case_2 <- list(
  hgb_s_pct = 49.1,
  hgb_a_pct = 0.0,
  hgb_c_pct = 47.9,
  hgb_f_pct = 1.5,
  hgb_a2_pct = 1.5,
  is_transfused = 0,
  total_hgb = 11.2,
  mcv = 78
)

rule_result_2 <- apply_clinical_rules(test_case_2)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_2$prediction, "\n")
cat("Confidence:", round(rule_result_2$confidence, 3), "\n")
cat("High Confidence:", rule_result_2$high_confidence, "\n")
cat("Rationale:", rule_result_2$rationale, "\n\n")

# Test Case 3: Ambiguous case requiring ML
cat("=== Test Case 3: Ambiguous Case (Requires ML) ===\n")
test_case_3 <- list(
  hgb_s_pct = 67.3,
  hgb_a_pct = 15.2,
  hgb_c_pct = 0.0,
  hgb_f_pct = 8.1,
  hgb_a2_pct = 9.4,
  is_transfused = 0,
  total_hgb = 9.8,
  mcv = 72
)

rule_result_3 <- apply_clinical_rules(test_case_3)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_3$prediction, "\n")
cat("High Confidence:", rule_result_3$high_confidence, "\n")

if (!rule_result_3$high_confidence) {
  cat("Applying machine learning...\n")
  ml_result_3 <- predict_with_ml(test_case_3, model)
  cat("ML Prediction:", ml_result_3$prediction, "\n")
  cat("ML Confidence:", round(ml_result_3$confidence, 3), "\n")
  
  # Test SHAP explanation
  cat("Generating SHAP explanation...\n")
  shap_result_3 <- safely_generate_shap(test_case_3, model)
  if (!is.null(shap_result_3$result)) {
    cat("SHAP analysis completed successfully!\n")
    cat("Top 3 contributing features:\n")
    top_features <- shap_result_3$result$feature_contributions[1:3, ]
    for (i in 1:nrow(top_features)) {
      cat(paste0(i, ". ", top_features$feature[i], " (SHAP: ", 
                round(top_features$shap_value[i], 3), ")\n"))
    }
  } else {
    cat("SHAP analysis failed:", shap_result_3$error, "\n")
  }
}
cat("\n")

# Test Case 4: AS Pattern (Sickle Cell Trait) - NEW 5-CLASS TEST
cat("=== Test Case 4: AS Pattern (Sickle Cell Trait) ===\n")
test_case_4 <- list(
  hgb_s_pct = 38.5,
  hgb_a_pct = 58.2,
  hgb_c_pct = 0.0,
  hgb_f_pct = 1.1,
  hgb_a2_pct = 2.2,
  is_transfused = 0,
  total_hgb = 13.5,
  mcv = 87
)

rule_result_4 <- apply_clinical_rules(test_case_4)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_4$prediction, "\n")
cat("Confidence:", round(rule_result_4$confidence, 3), "\n")
cat("High Confidence:", rule_result_4$high_confidence, "\n")
cat("Rationale:", rule_result_4$rationale, "\n")

if (!rule_result_4$high_confidence) {
  cat("Applying machine learning for AS classification...\n")
  ml_result_4 <- predict_with_ml(test_case_4, model)
  cat("ML Prediction:", ml_result_4$prediction, "\n")
  cat("ML Confidence:", round(ml_result_4$confidence, 3), "\n")
}
cat("\n")

# Test Case 5: Transfused patient (affects interpretation)
cat("=== Test Case 5: Transfused Patient ===\n")
test_case_5 <- list(
  hgb_s_pct = 42.1,
  hgb_a_pct = 35.8,
  hgb_c_pct = 0.0,
  hgb_f_pct = 1.2,
  hgb_a2_pct = 2.9,
  is_transfused = 1,  # Recently transfused
  total_hgb = 10.5,
  mcv = 88
)

rule_result_5 <- apply_clinical_rules(test_case_5)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_5$prediction, "\n")
cat("High Confidence:", rule_result_5$high_confidence, "\n")
cat("Rationale:", rule_result_5$rationale, "\n")

if (!rule_result_5$high_confidence) {
  cat("Applying machine learning (transfusion-aware)...\n")
  ml_result_5 <- predict_with_ml(test_case_5, model)
  cat("ML Prediction:", ml_result_5$prediction, "\n")
  cat("ML Confidence:", round(ml_result_5$confidence, 3), "\n")
}

# Test input validation
cat("=== Testing Input Validation ===\n")
invalid_input <- list(
  hgb_s_pct = 150,  # Invalid: > 100
  hgb_a_pct = -5,   # Invalid: < 0
  hgb_c_pct = 0,
  hgb_f_pct = NULL, # Invalid: missing
  hgb_a2_pct = 2.1,
  is_transfused = NULL  # Invalid: missing
)

# Create a mock input object for validation
mock_input <- as.list(invalid_input)
validation_result <- validate_inputs(mock_input)
cat("Validation Result:\n")
cat("Valid:", validation_result$valid, "\n")
if (!validation_result$valid) {
  cat("Errors:\n")
  for (error in validation_result$errors) {
    cat("-", error, "\n")
  }
}

# Test Case 6: User's Problem Case - Transfused AS pattern with HgbC = 0 but ML suggesting SC
cat("=== Test Case 6: Clinical Constraint Test - TRANSFUSED (HgbS=55%, HgbA=45%, HgbC=0%) ===\n")
test_case_6 <- list(
  hgb_s_pct = 55,
  hgb_a_pct = 45,
  hgb_c_pct = 0,  # This should eliminate SC as a possibility
  hgb_f_pct = 0,
  hgb_a2_pct = 2,
  is_transfused = 1,  # TRANSFUSED - this is key!
  total_hgb = 13.0,
  mcv = 85
)

rule_result_6 <- apply_clinical_rules(test_case_6)
cat("Clinical Rule Result:\n")
cat("Prediction:", rule_result_6$prediction, "\n")
cat("High Confidence:", rule_result_6$high_confidence, "\n")

if (!rule_result_6$high_confidence) {
  cat("Applying machine learning with clinical constraints...\n")
  ml_result_6 <- predict_with_ml(test_case_6, model)
  cat("ML Prediction:", ml_result_6$prediction, "\n")
  cat("ML Confidence:", round(ml_result_6$confidence, 3), "\n")
  
  # Show the probability breakdown
  cat("Constrained Probability Breakdown:\n")
  prob_df <- ml_result_6$probabilities %>%
    pivot_longer(everything(), names_to = "class", values_to = "prob") %>%
    mutate(
      class = str_remove(class, "^\\.pred_"),
      percentage = round(prob * 100, 1)
    ) %>%
    arrange(desc(prob))
  
  for(i in 1:nrow(prob_df)) {
    cat(paste0("  ", prob_df$class[i], ": ", prob_df$percentage[i], "%\n"))
  }
  
  cat("Note: SC should now be heavily penalized since HgbC = 0%\n")
}
cat("\n")

cat("\n=== All Tests Completed ===\n")
cat("The app is ready to use! Run launch_app.R to start the Shiny interface.\n") 