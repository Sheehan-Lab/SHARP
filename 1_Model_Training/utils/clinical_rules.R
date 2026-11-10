# In utils/clinical_rules.R

# ==============================================================================
# Clinical Rules Module
# ==============================================================================
# Implementation of clinical decision rules for hemoglobin genotype prediction
# VERSION WITH AA OVERRIDE
# ==============================================================================

library(dplyr)

#' Apply enhanced clinical rules for genotype prediction (5-class with AA override)
#' @param data A single-row list or dataframe containing patient data
#' @return A list with prediction results
apply_clinical_rules <- function(data) {
  
  # Convert all inputs to numeric, providing defaults for optional CBC/quality values
  hgb_s_pct <- as.numeric(data$hgb_s_pct)
  hgb_a_pct <- as.numeric(data$hgb_a_pct)
  hgb_c_pct <- as.numeric(data$hgb_c_pct)
  hgb_f_pct <- as.numeric(data$hgb_f_pct)
  hgb_a2_pct <- as.numeric(data$hgb_a2_pct)
  is_transfused <- as.numeric(data$is_transfused)
  mcv <- if(!is.null(data$mcv) && !is.na(data$mcv)) as.numeric(data$mcv) else NA
  total_hgb <- if(!is.null(data$total_hgb) && !is.na(data$total_hgb)) as.numeric(data$total_hgb) else NA
  
  # --- AA OVERRIDE: Check this condition first! ---
  # If a patient has dominant HbA and negligible variants, it's AA.
  if (!is.na(hgb_a_pct) && hgb_a_pct > 95 &&
      !is.na(hgb_s_pct) && hgb_s_pct < 2 &&
      !is.na(hgb_c_pct) && hgb_c_pct < 2) {
    return(list(
      prediction = "AA",
      high_confidence = TRUE,
      source = "Clinical Rule"
    ))
  }
  
  # --- If not AA, proceed with other rules ---
  # These rules are based on your original logic, simplified for clarity
  
  # Rule: SS (Sickle Cell Disease)
  if (hgb_s_pct > 75 && hgb_a_pct == 0 && hgb_c_pct == 0) {
    if (is.na(mcv) || mcv > 80) {
      return(list(prediction = "SS", high_confidence = TRUE, source = "Clinical Rule"))
    }
  }
  
  # Rule: SC Disease
  if (hgb_s_pct >= 45 && hgb_s_pct <= 55 && hgb_c_pct >= 40) {
    if (is.na(mcv) || (mcv >= 70 && mcv <= 90)) {
      return(list(prediction = "SC", high_confidence = TRUE, source = "Clinical Rule"))
    }
  }
  
  # Rule: S-beta0 Thalassemia
  if (hgb_s_pct > 70 && hgb_a_pct == 0 && hgb_a2_pct > 5.7) {
    if (is.na(mcv) || mcv < 80) {
      return(list(prediction = "S_beta_0", high_confidence = TRUE, source = "Clinical Rule"))
    }
  }
  
  # Rule: S-beta+ Thalassemia (non-transfused only)
  if (is_transfused == 0 && hgb_s_pct >= 60 && hgb_a_pct >= 5 && hgb_a_pct <= 30 && hgb_a2_pct > 5.7) {
    if (is.na(mcv) || mcv < 80) {
      return(list(prediction = "S_beta_plus", high_confidence = TRUE, source = "Clinical Rule"))
    }
  }
  
  # Rule: AS Trait (non-transfused only)
  if (is_transfused == 0 && hgb_a_pct >= 50 && hgb_s_pct >= 35) {
    return(list(prediction = "AS", high_confidence = TRUE, source = "Clinical Rule"))
  }
  
  # If no high-confidence rule was met, return AMBIGUOUS
  return(list(
    prediction = "AMBIGUOUS",
    high_confidence = FALSE,
    source = "Needs ML"
  ))
}