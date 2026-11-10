# In utils/clinical_rules.R

# ==============================================================================
# Clinical Rules Module - COMPLETE VERSION
# ==============================================================================
# Contains the main rules function with AA override AND the required helper functions.
# ==============================================================================

library(dplyr)

#' Apply enhanced clinical rules for genotype prediction (with AA override)
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
  
  # --- AA OVERRIDE: Check this condition first! ---
  if (!is.na(hgb_a_pct) && hgb_a_pct > 95 &&
      !is.na(hgb_s_pct) && hgb_s_pct < 2 &&
      !is.na(hgb_c_pct) && hgb_c_pct < 2) {
    return(list(
      prediction = "AA (Not a Sickle Cell Patient)",
      high_confidence = TRUE,
      source = "clinical_rules",
      rule_applied = "AA_Override",
      confidence = 1.0,
      rationale = "Dominant HgbA (>95%) with negligible HgbS and HgbC is definitive for a normal (AA) genotype."
    ))
  }
  
  # --- If not AA, proceed with other rules ---
  # Rule: SS (Sickle Cell Disease)
  if (hgb_s_pct > 75 && hgb_a_pct == 0 && hgb_c_pct == 0) {
    if (is.na(mcv) || mcv > 80) {
      return(list(prediction = "SS", high_confidence = TRUE, confidence = 0.95, rationale = "Classic SS pattern.", rule_applied = "SS_classic"))
    }
  }
  
  # Rule: SC Disease
  if (hgb_s_pct >= 45 && hgb_s_pct <= 55 && hgb_c_pct >= 40) {
    if (is.na(mcv) || (mcv >= 70 && mcv <= 90)) {
      return(list(prediction = "SC", high_confidence = TRUE, confidence = 0.90, rationale = "Classic SC pattern.", rule_applied = "SC_classic"))
    }
  }
  
  # Rule: S-beta0 Thalassemia
  if (hgb_s_pct > 70 && hgb_a_pct == 0 && hgb_a2_pct > 5.7) {
    if (is.na(mcv) || mcv < 80) {
      return(list(prediction = "S_beta_0", high_confidence = TRUE, confidence = 0.85, rationale = "Classic S-beta-0 pattern.", rule_applied = "S_beta_0_classic"))
    }
  }
  
  # Rule: S-beta+ Thalassemia (non-transfused only)
  if (is_transfused == 0 && hgb_s_pct >= 60 && hgb_a_pct >= 5 && hgb_a_pct <= 30 && hgb_a2_pct > 5.7) {
    if (is.na(mcv) || mcv < 80) {
      return(list(prediction = "S_beta_plus", high_confidence = TRUE, confidence = 0.80, rationale = "Classic S-beta-plus pattern.", rule_applied = "S_beta_plus_classic"))
    }
  }
  
  # Rule: AS Trait (non-transfused only)
  if (is_transfused == 0 && hgb_a_pct >= 50 && hgb_s_pct >= 35) {
    return(list(prediction = "AS", high_confidence = TRUE, confidence = 0.90, rationale = "Classic AS (trait) pattern.", rule_applied = "AS_classic"))
  }
  
  # If no high-confidence rule was met, return AMBIGUOUS
  return(list(
    prediction = "AMBIGUOUS",
    high_confidence = FALSE,
    source = "clinical_rules", # Source is still rules, but confidence is low
    rule_applied = NULL,
    confidence = 0,
    rationale = "No high-confidence clinical rule was met. Proceeding to machine learning analysis."
  ))
}


## HELPER FUNCTION - RESTORED FOR THE APP'S UI ##
#' Get detailed rule explanations for display
#' @param rule_name Character, name of the rule applied
#' @return Character explanation
get_rule_explanation <- function(rule_name) {
  # Helper function for null coalescing
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  explanations <- list(
    "AA_Override" = "The override rule for normal (AA) genotype was met. This pattern is not indicative of sickle cell disease or trait.",
    "SS_classic" = "Classic SS pattern: Predominantly HgbS with no HbA or HbC present. Elevated HgbF is common. Usually associated with severe disease.",
    "SC_classic" = "Classic SC pattern: Approximately equal amounts of HgbS and HgbC with no HbA. Generally milder than SS disease but can have complications.",
    "S_beta_0_classic" = "S-β⁰ thalassemia pattern: High HgbS with elevated HgbA2 but no HbA. Microcytic red cells typical. Clinically similar to SS disease.",
    "S_beta_plus_classic" = "S-β⁺ thalassemia pattern: Moderate HgbS with some HbA and elevated HgbA2. Microcytic red cells. Generally milder than SS or S-β⁰.",
    "AS_classic" = "Classic sickle cell trait pattern: Majority HbA (50-65%) with moderate HgbS (35-45%). Usually asymptomatic carrier state. Important for genetic counseling."
  )
  
  return(explanations[[rule_name]] %||% "Standard clinical interpretation based on electrophoresis pattern.")
}