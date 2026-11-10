# ==============================================================================
# Data Validation Module
# ==============================================================================
# Functions to validate user inputs for hemoglobin genotype prediction
# ==============================================================================

#' Validate all user inputs
#' @param input Shiny input object
#' @return List with valid (logical) and errors (character vector)
validate_inputs <- function(input) {
  errors <- c()
  
  # Check required electrophoresis values
  required_fields <- list(
    "hgb_s_pct" = "HgbS %",
    "hgb_a_pct" = "HgbA %", 
    "hgb_f_pct" = "HgbF %",
    "hgb_a2_pct" = "HgbA2 %"
  )
  
  for (field in names(required_fields)) {
    if (is.null(input[[field]]) || is.na(input[[field]]) || input[[field]] == "") {
      errors <- c(errors, paste("Missing required field:", required_fields[[field]]))
    }
  }
  
  # Check transfusion status
  if (is.null(input$is_transfused) || input$is_transfused == "") {
    errors <- c(errors, "Transfusion status is required")
  }
  
  # Validate numeric ranges
  if (!is.null(input$hgb_s_pct) && !is.na(input$hgb_s_pct)) {
    if (input$hgb_s_pct < 0 || input$hgb_s_pct > 100) {
      errors <- c(errors, "HgbS% must be between 0-100")
    }
  }
  
  if (!is.null(input$hgb_a_pct) && !is.na(input$hgb_a_pct)) {
    if (input$hgb_a_pct < 0 || input$hgb_a_pct > 100) {
      errors <- c(errors, "HgbA% must be between 0-100")
    }
  }
  
  if (!is.null(input$hgb_c_pct) && !is.na(input$hgb_c_pct)) {
    if (input$hgb_c_pct < 0 || input$hgb_c_pct > 100) {
      errors <- c(errors, "HgbC% must be between 0-100")
    }
  }
  
  if (!is.null(input$hgb_f_pct) && !is.na(input$hgb_f_pct)) {
    if (input$hgb_f_pct < 0 || input$hgb_f_pct > 100) {
      errors <- c(errors, "HgbF% must be between 0-100")
    }
  }
  
  if (!is.null(input$hgb_a2_pct) && !is.na(input$hgb_a2_pct)) {
    if (input$hgb_a2_pct < 0 || input$hgb_a2_pct > 100) {
      errors <- c(errors, "HgbA2% must be between 0-100")
    }
  }
  
  # Check that percentages sum reasonably (allowing some measurement error)
  if (all(!is.na(c(input$hgb_s_pct, input$hgb_a_pct, input$hgb_c_pct, 
                   input$hgb_f_pct, input$hgb_a2_pct)))) {
    total_pct <- sum(c(input$hgb_s_pct, input$hgb_a_pct, input$hgb_c_pct, 
                       input$hgb_f_pct, input$hgb_a2_pct), na.rm = TRUE)
    if (total_pct > 110) {
      errors <- c(errors, "Total hemoglobin percentages exceed 110% - please check values")
    }
    if (total_pct < 85) {
      errors <- c(errors, "Total hemoglobin percentages below 85% - please check values")
    }
  }
  
  # Validate CBC parameters if included
  if (!is.null(input$include_cbc) && input$include_cbc) {
    if (!is.null(input$total_hgb) && !is.na(input$total_hgb)) {
      if (input$total_hgb < 3 || input$total_hgb > 20) {
        errors <- c(errors, "Total hemoglobin must be between 3-20 g/dL")
      }
    }
    
    if (!is.null(input$mcv) && !is.na(input$mcv)) {
      if (input$mcv < 50 || input$mcv > 120) {
        errors <- c(errors, "MCV must be between 50-120 fL")
      }
    }
    
    if (!is.null(input$mchc) && !is.na(input$mchc)) {
      if (input$mchc < 25 || input$mchc > 40) {
        errors <- c(errors, "MCHC must be between 25-40 g/dL")
      }
    }
    
    if (!is.null(input$rdw) && !is.na(input$rdw)) {
      if (input$rdw < 10 || input$rdw > 30) {
        errors <- c(errors, "RDW must be between 10-30%")
      }
    }
    
    if (!is.null(input$hematocrit) && !is.na(input$hematocrit)) {
      if (input$hematocrit < 15 || input$hematocrit > 60) {
        errors <- c(errors, "Hematocrit must be between 15-60%")
      }
    }
    
    if (!is.null(input$rbc_count) && !is.na(input$rbc_count)) {
      if (input$rbc_count < 2 || input$rbc_count > 8) {
        errors <- c(errors, "RBC count must be between 2-8 M/μL")
      }
    }
  }
  
  # Validate clarity score if provided
  if (!is.null(input$clarity_score) && !is.na(input$clarity_score)) {
    if (input$clarity_score < 1 || input$clarity_score > 10) {
      errors <- c(errors, "Clarity score must be between 1-10")
    }
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors
  ))
}

#' Check for critical clinical combinations that need attention
#' @param input Shiny input object
#' @return List with warnings and clinical notes
check_clinical_warnings <- function(input) {
  warnings <- c()
  
  # Check for transfusion with HbA present (confounding)
  if (!is.null(input$is_transfused) && input$is_transfused == 1) {
    if (!is.null(input$hgb_a_pct) && !is.na(input$hgb_a_pct) && input$hgb_a_pct > 5) {
      warnings <- c(warnings, 
                   "Patient is transfused with HbA present - consider impact on genotype interpretation")
    }
  }
  
  # Check for extremely high HbF (possible HPFH)
  if (!is.null(input$hgb_f_pct) && !is.na(input$hgb_f_pct) && input$hgb_f_pct > 30) {
    warnings <- c(warnings, 
                 "Very high HbF level - consider hereditary persistence of fetal hemoglobin (HPFH)")
  }
  
  # Check for very low total hemoglobin with sickle syndromes
  if (!is.null(input$total_hgb) && !is.na(input$total_hgb) && input$total_hgb < 7) {
    if (!is.null(input$hgb_s_pct) && !is.na(input$hgb_s_pct) && input$hgb_s_pct > 50) {
      warnings <- c(warnings, 
                   "Very low hemoglobin with high HbS - consider severe sickle syndrome")
    }
  }
  
  # Check for borderline values that might affect classification
  if (!is.null(input$hgb_s_pct) && !is.na(input$hgb_s_pct)) {
    if (input$hgb_s_pct >= 45 && input$hgb_s_pct <= 55) {
      warnings <- c(warnings, 
                   "HbS in borderline range (45-55%) - prediction may have higher uncertainty")
    }
  }
  
  return(list(
    has_warnings = length(warnings) > 0,
    warnings = warnings
  ))
}

#' Get reference ranges for hemoglobin variants
#' @return Data frame with reference ranges
get_reference_ranges <- function() {
  data.frame(
    Parameter = c("HbA", "HbS", "HbC", "HbF", "HbA2", "Total Hgb", "MCV", "MCHC"),
    Normal_Range = c("95-98%", "0%", "0%", "<2%", "2.0-3.5%", "12-16 g/dL", "80-100 fL", "32-36 g/dL"),
    SS_Typical = c("0%", ">80%", "0%", "5-15%", "2-4%", "6-9 g/dL", ">80 fL", "32-36 g/dL"),
    SC_Typical = c("0%", "45-50%", "45-50%", "<5%", "2-4%", "9-12 g/dL", "70-90 fL", "32-36 g/dL"),
    AS_Typical = c("55-65%", "35-45%", "0%", "<2%", "2-4%", "Normal", "Normal", "Normal"),
    stringsAsFactors = FALSE
  )
} 