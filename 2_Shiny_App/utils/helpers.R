# ==============================================================================
# Helper Functions Module
# ==============================================================================
# UI rendering, data preparation, and utility functions for the Shiny app
# ==============================================================================

library(shiny)
library(bslib)
library(DT)
library(plotly)

#' Prepare input data from Shiny inputs
#' @param input Shiny input object
#' @return List with processed patient data
prepare_input_data <- function(input) {
  
  # Extract all input values
  data <- list(
    # Required electrophoresis values
    hgb_s_pct = input$hgb_s_pct,
    hgb_a_pct = input$hgb_a_pct,
    hgb_c_pct = input$hgb_c_pct,
    hgb_f_pct = input$hgb_f_pct,
    hgb_a2_pct = input$hgb_a2_pct,
    
    # Required transfusion status
    is_transfused = input$is_transfused,
    
    # Optional CBC parameters
    total_hgb = if(input$include_cbc) input$total_hgb else NULL,
    mcv = if(input$include_cbc) input$mcv else NULL,
    mchc = if(input$include_cbc) input$mchc else NULL,
    rdw = if(input$include_cbc) input$rdw else NULL,
    hematocrit = if(input$include_cbc) input$hematocrit else NULL,
    rbc_count = if(input$include_cbc) input$rbc_count else NULL,
    
    # Optional quality parameters
    test_confidence_level = if(input$test_confidence_level != "") input$test_confidence_level else NULL,
    clarity_score = input$clarity_score,
    
    # Patient info
    patient_id = if(input$patient_id != "") input$patient_id else NULL
  )
  
  return(data)
}

#' Run ML prediction on patient data (Enhanced 5-class with AS)
#' @param input_data List containing patient data
#' @param model Trained tidymodels workflow
#' @return List with ML prediction results
predict_with_ml <- function(input_data, model) {
  
  # Prepare data for ML prediction
  ml_data <- prepare_ml_data(input_data)
  
  # Get predictions and probabilities
  predictions <- predict(model, new_data = ml_data, type = "class")
  probabilities <- predict(model, new_data = ml_data, type = "prob")
  
  # Apply clinical constraints to probabilities
  probabilities_constrained <- apply_clinical_constraints(probabilities, input_data)
  
  # Recalculate prediction based on constrained probabilities
  max_prob_col <- probabilities_constrained %>%
    pivot_longer(everything(), names_to = "class", values_to = "prob") %>%
    slice_max(prob, n = 1) %>%
    pull(class)
  
  # Extract the class name from column name (remove .pred_ prefix)
  final_prediction <- str_remove(max_prob_col, "^\\.pred_")
  final_confidence <- probabilities_constrained %>%
    select(all_of(max_prob_col)) %>%
    pull(1)
  
  # Combine results and handle 5-class predictions
  result <- list(
    prediction = final_prediction,
    confidence = final_confidence,
    needs_expert_review = final_confidence <= 0.90,
    probabilities = probabilities_constrained,
    source = "machine_learning"
  )
  
  return(result)
}

#' Apply clinical constraints to ML probabilities
#' @param probabilities Data frame with raw ML probabilities
#' @param input_data List containing patient data
#' @return Data frame with clinically constrained probabilities
apply_clinical_constraints <- function(probabilities, input_data) {
  
  # Extract key values
  hgb_s_pct <- as.numeric(input_data$hgb_s_pct)
  hgb_a_pct <- as.numeric(input_data$hgb_a_pct)
  hgb_c_pct <- as.numeric(input_data$hgb_c_pct)
  hgb_f_pct <- as.numeric(input_data$hgb_f_pct)
  hgb_a2_pct <- as.numeric(input_data$hgb_a2_pct)
  
  # Create a copy to modify
  prob_constrained <- probabilities
  
  # Constraint 1: SC requires significant HgbC (>10% minimum, ideally >40%)
  if (hgb_c_pct < 10) {
    if (".pred_SC" %in% names(prob_constrained)) {
      # Heavily penalize SC when HgbC is very low
      penalty_factor <- if (hgb_c_pct == 0) 0.01 else (hgb_c_pct / 10) * 0.1
      prob_constrained$.pred_SC <- prob_constrained$.pred_SC * penalty_factor
    }
  }
  
  # Constraint 2: SS should have minimal HgbA (<5%)
  if (hgb_a_pct > 10) {
    if (".pred_SS" %in% names(prob_constrained)) {
      # Penalize SS when HgbA is high
      penalty_factor <- max(0.1, 1 - (hgb_a_pct - 5) / 50)
      prob_constrained$.pred_SS <- prob_constrained$.pred_SS * penalty_factor
    }
  }
  
  # Constraint 3: Boost AS when pattern is classic trait (HgbA 40-70%, HgbS 30-60%, HgbC ≈ 0)
  if (hgb_a_pct >= 40 && hgb_a_pct <= 70 && 
      hgb_s_pct >= 30 && hgb_s_pct <= 60 && 
      hgb_c_pct <= 2) {
    if (".pred_AS" %in% names(prob_constrained)) {
      # Boost AS for classic trait pattern
      boost_factor <- 1.5
      prob_constrained$.pred_AS <- prob_constrained$.pred_AS * boost_factor
    }
  }
  
  # Constraint 4: S-beta patterns require elevated HgbA2 (>3.5%, ideally >5.7%)
  if (hgb_a2_pct < 3.5) {
    if (".pred_S_beta_0" %in% names(prob_constrained)) {
      penalty_factor <- 0.2
      prob_constrained$.pred_S_beta_0 <- prob_constrained$.pred_S_beta_0 * penalty_factor
    }
    if (".pred_S_beta_plus" %in% names(prob_constrained)) {
      penalty_factor <- 0.2
      prob_constrained$.pred_S_beta_plus <- prob_constrained$.pred_S_beta_plus * penalty_factor
    }
  }
  
  # Constraint 5: S-beta+ requires some HgbA (>5%)
  if (hgb_a_pct < 5) {
    if (".pred_S_beta_plus" %in% names(prob_constrained)) {
      penalty_factor <- 0.1
      prob_constrained$.pred_S_beta_plus <- prob_constrained$.pred_S_beta_plus * penalty_factor
    }
  }
  
  # Renormalize probabilities to sum to 1
  total_prob <- rowSums(prob_constrained, na.rm = TRUE)
  if (total_prob > 0) {
    prob_constrained <- prob_constrained / total_prob
  }
  
  return(prob_constrained)
}

#' Prepare data for ML prediction (similar to training pipeline)
#' @param input_data List containing patient data
#' @return Data frame formatted for ML prediction
prepare_ml_data <- function(input_data) {
  
  # Create data frame with same structure as training data
  ml_df <- data.frame(
    # Primary electrophoresis values
    hgb_s_pct = as.numeric(input_data$hgb_s_pct),
    hgb_a_pct = as.numeric(input_data$hgb_a_pct),
    hgb_c_pct = as.numeric(input_data$hgb_c_pct),
    hgb_f_pct = as.numeric(input_data$hgb_f_pct),
    hgb_a2_pct = as.numeric(input_data$hgb_a2_pct),
    
    # CBC parameters (with defaults if missing)
    total_hgb = if(!is.null(input_data$total_hgb)) as.numeric(input_data$total_hgb) else 10,
    hematocrit = if(!is.null(input_data$hematocrit)) as.numeric(input_data$hematocrit) else 30,
    mcv = if(!is.null(input_data$mcv)) as.numeric(input_data$mcv) else 85,
    mchc = if(!is.null(input_data$mchc)) as.numeric(input_data$mchc) else 34,
    rdw = if(!is.null(input_data$rdw)) as.numeric(input_data$rdw) else 15,
    rbc_count = if(!is.null(input_data$rbc_count)) as.numeric(input_data$rbc_count) else 4.5,
    
    # Derived features
    hgb_a_present = as.numeric(as.numeric(input_data$hgb_a_pct) > 0),
    hgb_s_present = as.numeric(as.numeric(input_data$hgb_s_pct) > 0),
    hgb_c_present = as.numeric(as.numeric(input_data$hgb_c_pct) > 0),
    hgb_f_elevated = as.numeric(as.numeric(input_data$hgb_f_pct) > 2),
    hgb_a2_elevated = as.numeric(as.numeric(input_data$hgb_a2_pct) > 3.5),
    
    # Ratios
    hgb_s_to_total = as.numeric(input_data$hgb_s_pct) / 100,
    hgb_a_to_s_ratio = if(as.numeric(input_data$hgb_s_pct) > 0) {
      as.numeric(input_data$hgb_a_pct) / as.numeric(input_data$hgb_s_pct)
    } else { 0 },
    total_variant_pct = as.numeric(input_data$hgb_s_pct) + as.numeric(input_data$hgb_c_pct),
    
    # Quality scores
    completeness_score = 10,
    clarity_score = if(!is.null(input_data$clarity_score)) as.numeric(input_data$clarity_score) else 8,
    total_quality_score = if(!is.null(input_data$clarity_score)) as.numeric(input_data$clarity_score) else 8,
    
    # Historical indicator
    had_hba_below_1_ever = as.numeric(as.numeric(input_data$hgb_a_pct) < 1),
    
    # TRANSFUSION STATUS
    is_transfused = as.numeric(input_data$is_transfused)
  )
  
  return(ml_df)
}

#' Compile final prediction result from rules and ML
#' @param rule_result List with rule-based prediction
#' @param ml_result List with ML prediction (can be NULL)
#' @param input_data List with original input data
#' @return List with final prediction result
compile_final_result <- function(rule_result, ml_result, input_data) {
  
  if (rule_result$high_confidence) {
    # Use rule-based prediction
    result <- list(
      prediction = rule_result$prediction,
      confidence = rule_result$confidence,
      source = "clinical_rules",
      rule_applied = rule_result$rule_applied,
      rationale = rule_result$rationale,
      needs_expert_review = FALSE,
      probabilities = NULL,
      trait_flag = rule_result$trait_flag  # Preserve trait flag from rules
    )
  } else if (!is.null(ml_result)) {
    # Use ML prediction but PRESERVE trait flag from rules
    result <- list(
      prediction = ml_result$prediction,
      confidence = ml_result$confidence,
      source = "machine_learning",
      rule_applied = NULL,
      rationale = "Machine learning prediction based on all available features",
      needs_expert_review = ml_result$needs_expert_review,
      probabilities = ml_result$probabilities,
      trait_flag = rule_result$trait_flag  # IMPORTANT: Preserve trait flag even for ML predictions
    )
  } else {
    # Ambiguous case
    result <- list(
      prediction = "AMBIGUOUS",
      confidence = 0,
      source = "none",
      rule_applied = NULL,
      rationale = "Unable to classify - requires expert review",
      needs_expert_review = TRUE,
      probabilities = NULL,
      trait_flag = rule_result$trait_flag  # Preserve trait flag even for ambiguous cases
    )
  }
  
  # Add clinical warnings
  warnings <- check_clinical_warnings(input_data)
  result$clinical_warnings <- warnings
  
  return(result)
}

#' Render prediction results UI
#' @param prediction_result List with prediction information
#' @param input_data List with input data
#' @return Shiny UI element
render_prediction_results <- function(prediction_result, input_data) {
  
  if (is.null(prediction_result)) {
    return(div(
      class = "alert alert-info",
      icon("info-circle"),
      h4("Ready for Prediction"),
      p("Enter patient data and click 'Predict Genotype' to begin analysis.")
    ))
  }
  
  # Determine card styling based on confidence
  card_class <- if (prediction_result$source == "clinical_rules") {
    "border-success"
  } else if (prediction_result$confidence > 0.9) {
    "border-success"
  } else if (prediction_result$confidence > 0.7) {
    "border-warning"
  } else {
    "border-danger"
  }
  
  div(
    class = "row",
    
    # Main prediction card
    div(
      class = "col-12",
      div(
        class = paste("card", card_class),
        div(
          class = "card-header d-flex justify-content-between align-items-center",
          h4("Prediction Result", class = "mb-0"),
          span(
            class = paste("badge", if(prediction_result$needs_expert_review) "bg-warning" else "bg-success"),
            if(prediction_result$needs_expert_review) "Expert Review Recommended" else "Automated Classification"
          )
        ),
        div(
          class = "card-body",
          
          # Main prediction
          div(
            class = "row mb-3",
            div(
              class = "col-md-6",
              # Display genotype with better formatting
              h2(
                case_when(
                  prediction_result$prediction == "SS" ~ "SS (Sickle Cell Disease)",
                  prediction_result$prediction == "SC" ~ "SC (Sickle-C Disease)", 
                  prediction_result$prediction == "S_beta_0" ~ "Sβ⁰ (Sickle Beta-Zero Thalassemia)",
                  prediction_result$prediction == "S_beta_plus" ~ "Sβ⁺ (Sickle Beta-Plus Thalassemia)",
                  prediction_result$prediction == "AS" ~ "AS (Sickle Cell Trait)",
                  prediction_result$prediction == "AS_borderline" ~ "AS (Sickle Cell Trait - Borderline)",
                  TRUE ~ prediction_result$prediction
                ),
                class = "text-primary mb-2"
              ),
              p(class = "lead", paste("Confidence:", round(prediction_result$confidence * 100, 1), "%")),
              p(class = "text-muted", paste("Source:", 
                case_when(
                  prediction_result$source == "clinical_rules" ~ "Clinical Rules",
                  prediction_result$source == "machine_learning" ~ "Machine Learning",
                  TRUE ~ "Manual Review Required"
                )
              ))
            ),
            div(
              class = "col-md-6",
              h5("Clinical Rationale:"),
              p(prediction_result$rationale)
            )
          ),
          
          # Expert review warning if needed
          if (prediction_result$needs_expert_review) {
            div(
              class = "alert alert-warning mb-3",
              icon("exclamation-triangle"),
              strong(" Expert Review Recommended"),
              br(),
              "Low confidence prediction or ambiguous features require expert interpretation."
            )
          },
          
          # Rule information if applicable
          if (!is.null(prediction_result$rule_applied)) {
            div(
              class = "alert alert-info mb-3",
              icon("check-circle"),
              strong(" Clinical Rule Applied: "), prediction_result$rule_applied,
              br(),
              get_rule_explanation(prediction_result$rule_applied)
            )
          },
          
          # Clinical warnings if any
          if (!is.null(prediction_result$clinical_warnings) && 
              prediction_result$clinical_warnings$has_warnings) {
            div(
              class = "alert alert-warning mb-3",
              icon("exclamation-triangle"),
              strong(" Clinical Considerations:"),
              tags$ul(
                lapply(prediction_result$clinical_warnings$warnings, function(w) {
                  tags$li(w)
                })
              )
            )
          },
          
          # NEW: Trait flag warning - important for genetic counseling
          if (!is.null(prediction_result$trait_flag)) {
            div(
              class = "alert alert-info mb-3",
              style = "border-left: 4px solid #17a2b8; background-color: #f8f9fa;",
              icon("info-circle"),
              strong("Additional Clinical Consideration:"),
              br(),
              span("Electrophoresis pattern suggests possible sickle cell trait (AS):"),
              br(),
              span(class = "text-muted", prediction_result$trait_flag$message),
              br(),
              br(),
              div(
                class = "mt-2 p-2",
                style = "background-color: #e9ecef; border-radius: 4px; border-left: 3px solid #6c757d;",
                strong("Clinical Recommendations:"),
                tags$ul(
                  class = "mb-0 mt-1",
                  tags$li("Consider genetic counseling for family planning"),
                  tags$li("Review family history for sickle cell disease"),
                  tags$li("Document trait possibility in medical record"),
                  tags$li(paste0("Pattern confidence: ", round(prediction_result$trait_flag$confidence * 100, 1), "%"))
                )
              )
            )
          }
        )
      )
    ),
    
    # Probability table if ML was used
    if (!is.null(prediction_result$probabilities)) {
      div(
        class = "col-12 mt-3",
        div(
          class = "card",
          div(class = "card-header", h5("Prediction Probabilities")),
          div(
            class = "card-body",
            render_probability_table(prediction_result$probabilities)
          )
        )
      )
    }
  )
}

#' Render probability table
#' @param probabilities Data frame with prediction probabilities
#' @return DT datatable
render_probability_table <- function(probabilities) {
  
  # Convert to long format for display
  prob_df <- probabilities %>%
    pivot_longer(everything(), names_to = "Genotype", values_to = "Probability") %>%
    mutate(
      Genotype = str_remove(Genotype, "^\\.pred_"),
      # Convert genotype codes to readable names
      Genotype = case_when(
        Genotype == "SS" ~ "SS (Sickle Cell Disease)",
        Genotype == "SC" ~ "SC (Sickle-C Disease)",
        Genotype == "S_beta_0" ~ "Sβ⁰ (Sickle Beta-Zero Thalassemia)",
        Genotype == "S_beta_plus" ~ "Sβ⁺ (Sickle Beta-Plus Thalassemia)",
        Genotype == "AS" ~ "AS (Sickle Cell Trait)",
        TRUE ~ Genotype
      ),
      Probability = round(Probability, 3),
      Percentage = paste0(round(Probability * 100, 1), "%")
    ) %>%
    arrange(desc(Probability))
  
  datatable(
    prob_df %>% select(Genotype, Percentage),
    options = list(
      dom = 't',
      paging = FALSE,
      searching = FALSE,
      ordering = FALSE
    ),
    rownames = FALSE
  )
}

#' Render SHAP explanation UI
#' @param shap_data List with SHAP analysis results
#' @param prediction_result List with prediction information
#' @return Shiny UI element
render_shap_explanation <- function(shap_data, prediction_result) {
  
  if (is.null(shap_data)) {
    return(div(
      class = "alert alert-info",
      icon("info-circle"),
      h4("Explanation"),
      p("Explanations will appear here when machine learning predictions are made."),
      if (!is.null(prediction_result) && prediction_result$source == "clinical_rules") {
        p("This case was classified using clinical rules, so further analysis is not needed.")
      } else {
        p("Machine learning prediction was made but SHAP analysis is still processing or failed.")
      }
    ))
  }
  
  # Check if we have feature contributions (minimum requirement)
  if (is.null(shap_data$feature_contributions) || nrow(shap_data$feature_contributions) == 0) {
    return(div(
      class = "alert alert-warning",
      icon("exclamation-triangle"),
      h4("SHAP Analysis Unavailable"),
      p("SHAP analysis could not be completed for this prediction."),
      p("The machine learning prediction is still valid and reliable.")
    ))
  }
  
  div(
    # SHAP waterfall plot (only if available)
    if (!is.null(shap_data$waterfall_plot)) {
      div(
        class = "mb-4",
        h4("Feature Contributions"),
        plotOutput("shap_waterfall", height = "400px")
      )
    } else {
      div(
        class = "alert alert-info mb-4",
        h4("Waterfall Plot"),
        p("Waterfall plot could not be generated, but feature contributions are available below.")
      )
    },
    
    # SHAP insights (should always be available if we have contributions)
    div(
      class = "mb-4",
      h4("Key Insights"),
      div(
        class = "alert alert-light",
        if (!is.null(shap_data$feature_contributions) && nrow(shap_data$feature_contributions) > 0) {
          tags$ul(
            lapply(summarize_shap_insights(shap_data$feature_contributions, prediction_result), 
                   function(insight) tags$li(HTML(insight)))
          )
        } else {
          p("Feature insights are not available.")
        }
      )
    ),
    
    # Feature contributions table (should always be available)
    div(
      h4("Detailed Feature Contributions"),
      if (!is.null(shap_data$feature_contributions) && nrow(shap_data$feature_contributions) > 0) {
        DT::dataTableOutput("shap_contributions_table")
      } else {
        div(
          class = "alert alert-info",
          p("Feature contributions table is not available.")
        )
      }
    )
  )
}

#' Render clinical guidelines content
#' @return Shiny UI element
render_clinical_guidelines <- function() {
  div(
    h3("Clinical Decision Rules"),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("SS (Sickle Cell Disease)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbS > 75%"),
          tags$li("HgbA = 0%"),
          tags$li("HgbC = 0%"),
          tags$li("HgbF < 20%"),
          tags$li("HgbA2 < 3.5%"),
          tags$li("MCV typically > 80 fL")
        ),
        p(class = "text-muted", "Most severe form with frequent complications.")
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("SC (Sickle-Hemoglobin C Disease)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbS 45-55%"),
          tags$li("HgbC 40-55%"),
          tags$li("HgbA = 0%"),
          tags$li("HgbF < 5%"),
          tags$li("MCV 70-90 fL")
        ),
        p(class = "text-muted", "Generally milder than SS but can have complications.")
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("S-β⁰ (Sickle Beta-Zero Thalassemia)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbS > 70%"),
          tags$li("HgbA = 0%"),
          tags$li("HgbA2 > 5.7%"),
          tags$li("HgbF < 20%"),
          tags$li("MCV < 80 fL (microcytic)")
        ),
        p(class = "text-muted", "Clinically similar to SS disease.")
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("S-β⁺ (Sickle Beta-Plus Thalassemia)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbS 60-90%"),
          tags$li("HgbA 5-30%"),
          tags$li("HgbA2 > 5.7%"),
          tags$li("HgbF < 20%"),
          tags$li("MCV < 80 fL (microcytic)")
        ),
        p(class = "text-muted", "Generally milder than SS or S-β⁰."),
        div(
          class = "alert alert-warning",
          strong("Note: "), "Transfusion confounds HgbA levels in this pattern."
        )
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("AS (Sickle Cell Trait)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbA 50-65%"),
          tags$li("HgbS 35-45%"),
          tags$li("HgbA2 < 3.5%"),
          tags$li("HgbF < 2%"),
          tags$li("Normal MCV")
        ),
        p(class = "text-muted", "Usually asymptomatic carrier state."),
        div(
          class = "alert alert-warning",
          strong("Note: "), "Transfusion confounds HgbA levels in this pattern."
        )
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("AA (Normal Hemoglobin)")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("HgbA 92-98%"),
          tags$li("HgbS 0%"),
          tags$li("HgbC 0%"),
          tags$li("HgbA2 1.5-3.5%"),
          tags$li("HgbF < 2%")
        ),
        p(class = "text-muted", "Normal hemoglobin pattern - no sickle cell disease or trait.")
      )
    ),
    
    div(
      class = "alert alert-info",
      h5("Transfusion Considerations"),
      p("Recent transfusion (within 3 months) can significantly alter HgbA levels, making certain patterns unreliable for diagnosis. The system accounts for this by:"),
      tags$ul(
        tags$li("Flagging AS and S-β⁺ patterns as ambiguous when transfusion is present"),
        tags$li("Using transfusion status as a feature in machine learning predictions"),
        tags$li("Including transfusion interactions in the analysis")
      )
    )
  )
}

#' Render about content
#' @return Shiny UI element
render_about_content <- function() {
  div(
    h3("About This System"),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("Overview")),
      div(
        class = "card-body",
        p("This system combines clinical decision rules with machine learning to predict hemoglobin genotypes from electrophoresis and CBC data."),
        p("The system uses a two-tier approach:"),
        tags$ol(
          tags$li("Clinical rules for high-confidence cases"),
          tags$li("Machine learning with SHAP explanations for ambiguous cases")
        )
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("Model Information")),
      div(
        class = "card-body",
        tags$ul(
          tags$li("Algorithm: Enhanced XGBoost with tidymodels"),
          tags$li("Classes: 5-class prediction (SS, SC, S-β⁰, S-β⁺, AS)"),
          tags$li("Training data: Enhanced dataset with AS patients"),
          tags$li("Cross-validation accuracy: ~95% on 5-class problem"),
          tags$li("Transfusion-aware features and interactions"),
          tags$li("SHAP explanations for interpretability")
        )
      )
    ),
    
    div(
      class = "card mb-3",
      div(class = "card-header", h5("Important Disclaimers")),
      div(
        class = "card-body",
        div(
          class = "alert alert-warning",
          tags$ul(
            tags$li("This tool is for research and educational purposes only"),
            tags$li("Clinical decisions should always involve expert interpretation"),
            tags$li("Recent transfusion affects interpretation accuracy"),
            tags$li("Unusual patterns may require additional testing")
          )
        )
      )
    ),
    
    div(
      class = "card",
      div(class = "card-header", h5("SHAP Explanations")),
      div(
        class = "card-body",
        p("SHAP (SHapley Additive exPlanations) values show how each feature contributes to the prediction:"),
        tags$ul(
          tags$li("Positive SHAP values support the predicted class"),
          tags$li("Negative SHAP values oppose the predicted class"),
          tags$li("Larger absolute values indicate stronger influence"),
          tags$li("The waterfall plot shows the cumulative effect of all features")
        )
      )
    )
  )
}

#' Generate prediction report for download
#' @param prediction_result List with prediction information
#' @param shap_data List with SHAP analysis results
#' @param input_data List with input data
#' @param file Character, output file path
generate_prediction_report <- function(prediction_result, shap_data, input_data, file) {
  
  # Create HTML report
  html_content <- paste0(
    "<!DOCTYPE html>",
    "<html><head><title>Genotype Prediction Report</title>",
    "<style>body{font-family:Arial,sans-serif;margin:40px;}",
    ".result{background:#f8f9fa;padding:20px;border-radius:5px;margin:10px 0;}",
    ".warning{background:#fff3cd;padding:15px;border-radius:5px;margin:10px 0;}",
    "table{border-collapse:collapse;width:100%;}",
    "th,td{border:1px solid #ddd;padding:8px;text-align:left;}",
    "th{background-color:#f2f2f2;}</style></head><body>",
    
    "<h1>Hemoglobin Genotype Prediction Report</h1>",
    "<p><strong>Generated on:</strong> ", Sys.time(), "</p>",
    
    if (!is.null(input_data$patient_id)) paste0("<p><strong>Patient ID:</strong> ", input_data$patient_id, "</p>"),
    
    "<h2>Input Data</h2>",
    "<table>",
    "<tr><th>Parameter</th><th>Value</th></tr>",
    "<tr><td>HgbS %</td><td>", input_data$hgb_s_pct, "</td></tr>",
    "<tr><td>HgbA %</td><td>", input_data$hgb_a_pct, "</td></tr>",
    "<tr><td>HgbC %</td><td>", input_data$hgb_c_pct, "</td></tr>",
    "<tr><td>HgbF %</td><td>", input_data$hgb_f_pct, "</td></tr>",
    "<tr><td>HgbA2 %</td><td>", input_data$hgb_a2_pct, "</td></tr>",
    "<tr><td>Transfused</td><td>", ifelse(input_data$is_transfused == 1, "Yes", "No"), "</td></tr>",
    "</table>",
    
    "<h2>Prediction Result</h2>",
    "<div class='result'>",
    "<h3>", prediction_result$prediction, "</h3>",
    "<p><strong>Confidence:</strong> ", round(prediction_result$confidence * 100, 1), "%</p>",
    "<p><strong>Source:</strong> ", prediction_result$source, "</p>",
    "<p><strong>Rationale:</strong> ", prediction_result$rationale, "</p>",
    "</div>",
    
    if (prediction_result$needs_expert_review) {
      "<div class='warning'><strong>Expert Review Recommended:</strong> Low confidence prediction requires expert interpretation.</div>"
    } else { "" },
    
    "</body></html>"
  )
  
  writeLines(html_content, file)
}

#' Create a simple feature importance plot from input values when SHAP fails
#' @param input_data List with input data
#' @return ggplot object
create_simple_feature_plot <- function(input_data) {
  
  # Create a simple feature value plot
  feature_data <- data.frame(
    Feature = c("HgbS %", "HgbA %", "HgbC %", "HgbF %", "HgbA2 %", "Transfused"),
    Value = c(
      as.numeric(input_data$hgb_s_pct),
      as.numeric(input_data$hgb_a_pct),
      as.numeric(input_data$hgb_c_pct),
      as.numeric(input_data$hgb_f_pct),
      as.numeric(input_data$hgb_a2_pct),
      as.numeric(input_data$is_transfused) * 100  # Convert to percentage for display
    ),
    Type = c("Hemoglobin", "Hemoglobin", "Hemoglobin", "Hemoglobin", "Hemoglobin", "Clinical")
  )
  
  ggplot(feature_data, aes(x = reorder(Feature, Value), y = Value, fill = Type)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Input Feature Values",
      subtitle = "(SHAP importance not available)",
      x = "Features",
      y = "Value (%)",
      fill = "Feature Type"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Hemoglobin" = "#3498db", "Clinical" = "#e74c3c"))
}

#' Check for clinical warnings based on input data
#' @param input_data List containing patient data
#' @return List with clinical warnings
check_clinical_warnings <- function(input_data) {
  
  warnings <- list()
  has_warnings <- FALSE
  
  # Check for very low hemoglobin
  if (!is.null(input_data$total_hgb) && as.numeric(input_data$total_hgb) < 6) {
    warnings <- append(warnings, "Very low total hemoglobin - severe anemia present")
    has_warnings <- TRUE
  }
  
  # Check for unusual electrophoresis totals
  total_hgb_pct <- as.numeric(input_data$hgb_s_pct) + as.numeric(input_data$hgb_a_pct) + 
                   as.numeric(input_data$hgb_c_pct) + as.numeric(input_data$hgb_f_pct) + 
                   as.numeric(input_data$hgb_a2_pct)
  
  if (total_hgb_pct < 90 || total_hgb_pct > 110) {
    warnings <- append(warnings, paste0("Unusual electrophoresis total (", round(total_hgb_pct, 1), "%) - verify results"))
    has_warnings <- TRUE
  }
  
  # Check for very high HgbF
  if (as.numeric(input_data$hgb_f_pct) > 25) {
    warnings <- append(warnings, "Very high HgbF - consider hereditary persistence of fetal hemoglobin")
    has_warnings <- TRUE
  }
  
  # Check for transfusion with unusual patterns
  if (as.numeric(input_data$is_transfused) == 1 && as.numeric(input_data$hgb_a_pct) < 30) {
    warnings <- append(warnings, "Low HgbA despite recent transfusion - verify transfusion history")
    has_warnings <- TRUE
  }
  
  return(list(
    has_warnings = has_warnings,
    warnings = warnings
  ))
} 


# ==============================================================================
# BATCH SPREADSHEET ANALYSIS MODULE (WITH FULL SHAP TRANSPARENCY)
# ==============================================================================

SHARP_BATCH_COLUMN_MAP <- list(
  patient_id            = c("patient_id", "patient_mrn", "mrn", "id", "sample_id", "subject_id"),
  hgb_s_pct             = c("hgb_s_pct", "hgb s pct", "hgb s %", "hgb s", "hbs", "hbs %", "hemoglobin s", "hgb s quant"),
  hgb_a_pct             = c("hgb_a_pct", "hgb a pct", "hgb a %", "hgb a", "hba", "hba %", "hemoglobin a", "hgb a quant"),
  hgb_c_pct             = c("hgb_c_pct", "hgb c pct", "hgb c %", "hgb c", "hbc", "hbc %", "hemoglobin c", "hgb c quant"),
  hgb_f_pct             = c("hgb_f_pct", "hgb f pct", "hgb f %", "hgb f", "hbf", "hbf %", "hemoglobin f", "hgb f quant"),
  hgb_a2_pct            = c("hgb_a2_pct", "hgb a2 pct", "hgb a2 %", "hgb a2", "hba2", "hba2 %", "hemoglobin a2", "hgb a2 quant"),
  is_transfused         = c("is_transfused", "transfused", "transfusion_status", "recent_transfusion", "transfusion"),
  total_hgb             = c("total_hgb", "hgb", "hemoglobin", "hemoglobin-blood", "hb"),
  mcv                   = c("mcv", "mean cell volume", "mean corpuscular volume", "mcv-mean cell volume fl"),
  mchc                  = c("mchc", "mean cell hemoglobin concentration", "mean corpuscular hemoglobin concentration"),
  rdw                   = c("rdw", "red cell distribution width"),
  hematocrit            = c("hematocrit", "hct", "hematocrit-blood"),
  rbc_count             = c("rbc_count", "rbc", "red blood cell count", "rbc count-blood"),
  clarity_score         = c("clarity_score", "clarity"),
  test_confidence_level = c("test_confidence_level", "test_quality", "confidence_level")
)

#' Generate template dataframe for SHARP batch upload
#' @return Tibble with required and optional columns and realistic example patients
get_sharp_batch_template <- function() {
  tibble::tibble(
    patient_id    = c("PT_001_SS", "PT_002_Transfused_SS", "PT_003_SC", "PT_004_Sbeta0", "PT_005_SbetaPlus", "PT_006_AS_Trait", "PT_007_AA_Normal"),
    hgb_s_pct     = c(84.5, 58.0, 48.5, 78.0, 66.0, 38.5, 0.0),
    hgb_a_pct     = c(0.0,  28.5, 0.0,  0.0,  18.5, 57.5, 96.8),
    hgb_c_pct     = c(0.0,  0.0,  46.0, 0.0,  0.0,  0.0,  0.0),
    hgb_f_pct     = c(12.3, 10.2, 2.5,  15.2, 9.0,  0.8,  0.5),
    hgb_a2_pct    = c(3.2,  3.3,  3.0,  6.8,  6.5,  3.2,  2.7),
    is_transfused = c(0,    1,    0,    0,    0,    0,    0),
    total_hgb     = c(7.8,  9.2,  10.9, 8.1,  9.8,  13.6, 14.5),
    mcv           = c(89.0, 86.5, 79.5, 68.5, 71.0, 87.0, 90.0),
    mchc          = c(34.2, 33.8, 34.8, 32.1, 32.8, 34.0, 34.1),
    rdw           = c(20.5, 19.2, 16.0, 19.8, 17.5, 13.1, 12.6),
    hematocrit    = c(23.0, 27.5, 32.0, 25.0, 29.5, 40.5, 43.0),
    rbc_count     = c(2.6,  3.2,  4.0,  3.7,  4.2,  4.7,  4.8),
    clarity_score = c(9,    8,    9,    8,    8,    9,    10)
  )
}

#' Standardize uploaded batch dataframe column names
#' @param df Raw uploaded dataframe
#' @return List with `data` (standardized tibble), `mapped` (character vector), `missing_required` (character vector)
standardize_sharp_batch_columns <- function(df) {
  current_names <- names(df)
  names_clean <- trimws(tolower(gsub("[_\\.]+", " ", current_names)))
  names_exact <- trimws(tolower(current_names))

  std_df <- df
  mapped_info <- c()

  for (std_name in names(SHARP_BATCH_COLUMN_MAP)) {
    aliases <- tolower(SHARP_BATCH_COLUMN_MAP[[std_name]])
    aliases_clean <- gsub("[_\\.]+", " ", aliases)
    idx <- which(names_exact %in% aliases | names_clean %in% aliases_clean)
    if (length(idx) > 0) {
      orig_col <- current_names[idx[1]]
      std_df[[std_name]] <- df[[orig_col]]
      mapped_info <- c(mapped_info, sprintf("%s -> %s", orig_col, std_name))
    }
  }

  # Default hgb_c_pct to 0 and is_transfused to 0 if omitted, but require S, A, F, A2
  if (!("hgb_c_pct" %in% names(std_df))) {
    std_df$hgb_c_pct <- 0
  }
  if (!("is_transfused" %in% names(std_df))) {
    std_df$is_transfused <- 0
  }

  required_cols <- c("hgb_s_pct", "hgb_a_pct", "hgb_f_pct", "hgb_a2_pct")
  missing_req <- setdiff(required_cols, names(std_df))

  return(list(
    data = std_df,
    mapped = mapped_info,
    missing_required = missing_req
  ))
}

#' Run full SHARP two-tier prediction + exact TreeSHAP batch analysis
#' @param raw_df Uploaded dataframe
#' @param model Trained tidymodels workflow
#' @param progress_callback Optional function(value, detail) for Shiny progress bar
#' @return List with `results_df`, `summary_stats`, and `error`
run_sharp_batch_analysis <- function(raw_df, model, progress_callback = NULL) {
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    return(list(results_df = NULL, summary_stats = NULL, error = "Uploaded spreadsheet is empty."))
  }

  if (!is.null(progress_callback)) progress_callback(0.1, "Validating and standardizing columns...")
  std_res <- standardize_sharp_batch_columns(raw_df)
  if (length(std_res$missing_required) > 0) {
    return(list(
      results_df = NULL,
      summary_stats = NULL,
      error = paste0(
        "Missing required electrophoresis columns: ",
        paste(std_res$missing_required, collapse = ", "),
        ". Please download the Template CSV to see the exact required column names."
      )
    ))
  }

  df <- std_res$data
  n_rows <- nrow(df)

  # Parse numeric columns cleanly
  num_cols <- c("hgb_s_pct", "hgb_a_pct", "hgb_c_pct", "hgb_f_pct", "hgb_a2_pct",
                "total_hgb", "mcv", "mchc", "rdw", "hematocrit", "rbc_count", "clarity_score")
  for (col in intersect(num_cols, names(df))) {
    df[[col]] <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(df[[col]]))))
  }
  for (col in c("hgb_s_pct", "hgb_a_pct", "hgb_c_pct", "hgb_f_pct", "hgb_a2_pct")) {
    df[[col]] <- ifelse(is.na(df[[col]]), 0, df[[col]])
  }

  # Parse transfusion indicator (accepts 0/1, Yes/No, True/False)
  raw_tf <- tolower(trimws(as.character(df$is_transfused)))
  df$is_transfused <- ifelse(raw_tf %in% c("1", "yes", "y", "true", "t", "transfused"), 1, 0)

  if (!("patient_id" %in% names(df))) {
    df$patient_id <- paste0("Row_", seq_len(n_rows))
  }

  if (!is.null(progress_callback)) progress_callback(0.35, "Running Tier-1 Clinical Rules & Tier-2 ML + Biochemical Constraints...")

  final_pred       <- character(n_rows)
  pred_label_clean <- character(n_rows)
  pred_source      <- character(n_rows)
  rule_applied_vec <- character(n_rows)
  confidence_vec   <- numeric(n_rows)
  expert_rev_vec   <- logical(n_rows)
  rationale_vec    <- character(n_rows)

  prob_SS          <- numeric(n_rows)
  prob_SC          <- numeric(n_rows)
  prob_S_beta_0    <- numeric(n_rows)
  prob_S_beta_plus <- numeric(n_rows)
  prob_AS          <- numeric(n_rows)

  for (i in seq_len(n_rows)) {
    row_list <- as.list(df[i, ])
    rule_res <- apply_clinical_rules(row_list)
    ml_res   <- predict_with_ml(row_list, model)

    probs_i <- ml_res$probabilities
    prob_SS[i]          <- if (".pred_SS" %in% names(probs_i)) round(probs_i$.pred_SS[1], 4) else NA_real_
    prob_SC[i]          <- if (".pred_SC" %in% names(probs_i)) round(probs_i$.pred_SC[1], 4) else NA_real_
    prob_S_beta_0[i]    <- if (".pred_S_beta_0" %in% names(probs_i)) round(probs_i$.pred_S_beta_0[1], 4) else NA_real_
    prob_S_beta_plus[i] <- if (".pred_S_beta_plus" %in% names(probs_i)) round(probs_i$.pred_S_beta_plus[1], 4) else NA_real_
    prob_AS[i]          <- if (".pred_AS" %in% names(probs_i)) round(probs_i$.pred_AS[1], 4) else NA_real_

    if (isTRUE(rule_res$high_confidence)) {
      final_pred[i]     <- rule_res$prediction
      pred_source[i]    <- "Clinical Rules"
      rule_applied_vec[i] <- rule_res$rule_applied %||% "Clinical_Rule"
      confidence_vec[i] <- round(rule_res$confidence, 4)
      expert_rev_vec[i] <- FALSE
      rationale_vec[i]  <- rule_res$rationale
    } else {
      final_pred[i]     <- ml_res$prediction
      pred_source[i]    <- "Machine Learning"
      rule_applied_vec[i] <- "None (Tier-2 ML)"
      confidence_vec[i] <- round(ml_res$confidence, 4)
      expert_rev_vec[i] <- isTRUE(ml_res$needs_expert_review)
      rationale_vec[i]  <- "Tier-2 5-class XGBoost prediction with biochemical plausibility constraints"
    }

    pred_label_clean[i] <- dplyr::case_when(
      final_pred[i] == "SS" ~ "HbSS (Sickle Cell Anemia)",
      final_pred[i] == "SC" ~ "HbSC (Sickle-Hemoglobin C Disease)",
      final_pred[i] == "S_beta_0" ~ "HbSβ⁰ (Sickle Beta-Zero Thalassemia)",
      final_pred[i] == "S_beta_plus" ~ "HbSβ⁺ (Sickle Beta-Plus Thalassemia)",
      final_pred[i] == "AS" ~ "HbAS (Sickle Cell Trait)",
      grepl("^AA", final_pred[i]) ~ "HbAA (Normal Hemoglobin)",
      TRUE ~ final_pred[i]
    )
  }

  if (!is.null(progress_callback)) progress_callback(0.70, "Computing exact multi-class TreeSHAP feature contributions for all rows...")

  shap_cols_tbl <- compute_batch_shap_sharp(df, model, final_pred)

  if (!is.null(progress_callback)) progress_callback(0.90, "Assembling output spreadsheet with SHAP transparency columns...")

  pred_summary_tbl <- tibble::tibble(
    predicted_genotype        = final_pred,
    predicted_genotype_label  = pred_label_clean,
    prediction_confidence     = confidence_vec,
    prediction_confidence_pct = round(confidence_vec * 100, 1),
    prediction_source         = pred_source,
    clinical_rule_applied     = rule_applied_vec,
    needs_expert_review       = expert_rev_vec,
    clinical_rationale        = rationale_vec,
    prob_SS                   = prob_SS,
    prob_SC                   = prob_SC,
    prob_S_beta_0             = prob_S_beta_0,
    prob_S_beta_plus          = prob_S_beta_plus,
    prob_AS                   = prob_AS
  )

  results_df <- dplyr::bind_cols(raw_df, pred_summary_tbl, shap_cols_tbl)

  summary_stats <- list(
    total_rows          = n_rows,
    rule_classified     = sum(pred_source == "Clinical Rules"),
    ml_classified       = sum(pred_source == "Machine Learning"),
    automated_calls     = sum(!expert_rev_vec),
    expert_review_calls = sum(expert_rev_vec),
    coverage_pct        = round(mean(!expert_rev_vec) * 100, 1),
    mean_confidence_pct = round(mean(confidence_vec, na.rm = TRUE) * 100, 1),
    genotype_counts     = table(final_pred)
  )

  if (!is.null(progress_callback)) progress_callback(1.0, "Batch analysis complete!")

  return(list(
    results_df = results_df,
    summary_stats = summary_stats,
    error = NULL
  ))
}
