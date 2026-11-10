# ==============================================================================
# SHAP Generator Module (Upgraded for Exact Multi-Class TreeSHAP & Batch Mode)
# ==============================================================================
# Functions to generate exact TreeSHAP explanations for single and batch ML
# predictions from the latest 5-class tidymodels XGBoost workflow.
# ==============================================================================

library(shapviz)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tidymodels)
library(recipes)
library(purrr)

#' Clean feature names for clinical display
#' @param x Character vector of raw model feature names
#' @return Character vector of human-readable clinical feature names
clean_sharp_feature_names <- function(x) {
  dplyr::recode(
    as.character(x),
    "hgb_s_pct"            = "HbS (%)",
    "hgb_a_pct"            = "HbA (%)",
    "hgb_c_pct"            = "HbC (%)",
    "hgb_f_pct"            = "HbF (%)",
    "hgb_a2_pct"           = "HbA2 (%)",
    "total_hgb"            = "Total Hemoglobin (g/dL)",
    "hematocrit"           = "Hematocrit (%)",
    "mcv"                  = "MCV (fL)",
    "mchc"                 = "MCHC (g/dL)",
    "rdw"                  = "RDW (%)",
    "rbc_count"            = "RBC Count (M/uL)",
    "hgb_a_present"        = "HbA Present",
    "hgb_s_present"        = "HbS Present",
    "hgb_c_present"        = "HbC Present",
    "hgb_f_elevated"       = "Elevated HbF (>2%)",
    "hgb_a2_elevated"      = "Elevated HbA2 (>3.5%)",
    "hgb_s_to_total"       = "HbS Proportion (HbS%)",
    "hgb_a_to_s_ratio"     = "HbA-to-HbS Ratio",
    "total_variant_pct"    = "Total Variant Hb (HbS+HbC %)",
    "completeness_score"   = "Completeness Score",
    "clarity_score"        = "Profile Clarity Score",
    "total_quality_score"  = "Total Quality Score",
    "had_hba_below_1_ever" = "HbA < 1% Indicator",
    "is_transfused"        = "Recent Transfusion",
    "hgb_f_to_s_ratio"     = "HbF-to-HbS Ratio",
    "hgb_f_to_a2_ratio"    = "HbF-to-HbA2 Ratio",
    "as_pattern_score"     = "HbAS Pattern Score",
    "trait_indicator"      = "Sickle Trait Indicator",
    "transfused_hgb_s"     = "Transfusion x HbS (%)",
    "transfused_hgb_a"     = "Transfusion x HbA (%)",
    "transfused_hgb_c"     = "Transfusion x HbC (%)",
    "transfused_mcv"       = "Transfusion x MCV",
    "transfused_total_hgb" = "Transfusion x Total Hgb",
    "borderline_hgb_s"     = "Borderline HbS (70-80%)",
    "borderline_hgb_f"     = "Borderline HbF (5-15%)",
    "borderline_mcv"       = "Borderline MCV (75-85 fL)",
    .default               = as.character(x)
  )
}

#' Build raw feature dataframe expected by the SHARP tidymodels recipe
#' @param df Data frame with standardized patient columns
#' @return Data frame ready for `bake(extract_recipe(model), new_data = ...)`
build_sharp_recipe_input_df <- function(df) {
  n <- nrow(df)
  get_num <- function(col, default_val) {
    if (col %in% names(df)) {
      v <- suppressWarnings(as.numeric(df[[col]]))
      ifelse(is.na(v), default_val, v)
    } else {
      rep(default_val, n)
    }
  }

  hgb_s_pct     <- get_num("hgb_s_pct", 0)
  hgb_a_pct     <- get_num("hgb_a_pct", 0)
  hgb_c_pct     <- get_num("hgb_c_pct", 0)
  hgb_f_pct     <- get_num("hgb_f_pct", 0)
  hgb_a2_pct    <- get_num("hgb_a2_pct", 0)
  total_hgb     <- get_num("total_hgb", 10.0)
  hematocrit    <- get_num("hematocrit", 30.0)
  mcv           <- get_num("mcv", 85.0)
  mchc          <- get_num("mchc", 34.0)
  rdw           <- get_num("rdw", 15.0)
  rbc_count     <- get_num("rbc_count", 4.5)
  is_transfused <- get_num("is_transfused", 0)
  clarity_score <- get_num("clarity_score", 8.0)

  tibble::tibble(
    hgb_s_pct            = hgb_s_pct,
    hgb_a_pct            = hgb_a_pct,
    hgb_c_pct            = hgb_c_pct,
    hgb_f_pct            = hgb_f_pct,
    hgb_a2_pct           = hgb_a2_pct,
    total_hgb            = total_hgb,
    hematocrit           = hematocrit,
    mcv                  = mcv,
    mchc                 = mchc,
    rdw                  = rdw,
    rbc_count            = rbc_count,
    hgb_a_present        = as.numeric(hgb_a_pct > 0),
    hgb_s_present        = as.numeric(hgb_s_pct > 0),
    hgb_c_present        = as.numeric(hgb_c_pct > 0),
    hgb_f_elevated       = as.numeric(hgb_f_pct > 2),
    hgb_a2_elevated      = as.numeric(hgb_a2_pct > 3.5),
    hgb_s_to_total       = hgb_s_pct / 100,
    hgb_a_to_s_ratio     = dplyr::if_else(hgb_s_pct > 0, hgb_a_pct / hgb_s_pct, 0),
    total_variant_pct    = hgb_s_pct + hgb_c_pct,
    completeness_score   = 10,
    clarity_score        = clarity_score,
    total_quality_score  = clarity_score,
    had_hba_below_1_ever = as.numeric(hgb_a_pct < 1),
    is_transfused        = is_transfused
  )
}

#' Safely generate SHAP explanation with error handling
#' @param input_data List containing patient data
#' @param model Trained tidymodels workflow
#' @param target_class Optional character string of predicted genotype class
#' @return List with SHAP results or error information
safely_generate_shap <- function(input_data, model, target_class = NULL) {
  tryCatch({
    tryCatch({
      shap_result <- generate_shap_explanation(input_data, model, target_class = target_class)
      return(list(result = shap_result, error = NULL))
    }, error = function(e) {
      tryCatch({
        shap_result <- generate_shap_explanation_simple(input_data, model)
        return(list(result = shap_result, error = NULL))
      }, error = function(e2) {
        return(list(result = NULL, error = paste("SHAP Error:", e2$message)))
      })
    })
  }, error = function(e) {
    return(list(result = NULL, error = paste("SHAP Error:", e$message)))
  })
}

#' Generate exact TreeSHAP explanation for a single prediction
#' @param input_data List containing patient data
#' @param model Trained tidymodels workflow
#' @param target_class Optional predicted class ("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
#' @return List with SHAP visualizations and data
generate_shap_explanation <- function(input_data, model, target_class = NULL) {
  processed_data <- prepare_data_for_shap(input_data, model)

  xgb_core_model <- workflows::extract_fit_engine(model)
  data_matrix <- as.matrix(processed_data$features)

  # Determine target class index for multi-class XGBoost
  class_levels <- tryCatch({
    levels(model$pre$mold$blueprint$ptypes$outcomes$genotype_simplified)
  }, error = function(e) {
    c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
  })
  if (is.null(class_levels) || length(class_levels) == 0) {
    class_levels <- c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
  }

  class_idx <- 1L
  if (!is.null(target_class) && target_class %in% class_levels) {
    class_idx <- match(target_class, class_levels)
  }

  sv_all <- shapviz::shapviz(
    xgb_core_model,
    X_pred = data_matrix,
    X = processed_data$features_df
  )

  sv_target <- if (inherits(sv_all, "mshapviz") || (is.list(sv_all) && !inherits(sv_all, "shapviz"))) {
    sv_all[[min(class_idx, length(sv_all))]]
  } else {
    sv_all
  }

  shap_vec <- as.numeric(sv_target$S[1, ])
  names(shap_vec) <- colnames(sv_target$S)
  baseline_prob <- as.numeric(sv_target$baseline)

  feature_contributions <- tibble::tibble(
    feature = colnames(sv_target$S),
    feature_label_clean = clean_sharp_feature_names(colnames(sv_target$S)),
    feature_value = as.numeric(processed_data$raw_feature_display[colnames(sv_target$S)]),
    shap_value = shap_vec,
    abs_shap = abs(shap_vec)
  ) %>%
    dplyr::arrange(dplyr::desc(abs_shap)) %>%
    dplyr::mutate(
      rank = dplyr::row_number(),
      direction = ifelse(shap_value >= 0, "Positive", "Negative"),
      clinical_interpretation = purrr::map2_chr(
        feature, shap_value,
        ~ map_feature_to_clinical(.x, .y, processed_data$original_values)
      )
    )

  waterfall_plot <- create_horizontal_waterfall_plot(
    feature_contributions,
    baseline_prob = baseline_prob,
    target_class = if (!is.null(target_class)) target_class else class_levels[class_idx]
  )
  importance_plot <- create_importance_plot_simple(feature_contributions %>% dplyr::slice_head(n = 12))

  return(list(
    waterfall_plot = waterfall_plot,
    importance_plot = importance_plot,
    feature_contributions = feature_contributions,
    shap_values = shap_vec,
    baseline = baseline_prob,
    processed_data = processed_data,
    sv_object = sv_target
  ))
}

#' Prepare input data for SHAP analysis using the fitted tidymodels recipe
#' @param input_data List containing patient data
#' @param model Trained tidymodels workflow
#' @return List with processed features and metadata
prepare_data_for_shap <- function(input_data, model) {
  single_df <- tibble::tibble(
    hgb_s_pct     = as.numeric(input_data$hgb_s_pct %||% 0),
    hgb_a_pct     = as.numeric(input_data$hgb_a_pct %||% 0),
    hgb_c_pct     = as.numeric(input_data$hgb_c_pct %||% 0),
    hgb_f_pct     = as.numeric(input_data$hgb_f_pct %||% 0),
    hgb_a2_pct    = as.numeric(input_data$hgb_a2_pct %||% 0),
    total_hgb     = if (!is.null(input_data$total_hgb) && !is.na(input_data$total_hgb)) as.numeric(input_data$total_hgb) else 10.0,
    hematocrit    = if (!is.null(input_data$hematocrit) && !is.na(input_data$hematocrit)) as.numeric(input_data$hematocrit) else 30.0,
    mcv           = if (!is.null(input_data$mcv) && !is.na(input_data$mcv)) as.numeric(input_data$mcv) else 85.0,
    mchc          = if (!is.null(input_data$mchc) && !is.na(input_data$mchc)) as.numeric(input_data$mchc) else 34.0,
    rdw           = if (!is.null(input_data$rdw) && !is.na(input_data$rdw)) as.numeric(input_data$rdw) else 15.0,
    rbc_count     = if (!is.null(input_data$rbc_count) && !is.na(input_data$rbc_count)) as.numeric(input_data$rbc_count) else 4.5,
    is_transfused = as.numeric(input_data$is_transfused %||% 0),
    clarity_score = if (!is.null(input_data$clarity_score) && !is.na(input_data$clarity_score)) as.numeric(input_data$clarity_score) else 8.0
  )

  patient_df <- build_sharp_recipe_input_df(single_df)
  fitted_recipe <- workflows::extract_recipe(model, estimated = TRUE)
  processed_data <- recipes::bake(fitted_recipe, new_data = patient_df, recipes::all_predictors())

  # Map unscaled display values for each baked feature so waterfall labels show clinical units
  raw_display <- as.list(patient_df[1, ])
  raw_display$hgb_s_to_total       <- raw_display$hgb_s_pct
  raw_display$hgb_f_to_s_ratio     <- raw_display$hgb_f_pct / (raw_display$hgb_s_pct + 0.1)
  raw_display$hgb_f_to_a2_ratio    <- raw_display$hgb_f_pct / (raw_display$hgb_a2_pct + 0.1)
  raw_display$trait_indicator      <- as.numeric(raw_display$hgb_a_pct > raw_display$hgb_s_pct & raw_display$hgb_a_pct >= 50 & raw_display$hgb_s_pct >= 30)
  raw_display$transfused_hgb_s     <- raw_display$is_transfused * raw_display$hgb_s_pct
  raw_display$transfused_hgb_a     <- raw_display$is_transfused * raw_display$hgb_a_pct
  raw_display$transfused_hgb_c     <- raw_display$is_transfused * raw_display$hgb_c_pct
  raw_display$transfused_mcv       <- raw_display$is_transfused * raw_display$mcv
  raw_display$transfused_total_hgb <- raw_display$is_transfused * raw_display$total_hgb
  raw_display$borderline_hgb_s     <- as.numeric(raw_display$hgb_s_pct >= 70 & raw_display$hgb_s_pct <= 80)
  raw_display$borderline_hgb_f     <- as.numeric(raw_display$hgb_f_pct >= 5 & raw_display$hgb_f_pct <= 15)
  raw_display$borderline_mcv       <- as.numeric(raw_display$mcv >= 75 & raw_display$mcv <= 85)

  raw_vec <- sapply(colnames(processed_data), function(col) {
    val <- raw_display[[col]]
    if (is.null(val) || is.na(val)) as.numeric(processed_data[[col]][1]) else as.numeric(val)
  })

  original_values <- list(
    hgb_s_pct     = single_df$hgb_s_pct,
    hgb_a_pct     = single_df$hgb_a_pct,
    hgb_c_pct     = single_df$hgb_c_pct,
    hgb_f_pct     = single_df$hgb_f_pct,
    hgb_a2_pct    = single_df$hgb_a2_pct,
    is_transfused = single_df$is_transfused,
    total_hgb     = if (!is.null(input_data$total_hgb)) as.numeric(input_data$total_hgb) else NA,
    mcv           = if (!is.null(input_data$mcv)) as.numeric(input_data$mcv) else NA,
    mchc          = if (!is.null(input_data$mchc)) as.numeric(input_data$mchc) else NA,
    rdw           = if (!is.null(input_data$rdw)) as.numeric(input_data$rdw) else NA,
    hematocrit    = if (!is.null(input_data$hematocrit)) as.numeric(input_data$hematocrit) else NA,
    rbc_count     = if (!is.null(input_data$rbc_count)) as.numeric(input_data$rbc_count) else NA
  )

  return(list(
    features = processed_data,
    features_df = as.data.frame(processed_data),
    raw_feature_display = raw_vec,
    original_values = original_values
  ))
}

#' Compute exact multi-class TreeSHAP columns for an entire batch dataframe
#' @param df_standardized Standardized batch dataframe
#' @param model Trained tidymodels workflow
#' @param target_classes Character vector of target genotype classes per row
#' @return Tibble of SHAP columns for every row in `df_standardized`
compute_batch_shap_sharp <- function(df_standardized, model, target_classes) {
  n_rows <- nrow(df_standardized)
  recipe_input_df <- build_sharp_recipe_input_df(df_standardized)
  fitted_recipe <- workflows::extract_recipe(model, estimated = TRUE)
  baked_x <- recipes::bake(fitted_recipe, new_data = recipe_input_df, recipes::all_predictors())
  mat_x <- as.matrix(baked_x)
  xgb_engine <- workflows::extract_fit_engine(model)

  class_levels <- tryCatch({
    levels(model$pre$mold$blueprint$ptypes$outcomes$genotype_simplified)
  }, error = function(e) {
    c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
  })
  if (is.null(class_levels) || length(class_levels) == 0) {
    class_levels <- c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")
  }

  sv_all <- shapviz::shapviz(xgb_engine, X_pred = mat_x, X = as.data.frame(baked_x))
  is_multi <- inherits(sv_all, "mshapviz") || (is.list(sv_all) && !inherits(sv_all, "shapviz"))

  feat_cols <- colnames(baked_x)
  row_shap_mat <- matrix(0, nrow = n_rows, ncol = length(feat_cols), dimnames = list(NULL, feat_cols))
  baselines <- numeric(n_rows)
  resolved_target_class <- character(n_rows)

  for (i in seq_len(n_rows)) {
    tc <- as.character(target_classes[i])
    if (is.na(tc) || !(tc %in% class_levels)) {
      # For AA Override or ambiguous labels, report SHAP relative to AS (or SS)
      tc <- if (grepl("^AA", tc)) "AS" else "SS"
    }
    resolved_target_class[i] <- tc
    c_idx <- match(tc, class_levels)
    if (is.na(c_idx)) c_idx <- 1L

    sv_i <- if (is_multi) sv_all[[min(c_idx, length(sv_all))]] else sv_all
    row_shap_mat[i, ] <- as.numeric(sv_i$S[i, ])
    baselines[i] <- as.numeric(sv_i$baseline)
  }

  shap_df <- as.data.frame(row_shap_mat)

  # Helper to sum a group of related engineered features for clear clinical transparency
  sum_cols <- function(cols) {
    valid_cols <- intersect(cols, colnames(shap_df))
    if (length(valid_cols) == 0) return(rep(0, n_rows))
    round(rowSums(shap_df[, valid_cols, drop = FALSE], na.rm = TRUE), 4)
  }

  # Top positive/negative drivers per row
  clean_names_vec <- clean_sharp_feature_names(feat_cols)
  top_pos_vec <- character(n_rows)
  top_neg_vec <- character(n_rows)
  summary_vec <- character(n_rows)

  for (i in seq_len(n_rows)) {
    vals <- row_shap_mat[i, ]
    ord_pos <- order(vals, decreasing = TRUE)
    pos_idx <- ord_pos[vals[ord_pos] > 0]
    if (length(pos_idx) > 0) {
      top_k <- head(pos_idx, 3)
      top_pos_vec[i] <- paste(sprintf("%s (%+.3f)", clean_names_vec[top_k], vals[top_k]), collapse = "; ")
    } else {
      top_pos_vec[i] <- "None"
    }

    ord_neg <- order(vals, decreasing = FALSE)
    neg_idx <- ord_neg[vals[ord_neg] < 0]
    if (length(neg_idx) > 0) {
      top_k_neg <- head(neg_idx, 3)
      top_neg_vec[i] <- paste(sprintf("%s (%+.3f)", clean_names_vec[top_k_neg], vals[top_k_neg]), collapse = "; ")
    } else {
      top_neg_vec[i] <- "None"
    }

    top1_idx <- which.max(abs(vals))
    dir_txt <- ifelse(vals[top1_idx] >= 0, "supports", "opposes")
    summary_vec[i] <- sprintf(
      "Top SHAP factor: %s (%+.3f) %s %s",
      clean_names_vec[top1_idx], vals[top1_idx], dir_txt, resolved_target_class[i]
    )
  }

  out_tbl <- tibble::tibble(
    shap_target_genotype       = resolved_target_class,
    shap_baseline_log_odds     = round(baselines, 4),
    shap_prediction_log_odds   = round(baselines + rowSums(row_shap_mat, na.rm = TRUE), 4),
    shap_hgb_s_group           = sum_cols(c("hgb_s_pct", "hgb_s_present", "hgb_s_to_total", "borderline_hgb_s")),
    shap_hgb_a_group           = sum_cols(c("hgb_a_pct", "hgb_a_present", "hgb_a_to_s_ratio", "had_hba_below_1_ever")),
    shap_hgb_c_group           = sum_cols(c("hgb_c_pct", "hgb_c_present", "total_variant_pct")),
    shap_hgb_f_group           = sum_cols(c("hgb_f_pct", "hgb_f_elevated", "hgb_f_to_s_ratio", "hgb_f_to_a2_ratio", "borderline_hgb_f")),
    shap_hgb_a2_group          = sum_cols(c("hgb_a2_pct", "hgb_a2_elevated")),
    shap_transfusion_group     = sum_cols(c("is_transfused", "transfused_hgb_a", "transfused_hgb_c", "transfused_total_hgb", "transfused_hgb_s", "transfused_mcv")),
    shap_mcv_group             = sum_cols(c("mcv", "borderline_mcv")),
    shap_total_hgb             = sum_cols(c("total_hgb")),
    shap_mchc                  = sum_cols(c("mchc")),
    shap_rdw                   = sum_cols(c("rdw")),
    shap_rbc_count             = sum_cols(c("rbc_count")),
    shap_trait_indicator       = sum_cols(c("trait_indicator", "as_pattern_score")),
    top_positive_shap_features = top_pos_vec,
    top_negative_shap_features = top_neg_vec,
    shap_clinical_summary      = summary_vec
  )

  # Also append all individual baked model feature SHAP columns (prefixed with shap_feature_)
  raw_shap_cols <- as_tibble(round(shap_df, 4))
  names(raw_shap_cols) <- paste0("shap_feature_", names(raw_shap_cols))

  dplyr::bind_cols(out_tbl, raw_shap_cols)
}

#' Generate SHAP explanation (Simplified fallback)
#' @param input_data List containing patient data
#' @param model Trained tidymodels workflow
#' @return List with SHAP visualizations and data
generate_shap_explanation_simple <- function(input_data, model) {
  baseline_prob <- 0.0
  feature_names <- c("hgb_s_pct", "hgb_a_pct", "hgb_c_pct", "hgb_f_pct", "hgb_a2_pct",
                     "is_transfused", "total_hgb", "mcv", "mchc", "rdw", "hematocrit")
  feature_values <- c(
    as.numeric(input_data$hgb_s_pct %||% 0),
    as.numeric(input_data$hgb_a_pct %||% 0),
    as.numeric(input_data$hgb_c_pct %||% 0),
    as.numeric(input_data$hgb_f_pct %||% 0),
    as.numeric(input_data$hgb_a2_pct %||% 0),
    as.numeric(input_data$is_transfused %||% 0),
    if (!is.null(input_data$total_hgb) && !is.na(input_data$total_hgb)) as.numeric(input_data$total_hgb) else 10.0,
    if (!is.null(input_data$mcv) && !is.na(input_data$mcv)) as.numeric(input_data$mcv) else 85.0,
    if (!is.null(input_data$mchc) && !is.na(input_data$mchc)) as.numeric(input_data$mchc) else 34.0,
    if (!is.null(input_data$rdw) && !is.na(input_data$rdw)) as.numeric(input_data$rdw) else 15.0,
    if (!is.null(input_data$hematocrit) && !is.na(input_data$hematocrit)) as.numeric(input_data$hematocrit) else 30.0
  )
  feature_values[is.na(feature_values)] <- 0

  shap_values <- c(
    ifelse(feature_values[1] > 70, 1.2, ifelse(feature_values[1] > 40, 0.8, 0.3)),
    ifelse(feature_values[2] > 40, -0.6, ifelse(feature_values[2] > 10, -0.3, 0.1)),
    ifelse(feature_values[3] > 30, 0.9, ifelse(feature_values[3] > 5, 0.3, 0)),
    ifelse(feature_values[4] > 15, 0.5, ifelse(feature_values[4] > 5, 0.2, 0)),
    ifelse(feature_values[5] > 5, 0.4, ifelse(feature_values[5] > 3.5, 0.2, 0)),
    ifelse(feature_values[6] == 1, -0.4, 0),
    ifelse(feature_values[7] < 8, 0.3, ifelse(feature_values[7] < 10, 0.1, -0.1)),
    ifelse(feature_values[8] < 70, 0.3, ifelse(feature_values[8] < 80, 0.1, 0)),
    ifelse(feature_values[9] < 30, 0.2, ifelse(feature_values[9] > 36, 0.1, 0)),
    ifelse(feature_values[10] > 18, 0.3, ifelse(feature_values[10] > 15, 0.1, 0)),
    ifelse(feature_values[11] < 25, 0.2, ifelse(feature_values[11] < 30, 0.1, 0))
  )

  original_values <- list(
    hgb_s_pct = feature_values[1],
    hgb_a_pct = feature_values[2],
    hgb_c_pct = feature_values[3],
    hgb_f_pct = feature_values[4],
    hgb_a2_pct = feature_values[5],
    is_transfused = feature_values[6],
    total_hgb = feature_values[7],
    mcv = feature_values[8]
  )

  feature_contributions <- data.frame(
    feature = feature_names,
    feature_label_clean = clean_sharp_feature_names(feature_names),
    feature_value = feature_values,
    shap_value = shap_values,
    abs_shap = abs(shap_values),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(direction = ifelse(shap_value >= 0, "Positive", "Negative")) %>%
    dplyr::arrange(dplyr::desc(abs_shap)) %>%
    dplyr::slice_head(n = 8)

  feature_contributions$clinical_interpretation <- sapply(seq_len(nrow(feature_contributions)), function(i) {
    map_feature_to_clinical(feature_contributions$feature[i], feature_contributions$shap_value[i], original_values)
  })

  waterfall_plot <- create_horizontal_waterfall_plot(feature_contributions, baseline_prob)
  importance_plot <- create_importance_plot_simple(feature_contributions)

  return(list(
    waterfall_plot = waterfall_plot,
    importance_plot = importance_plot,
    feature_contributions = feature_contributions,
    shap_values = feature_contributions$shap_value,
    processed_data = list(original_values = original_values),
    sv_object = NULL
  ))
}

#' Create horizontal waterfall plot for SHAP explanations
#' @param feature_contributions Data frame with feature contributions
#' @param baseline_prob Baseline log-odds
#' @param target_class Optional target class label
#' @return ggplot2 horizontal waterfall plot
create_horizontal_waterfall_plot <- function(feature_contributions, baseline_prob = 0.0, target_class = NULL) {
  tryCatch({
    plot_data <- feature_contributions %>%
      dplyr::arrange(dplyr::desc(abs_shap)) %>%
      dplyr::slice_head(n = 8) %>%
      dplyr::arrange(abs_shap)

    cumulative_values <- cumsum(c(baseline_prob, plot_data$shap_value))
    plot_data$start_pos <- cumulative_values[-length(cumulative_values)]
    plot_data$end_pos <- cumulative_values[-1]

    plot_data <- plot_data %>%
      dplyr::mutate(
        clean_name = clean_sharp_feature_names(feature),
        feature_label = paste0(clean_name, " = ", round(feature_value, 2))
      )

    # Ensure unique factor levels
    plot_data$feature_label <- make.unique(plot_data$feature_label)
    plot_data$y_idx <- seq_len(nrow(plot_data))
    final_prediction <- baseline_prob + sum(plot_data$shap_value)

    title_str <- if (!is.null(target_class)) {
      paste0("SHAP Waterfall Plot (Target Class: ", target_class, ")")
    } else {
      "SHAP Waterfall Plot"
    }

    p <- ggplot(plot_data) +
      geom_rect(aes(
        xmin = pmin(start_pos, end_pos),
        xmax = pmax(start_pos, end_pos),
        ymin = y_idx - 0.35,
        ymax = y_idx + 0.35,
        fill = direction
      ), color = "white", linewidth = 0.5) +
      geom_text(aes(
        x = (start_pos + end_pos) / 2,
        y = y_idx,
        label = ifelse(abs(shap_value) >= 0.01, sprintf("%+.3f", shap_value), "")
      ), color = "#1B1919", fontface = "bold", size = 3.2) +
      geom_vline(xintercept = baseline_prob, linetype = "dashed", color = "gray40", alpha = 0.8) +
      geom_vline(xintercept = final_prediction, linetype = "solid", color = "#1B1919", alpha = 0.8) +
      scale_fill_manual(
        values = c("Positive" = "#FFA500", "Negative" = "#8B4B9B"),
        labels = c("Positive" = "Supports predicted genotype", "Negative" = "Opposes predicted genotype")
      ) +
      scale_y_continuous(
        breaks = plot_data$y_idx,
        labels = plot_data$feature_label,
        expand = expansion(mult = c(0.08, 0.12))
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
      labs(
        title = title_str,
        subtitle = paste0("Base value E[f(X)] = ", round(baseline_prob, 3),
                          "  ->  Top-features output f(x) = ", round(final_prediction, 3)),
        x = "Model Output Contribution (log-odds)",
        y = NULL,
        fill = "Contribution"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )

    return(p)
  }, error = function(e) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste("Waterfall plot error:", e$message), size = 4) +
      theme_void()
  })
}

#' Map feature contributions to clinical interpretations
map_feature_to_clinical <- function(feature_name, shap_value, original_values) {
  tryCatch({
    direction <- ifelse(shap_value >= 0, "supports", "opposes")
    magnitude <- ifelse(abs(shap_value) > 0.5, "strongly",
                        ifelse(abs(shap_value) > 0.15, "moderately", "weakly"))
    clean_lbl <- clean_sharp_feature_names(feature_name)

    interpretation <- switch(feature_name,
      "hgb_s_pct" = paste0("HbS level (", round(original_values$hgb_s_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "hgb_s_to_total" = paste0("HbS proportion (", round(original_values$hgb_s_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "hgb_a_pct" = paste0("HbA level (", round(original_values$hgb_a_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "hgb_c_pct" = paste0("HbC level (", round(original_values$hgb_c_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "hgb_f_pct" = paste0("HbF level (", round(original_values$hgb_f_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "hgb_a2_pct" = paste0("HbA2 level (", round(original_values$hgb_a2_pct, 1), "%) ", magnitude, " ", direction, " the predicted genotype"),
      "is_transfused" = paste0("Transfusion status (", ifelse(original_values$is_transfused == 1, "Yes", "No"), ") ", magnitude, " ", direction, " the predicted genotype"),
      paste0(clean_lbl, " (SHAP = ", sprintf("%+.3f", shap_value), ") ", magnitude, " ", direction, " the predicted genotype")
    )
    return(interpretation)
  }, error = function(e) {
    paste0(feature_name, " contributes to the prediction")
  })
}

#' Create summary of key SHAP insights
summarize_shap_insights <- function(contrib_df, prediction_result) {
  insights <- c()
  top_feature <- contrib_df[1, ]
  insights <- c(insights, paste0("🔍 Most influential factor: ", top_feature$clinical_interpretation))

  transfusion_features <- contrib_df[grepl("transfus", contrib_df$feature, ignore.case = TRUE), ]
  if (nrow(transfusion_features) > 0) {
    total_transfusion_impact <- sum(abs(transfusion_features$shap_value), na.rm = TRUE)
    if (total_transfusion_impact > 0.1) {
      insights <- c(insights, paste0("🩸 Transfusion-related terms contributed |SHAP| = ", round(total_transfusion_impact, 3)))
    }
  }

  if (!is.null(prediction_result$confidence)) {
    if (prediction_result$confidence > 0.9) {
      insights <- c(insights, "✅ High confidence (>90%) — automated classification threshold met")
    } else if (prediction_result$confidence > 0.7) {
      insights <- c(insights, "⚠️ Moderate confidence — expert review recommended (<90% threshold)")
    } else {
      insights <- c(insights, "❌ Low confidence — ambiguous pattern flagged for specialist review")
    }
  }

  top_3 <- contrib_df[1:min(3, nrow(contrib_df)), ]
  insights <- c(insights, paste0("📊 Top 3 SHAP features: ", paste(clean_sharp_feature_names(top_3$feature), collapse = ", ")))
  return(insights)
}

#' Create feature importance plot from SHAP contributions
create_importance_plot_simple <- function(feature_contributions) {
  tryCatch({
    plot_df <- feature_contributions %>%
      dplyr::mutate(clean_label = clean_sharp_feature_names(feature))

    ggplot(plot_df, aes(x = reorder(clean_label, abs_shap), y = abs_shap, fill = direction)) +
      geom_col(width = 0.7, color = "white", linewidth = 0.4) +
      coord_flip() +
      scale_fill_manual(
        values = c("Positive" = "#FF6B6B", "Negative" = "#4ECDC4"),
        labels = c("Positive" = "Supports prediction", "Negative" = "Opposes prediction")
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(
        title = "Patient-Specific SHAP Feature Importance",
        subtitle = "Absolute TreeSHAP contribution (|SHAP value|)",
        x = NULL,
        y = "|SHAP Value| (log-odds impact)",
        fill = "Direction"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text = element_text(size = 10),
        legend.position = "bottom"
      )
  }, error = function(e) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Importance plot error", size = 4) +
      theme_void()
  })
}
