# ===============================================================================
# Process AS Data and Update Pipeline for 5-Class All-Cases Training
# Author: Enhanced ML Pipeline  
# Date: January 2025
# ===============================================================================
# This script:
# 1. Processes AS (sickle trait) data from CSV to match existing format
# 2. Integrates AS data with main dataset 
# 3. Updates pipeline to train on ALL cases (not just ambiguous)
# 4. Implements 5-class classification: SS, SC, S_beta_0, S_beta_plus, AS
# ===============================================================================

library(tidyverse)
library(readr)
library(janitor)

# Set options for cleaner output
options(scipen = 999, digits = 3)

cat("🚀 STARTING ENHANCED 5-CLASS ALL-CASES ML PIPELINE\n")
cat("===============================================================================\n")

# ===============================================================================
# 1. PROCESS AS (SICKLE TRAIT) DATA
# ===============================================================================

cat("\n=== STEP 1: PROCESSING AS (SICKLE TRAIT) DATA ===\n")

# Load AS data
cat("Loading AS data from Sickle Trait Patients.csv...\n")
as_raw <- read_csv("Sickle Trait Patients.csv", col_types = cols(.default = "c"))

cat("AS raw data dimensions:", nrow(as_raw), "rows x", ncol(as_raw), "columns\n")

# Clean column names
as_raw <- as_raw %>% clean_names()
print("AS data columns:")
print(colnames(as_raw))

# Process AS data to match existing format
process_as_data <- function(df) {
  cat("Processing AS data to match main dataset format...\n")
  
  # Convert to long format and extract lab values
  as_long <- df %>%
    rename(
      test_id = order_id,
      patient_mrn = patient_mrn,
      component_name = lab_results_order_component_name,
      result_value = lab_results_order_result_value
    ) %>%
    filter(!is.na(result_value)) %>%
    mutate(
      result_value_numeric = suppressWarnings(as.numeric(result_value)),
      # Map component names to our standard format
      component_clean = case_when(
        component_name == "HGB S QUANT" ~ "hgb_s_pct",
        component_name == "HGB A QUANT" ~ "hgb_a_pct", 
        component_name == "HGB A2 QUANT" ~ "hgb_a2_pct",
        component_name == "HGB F QUANT" ~ "hgb_f_pct",
        component_name == "HGB C QUANT" ~ "hgb_c_pct",
        component_name == "HEMOGLOBIN-BLOOD" ~ "total_hgb",
        component_name == "HEMATOCRIT-BLOOD" ~ "hematocrit",
        component_name == "RBC COUNT-BLOOD" ~ "rbc_count",
        component_name == "MCV" ~ "mcv",
        component_name == "MCHC" ~ "mchc",
        component_name == "RDW" ~ "rdw",
        component_name == "HGB PHENOTYPE" ~ "phenotype",
        TRUE ~ str_to_lower(component_name)
      ),
      genotype = "AS"  # All are AS/sickle trait
    ) %>%
    filter(!is.na(result_value_numeric), !is.na(component_clean))
  
  # Score test quality (similar to main processing)
  test_scores <- as_long %>%
    filter(component_clean != "phenotype") %>%
    group_by(patient_mrn, test_id) %>%
    summarise(
      n_components = n(),
      has_hgb_s = any(component_clean == "hgb_s_pct"),
      has_hgb_a = any(component_clean == "hgb_a_pct"),
      has_hgb_c = any(component_clean == "hgb_c_pct"),
      has_hgb_f = any(component_clean == "hgb_f_pct"),
      has_hgb_a2 = any(component_clean == "hgb_a2_pct"),
      has_mcv = any(component_clean == "mcv"),
      completeness_score = has_hgb_s + has_hgb_a + has_hgb_c + has_hgb_f + has_hgb_a2 + has_mcv,
      hgb_s_val = if_else(has_hgb_s, result_value_numeric[component_clean == "hgb_s_pct"][1], NA_real_),
      hgb_a_val = if_else(has_hgb_a, result_value_numeric[component_clean == "hgb_a_pct"][1], 0),
      .groups = "drop"
    ) %>%
    mutate(
      clarity_score = case_when(
        # Clear AS pattern: HbS 35-45%, HbA 55-65%
        !is.na(hgb_s_val) & hgb_s_val >= 30 & hgb_s_val <= 50 & 
        hgb_a_val >= 50 & hgb_a_val <= 70 ~ 10,
        !is.na(hgb_s_val) & hgb_s_val >= 25 & hgb_s_val <= 55 ~ 8,
        !is.na(hgb_s_val) & hgb_s_val > 10 ~ 6,
        TRUE ~ 5
      ),
      total_quality_score = completeness_score * 2 + clarity_score
    )
  
  # Select best test per patient
  best_as_tests <- test_scores %>%
    group_by(patient_mrn) %>%
    slice_max(total_quality_score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Get lab values for best tests
  best_as_values <- as_long %>%
    inner_join(
      best_as_tests %>% select(patient_mrn, test_id),
      by = c("patient_mrn", "test_id")
    ) %>%
    filter(component_clean != "phenotype") %>%
    group_by(patient_mrn, component_clean) %>%
    summarise(value = median(result_value_numeric, na.rm = TRUE), .groups = "drop")
  
  # Pivot to wide format
  as_wide <- best_as_values %>%
    pivot_wider(
      names_from = component_clean,
      values_from = value,
      values_fill = NA
    ) %>%
    mutate(genotype = "AS", genotype_simplified = "AS")
  
  # Add quality scores
  as_final <- as_wide %>%
    left_join(
      best_as_tests %>% select(patient_mrn, completeness_score, clarity_score, total_quality_score),
      by = "patient_mrn"
    ) %>%
    mutate(
      # Ensure all expected columns exist with defaults first
      hgb_s_pct = if_else(is.na(hgb_s_pct), 0, hgb_s_pct),
      hgb_a_pct = if_else(is.na(hgb_a_pct), 0, hgb_a_pct),
      hgb_a2_pct = if_else(is.na(hgb_a2_pct), 3.2, hgb_a2_pct),
      hgb_f_pct = if_else(is.na(hgb_f_pct), 0, hgb_f_pct)
    ) %>%
    # Add missing columns that might not exist in AS data
    {
      if(!"hgb_c_pct" %in% colnames(.)) {
        mutate(., hgb_c_pct = 0)
      } else {
        mutate(., hgb_c_pct = if_else(is.na(hgb_c_pct), 0, hgb_c_pct))
      }
    } %>%
    {
      if(!"total_hgb" %in% colnames(.)) {
        mutate(., total_hgb = NA_real_)
      } else {
        .
      }
    } %>%
    {
      if(!"hematocrit" %in% colnames(.)) {
        mutate(., hematocrit = NA_real_)
      } else {
        .
      }
    } %>%
    {
      if(!"rbc_count" %in% colnames(.)) {
        mutate(., rbc_count = NA_real_)
      } else {
        .
      }
    } %>%
    {
      if(!"mcv" %in% colnames(.)) {
        mutate(., mcv = NA_real_)
      } else {
        .
      }
    } %>%
    {
      if(!"mchc" %in% colnames(.)) {
        mutate(., mchc = NA_real_)
      } else {
        .
      }
    } %>%
    {
      if(!"rdw" %in% colnames(.)) {
        mutate(., rdw = NA_real_)
      } else {
        .
      }
    } %>%
    mutate(
      # Create derived features
      hgb_a_present = if_else(hgb_a_pct > 0, 1, 0),
      hgb_s_present = if_else(hgb_s_pct > 0, 1, 0),
      hgb_c_present = if_else(hgb_c_pct > 0, 1, 0),
      hgb_f_elevated = if_else(hgb_f_pct > 8.6, 1, 0),
      hgb_a2_elevated = if_else(hgb_a2_pct > 4.0, 1, 0),
      hgb_s_to_total = hgb_s_pct / 100,
      hgb_a_to_s_ratio = if_else(hgb_s_pct > 0, hgb_a_pct / hgb_s_pct, NA_real_),
      total_variant_pct = hgb_s_pct + hgb_c_pct,
      
      # Quality indicators
      test_confidence_level = case_when(
        completeness_score >= 5 & clarity_score >= 8 ~ "High",
        completeness_score >= 4 & clarity_score >= 6 ~ "Medium", 
        TRUE ~ "Low"
      ),
      high_quality_test = (completeness_score >= 5 & clarity_score >= 8),
      had_hba_below_1_ever = FALSE,  # AS patients shouldn't have HbA < 1%
      phenotype_interpretation = "Hgb AS"
    )
  
  return(as_final)
}

# Process the AS data
as_processed <- process_as_data(as_raw)

cat("AS data processing results:\n")
cat("- Raw AS records:", nrow(as_raw), "\n")
cat("- Processed AS patients:", nrow(as_processed), "\n")
cat("- Average completeness score:", round(mean(as_processed$completeness_score), 1), "\n")
cat("- Average clarity score:", round(mean(as_processed$clarity_score), 1), "\n")

# Quality distribution of AS data
cat("\nAS data quality distribution:\n")
print(as_processed %>% count(test_confidence_level, sort = TRUE))

# ===============================================================================
# 2. LOAD AND COMBINE WITH MAIN DATASET 
# ===============================================================================

cat("\n=== STEP 2: LOADING AND COMBINING WITH MAIN DATASET ===\n")

# Load main processed data
if (file.exists("processed_hemoglobin_data_improved.csv")) {
  main_data <- read_csv("processed_hemoglobin_data_improved.csv", col_types = cols(patient_mrn = "c"))
  cat("Loaded main dataset:", nrow(main_data), "patients\n")
} else {
  cat("Main dataset not found. Please run data_preparation_v2.R first.\n")
  stop("processed_hemoglobin_data_improved.csv not found")
}

# Check for overlapping patients
overlap_patients <- intersect(main_data$patient_mrn, as_processed$patient_mrn)
cat("Overlapping patients between main and AS datasets:", length(overlap_patients), "\n")

if (length(overlap_patients) > 0) {
  cat("Removing", length(overlap_patients), "overlapping patients from AS data to avoid duplicates\n")
  as_processed <- as_processed %>%
    filter(!patient_mrn %in% overlap_patients)
}

# Check column compatibility
main_cols <- colnames(main_data)
as_cols <- colnames(as_processed) 

missing_in_as <- setdiff(main_cols, as_cols)
missing_in_main <- setdiff(as_cols, main_cols)

cat("Columns missing in AS data:", length(missing_in_as), "\n")
cat("Columns missing in main data:", length(missing_in_main), "\n")

# Add missing columns to AS data with appropriate defaults
for (col in missing_in_as) {
  as_processed[[col]] <- NA
}

# Add missing columns to main data with appropriate defaults  
for (col in missing_in_main) {
  main_data[[col]] <- NA
}

# Combine datasets
combined_data <- bind_rows(main_data, as_processed) %>%
  mutate(patient_mrn = as.character(patient_mrn))

cat("Combined dataset:\n")
cat("- Main data patients:", nrow(main_data), "\n")
cat("- Additional AS patients:", nrow(as_processed), "\n") 
cat("- Total combined patients:", nrow(combined_data), "\n")

# Check final genotype distribution
cat("\nFinal genotype distribution:\n")
final_genotype_dist <- combined_data %>% count(genotype_simplified, sort = TRUE)
print(final_genotype_dist)

# ===============================================================================
# 3. SAVE ENHANCED DATASET
# ===============================================================================

cat("\n=== STEP 3: SAVING ENHANCED DATASET ===\n")

write_csv(combined_data, "enhanced_hemoglobin_data_with_as.csv")
cat("✅ Enhanced dataset saved to: enhanced_hemoglobin_data_with_as.csv\n")

# Summary report
enhancement_report <- list(
  original_patients = nrow(main_data),
  additional_as_patients = nrow(as_processed),
  total_patients = nrow(combined_data),
  genotype_distribution = combined_data %>% count(genotype_simplified, sort = TRUE),
  quality_distribution = combined_data %>% count(test_confidence_level, sort = TRUE)
)

writeLines(capture.output(enhancement_report), "dataset_enhancement_report.txt")
cat("📊 Enhancement report saved to: dataset_enhancement_report.txt\n")

cat("\n🎯 STEP 1 COMPLETE: Enhanced dataset ready for 5-class all-cases training!\n")
cat("Next: Run the updated tidymodels pipeline with 5-class classification\n") 