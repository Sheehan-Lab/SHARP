# ===============================================================================
# Hemoglobin Genotype Prediction - Data Preparation Script (Updated for New-raw-query.csv)
# Author: Jonathan Wade
# Date: May 29 2025
# Updated: June 4 2025
# ===============================================================================
# This script processes raw lab data to create a clean dataset for ML modeling
# Now uses New-raw-query.csv as the source data, with updated column mapping
# ===============================================================================

library(tidyverse)
library(readr)
library(janitor)

# Set options for cleaner output
options(scipen = 999, digits = 3)

# ===============================================================================
# 1. LOAD AND EXAMINE RAW DATA (NEW FILE & COLUMN MAPPING)
# ===============================================================================

cat("Loading raw training data from New-raw-query.csv...\n")
raw_data <- read_csv("New-raw-query.csv", col_types = cols(.default = "c"))

cat("Raw data dimensions:", nrow(raw_data), "rows x", ncol(raw_data), "columns\n")

# Clean column names for easier handling
raw_data <- raw_data %>% clean_names()

# Ensure patient_mrn is character from the start
raw_data <- raw_data %>% mutate(patient_mrn = as.character(patient_mrn))

# Print column names to debug naming issues
print(colnames(raw_data))

# Use lab_orders_patient_mrn_order_id as the unique test identifier throughout
# (replaces serial_number everywhere)

# Map new column names to expected names for downstream code
raw_data <- raw_data %>%
  rename(
    patient_mrn = patient_mrn,  # Already matches
    test_id = lab_orders_patient_mrn_order_id,  # Use as unique test identifier
    lab_results_order_component_name = lab_results_order_component_name,  # Already matches
    lab_results_order_result_value = lab_results_order_result_value,      # Already matches
    genotype = samples_subject_genotype # Use the new genotype column as primary source
  )

# Create a numeric version of lab_results_order_result_value for quantitative fields
raw_data <- raw_data %>%
  mutate(lab_results_order_result_value_numeric = suppressWarnings(as.numeric(lab_results_order_result_value)))

# The genotype column is now directly available; no need to extract from HGB PHENOTYPE rows
# Remove the previous logic for extracting genotype from phenotype rows

# After initial load
cat('After initial load: Rows =', nrow(raw_data), ' Unique patient_mrn =', length(unique(raw_data$patient_mrn)), ' Missing patient_mrn =', sum(is.na(raw_data$patient_mrn)), '\n')

# ===============================================================================
# 2. FIRST FIX: CLEAN UP GENOTYPE TYPOS BEFORE PROCESSING
# ===============================================================================

cat("\n=== FIXING GENOTYPE TYPOS ===\n")

# Check original genotype distribution
cat("Original genotype distribution (focusing on Sbeta variants):\n")
original_genotypes <- raw_data %>% 
  count(genotype, sort = TRUE) %>%
  filter(str_detect(genotype, "Sbet|Sbeta|SBeta"))

print(original_genotypes)

# Fix genotype typos and consolidate variants
raw_data <- raw_data %>%
  mutate(
    genotype_cleaned = case_when(
      # Fix Sbet0 typo → Sbeta0
      genotype == "Sbet0" ~ "Sbeta0",
      # Standardize all Sbeta0 variants
      genotype %in% c("SBeta0", "Sbeta 0", "SBet0/HPFH") ~ "Sbeta0",
      # Standardize Sbeta+ variants  
      genotype %in% c("\"Sbeta\"") ~ "Sbeta+",
      # Keep everything else as is
      TRUE ~ genotype
    )
  ) %>%
  # Replace original genotype column
  select(-genotype) %>%
  rename(genotype = genotype_cleaned)

# Check corrected genotype distribution
cat("\nCorrected genotype distribution (Sbeta variants):\n")
corrected_genotypes <- raw_data %>% 
  count(genotype, sort = TRUE) %>%
  filter(str_detect(genotype, "Sbet|Sbeta|SBeta"))

print(corrected_genotypes)

# ===============================================================================
# 3. EXPLORE GENOTYPE AND LAB COMPONENT STRUCTURE
# ===============================================================================

cat("\n=== GENOTYPE DISTRIBUTION ===\n")
genotype_counts <- raw_data %>% 
  count(genotype, sort = TRUE) %>%
  filter(!is.na(genotype) & genotype != "" & !str_detect(genotype, "^\\d+$"))

print(genotype_counts)

cat("\n=== LAB COMPONENTS AVAILABLE ===\n")
lab_components <- raw_data %>% 
  count(lab_results_order_component_name, sort = TRUE) %>%
  filter(!is.na(lab_results_order_component_name))

print(lab_components)

# ===============================================================================
# 4. IDENTIFY AND CHARACTERIZE TRANSFUSION CASES (KEEP FOR ML TRAINING)
# ===============================================================================

cat("\n=== IDENTIFYING TRANSFUSION STATUS FOR ML FEATURE ===\n")

# Extract HGB PHENOTYPE interpretations
phenotype_data <- raw_data %>%
  filter(lab_results_order_component_name == "HGB PHENOTYPE") %>%
  select(test_id, patient_mrn, genotype, phenotype_text = lab_results_order_result_value) %>%
  distinct()

# Define transfusion indicators
transfusion_patterns <- c(
  "Patient Recently Transfused",
  "Suggestive of Recent Transfusion",
  "Recent Transfusion",
  "Transfused",
  "transfus"  # catch any case variations
)

# Create transfusion indicator
phenotype_data <- phenotype_data %>%
  mutate(
    is_transfused = str_detect(phenotype_text, 
                              regex(paste(transfusion_patterns, collapse = "|"), 
                                   ignore_case = TRUE)),
    clean_phenotype = str_trim(phenotype_text)
  )

# Summary of transfusion cases
transfusion_summary <- phenotype_data %>%
  count(is_transfused, clean_phenotype, sort = TRUE)

cat("Transfusion cases found:\n")
print(transfusion_summary %>% filter(is_transfused))

# Create transfusion lookup for including in ML features
transfusion_lookup <- phenotype_data %>%
  select(patient_mrn, is_transfused) %>%
  distinct() %>%
  # Handle patients with multiple tests - if ANY test shows transfusion, mark as transfused
  group_by(patient_mrn) %>%
  summarise(is_transfused = any(is_transfused, na.rm = TRUE), .groups = "drop")

cat("\nNumber of transfused patients identified for ML feature:", sum(transfusion_lookup$is_transfused), "\n")
cat("Number of non-transfused patients:", sum(!transfusion_lookup$is_transfused), "\n")

# ===============================================================================
# 5. EXTRACT QUANTITATIVE ELECTROPHORESIS VALUES
# ===============================================================================

cat("\n=== EXTRACTING ELECTROPHORESIS VALUES ===\n")

# Define the key electrophoresis components we want
key_components <- c(
  "HGB S QUANT",
  "HGB A QUANT", 
  "HGB A2 QUANT",
  "HGB F QUANT",
  "HGB C QUANT",
  "HEMOGLOBIN-BLOOD",    # Total Hgb
  "HEMATOCRIT-BLOOD",    # Hct
  "RBC COUNT-BLOOD",     # RBC count
  "MCV",                 # Mean corpuscular volume
  "MCHC",                # Mean corpuscular Hgb concentration
  "RDW"                  # Red cell distribution width
)

# Extract quantitative lab values (NOW INCLUDING TRANSFUSED PATIENTS)
quant_data <- raw_data %>%
  filter(
    lab_results_order_component_name %in% key_components,
    !is.na(lab_results_order_result_value_numeric)
  ) %>%
  select(
    test_id,
    patient_mrn, 
    genotype,
    component = lab_results_order_component_name,
    value = lab_results_order_result_value_numeric
  ) %>%
  mutate(
    value = as.numeric(value),
    # Clean component names for easier handling
    component_clean = case_when(
      component == "HGB S QUANT" ~ "hgb_s_pct",
      component == "HGB A QUANT" ~ "hgb_a_pct", 
      component == "HGB A2 QUANT" ~ "hgb_a2_pct",
      component == "HGB F QUANT" ~ "hgb_f_pct",
      component == "HGB C QUANT" ~ "hgb_c_pct",
      component == "HEMOGLOBIN-BLOOD" ~ "total_hgb",
      component == "HEMATOCRIT-BLOOD" ~ "hematocrit",
      component == "RBC COUNT-BLOOD" ~ "rbc_count",
      component == "MCV" ~ "mcv",
      component == "MCHC" ~ "mchc",
      component == "RDW" ~ "rdw",
      TRUE ~ str_to_lower(component)
    )
  )

# After including all patients (after quant_data creation)
cat('After including all patients: Rows =', nrow(quant_data), ' Unique patient_mrn =', length(unique(quant_data$patient_mrn)), ' Missing patient_mrn =', sum(is.na(quant_data$patient_mrn)), '\n')

# Check value ranges for each component
cat("Value ranges by component:\n")
value_summary <- quant_data %>%
  group_by(component_clean) %>%
  summarise(
    n = n(),
    min_val = min(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    max_val = max(value, na.rm = TRUE),
    .groups = "drop"
  )

print(value_summary)

# Identify patients who ever had HBA < 1
hba_low_patients <- quant_data %>%
  filter(component_clean == "hgb_a_pct", value < 1) %>%
  distinct(patient_mrn) %>%
  mutate(had_hba_below_1_ever = TRUE)

# ===============================================================================
# 6. INTELLIGENT TEST SELECTION (ONE BEST TEST PER PATIENT)
# ===============================================================================

cat("\n=== INTELLIGENT TEST SELECTION PER PATIENT ===\n")

# Create test quality scoring system
score_test_quality <- function(df) {
  df %>%
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
      hgb_s_val = if_else(has_hgb_s, value[component_clean == "hgb_s_pct"][1], NA_real_),
      hgb_a_val = if_else(has_hgb_a, value[component_clean == "hgb_a_pct"][1], 0),
      hgb_c_val = if_else(has_hgb_c, value[component_clean == "hgb_c_pct"][1], 0),
      hgb_a2_val = if_else(has_hgb_a2, value[component_clean == "hgb_a2_pct"][1], NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      clarity_score = case_when(
        !is.na(hgb_s_val) & hgb_s_val > 85 & hgb_a_val == 0 & hgb_c_val == 0 ~ 10,
        !is.na(hgb_s_val) & !is.na(hgb_c_val) & hgb_s_val >= 40 & hgb_s_val <= 60 & hgb_c_val >= 40 & hgb_c_val <= 60 ~ 10,
        !is.na(hgb_a2_val) & hgb_a2_val > 3.5 & hgb_a_val > 0 ~ 9,
        !is.na(hgb_s_val) & hgb_s_val > 75 & hgb_a_val == 0 ~ 8,
        !is.na(hgb_s_val) & hgb_s_val > 50 ~ 6,
        TRUE ~ 5
      ),
      total_quality_score = completeness_score * 2 + clarity_score,
      recency_score = 0.001, # No recency info, keep as default
      final_score = total_quality_score + recency_score
    )
}

# Score all tests
test_scores <- score_test_quality(quant_data)

cat("Test quality scoring summary:\n")
score_summary <- test_scores %>%
  summarise(
    total_tests = n(),
    avg_completeness = round(mean(completeness_score), 1),
    avg_clarity = round(mean(clarity_score), 1),
    avg_total_quality = round(mean(total_quality_score), 1)
  )
print(score_summary)

# Select best test per patient
best_tests_per_patient <- test_scores %>%
  group_by(patient_mrn) %>%
  slice_max(final_score, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Selected best tests:\n")
cat("- Patients with multiple tests:", 
    length(unique(test_scores$patient_mrn)) - nrow(best_tests_per_patient), "\n")
cat("- Patients after selection:", nrow(best_tests_per_patient), "\n")

# After intelligent test selection (after best_tests_per_patient creation)
cat('After test selection: Rows =', nrow(best_tests_per_patient), ' Unique patient_mrn =', length(unique(best_tests_per_patient$patient_mrn)), ' Missing patient_mrn =', sum(is.na(best_tests_per_patient$patient_mrn)), '\n')

# Get the actual lab values for the selected best tests
selected_test_data <- quant_data %>%
  inner_join(
    best_tests_per_patient %>% select(patient_mrn, test_id),
    by = c("patient_mrn", "test_id")
  )

cat("Lab measurements from selected tests:", nrow(selected_test_data), "\n")

# Aggregate duplicate measurements within the same test before pivoting
patient_values <- selected_test_data %>%
  group_by(patient_mrn, genotype, component_clean) %>%
  summarise(
    value = median(value, na.rm = TRUE),  # Use median if multiple values within same test
    .groups = "drop"
  )

cat("Aggregated measurements:", nrow(patient_values), "\n")

# Pivot to wide format
wide_data <- patient_values %>%
  pivot_wider(
    names_from = component_clean,
    values_from = value,
    values_fill = NA
  )

cat("Wide format dimensions (best tests only):", nrow(wide_data), "rows x", ncol(wide_data), "columns\n")

# After pivot_wider (after wide_data creation)
cat('After pivot_wider: Rows =', nrow(wide_data), ' Unique patient_mrn =', length(unique(wide_data$patient_mrn)), ' Missing patient_mrn =', sum(is.na(wide_data$patient_mrn)), '\n')

# Add test quality metadata to the final dataset
wide_data <- wide_data %>%
  left_join(
    best_tests_per_patient %>% 
      select(patient_mrn, completeness_score, clarity_score, total_quality_score),
    by = "patient_mrn"
  )

# After merging phenotype_clean (after wide_data <- wide_data %>% left_join(...))
cat('After merging phenotype_clean: Rows =', nrow(wide_data), ' Unique patient_mrn =', length(unique(wide_data$patient_mrn)), ' Missing patient_mrn =', sum(is.na(wide_data$patient_mrn)), '\n')

# ===============================================================================
# 7. ADD CLINICAL PHENOTYPE DATA
# ===============================================================================

# Add the clean phenotype interpretations (INCLUDING TRANSFUSED PATIENTS)
phenotype_clean <- phenotype_data %>%
  select(patient_mrn, phenotype_interpretation = clean_phenotype) %>%
  distinct()

# Merge with wide data and add transfusion status
wide_data <- wide_data %>%
  left_join(phenotype_clean, by = "patient_mrn") %>%
  left_join(transfusion_lookup, by = "patient_mrn") %>%
  mutate(
    # Default to not transfused if no phenotype data available
    is_transfused = if_else(is.na(is_transfused), FALSE, is_transfused),
    # Convert to numeric for ML (0 = not transfused, 1 = transfused)
    is_transfused_numeric = as.numeric(is_transfused)
  )

# ===============================================================================
# 8. QUALITY-BASED DATA FILTERING AND VALIDATION
# ===============================================================================

cat("\n=== QUALITY-BASED DATA FILTERING (INCLUDING TRANSFUSED PATIENTS) ===\n")

# Check for patients with valid genotypes (updated to include fixed variants)
valid_genotypes <- c("SS", "SC", "AS", "AA", "Sbeta0", "Sbeta+", "CC", "AC", "SF", 
                    "SS/HPFH", "SS- Lepore")

# Apply more stringent quality filters but keep transfused patients
quality_filtered_patients <- wide_data %>%
  filter(
    # Must have a valid genotype
    genotype %in% valid_genotypes | str_detect(genotype, "^(SS|SC|AS|AA|S)"),
    # Filter out obvious data quality issues for hgb_s_pct if it's present, allow NAs to pass for now.
    (is.na(hgb_s_pct) | (hgb_s_pct >= 0 & hgb_s_pct <= 100)),
    # Require minimum diagnostic clarity
    clarity_score >= 5
  )

cat("Quality filtering results (including transfused patients):\n")
cat("- Before quality filtering:", nrow(wide_data), "patients\n")
cat("- After quality filtering:", nrow(quality_filtered_patients), "patients\n")
cat("- Patients removed due to low quality:", nrow(wide_data) - nrow(quality_filtered_patients), "\n")

# Show transfusion distribution in final dataset
cat("\nTransfusion status in final dataset:\n")
transfusion_final_summary <- quality_filtered_patients %>% 
  count(is_transfused, name = "count") %>%
  mutate(percentage = round(count / sum(count) * 100, 1))
print(transfusion_final_summary)

# Quality distribution analysis
quality_analysis <- quality_filtered_patients %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    avg_completeness = round(mean(completeness_score, na.rm = TRUE), 1),
    avg_clarity = round(mean(clarity_score, na.rm = TRUE), 1),
    avg_total_quality = round(mean(total_quality_score, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

cat("\nQuality analysis by genotype:\n")
print(quality_analysis)

# Update the complete_patients to use quality_filtered_patients
complete_patients <- quality_filtered_patients

cat("Final patients with high-quality data:", nrow(complete_patients), "\n")

# Check genotype distribution in quality-filtered data
cat("\nGenotype distribution in quality-filtered data:\n")
clean_genotype_dist <- complete_patients %>% 
  count(genotype, sort = TRUE)
print(clean_genotype_dist)

# ===============================================================================
# 8A. INCORPORATE MANUALLY IDENTIFIED SBETA CASES
# ===============================================================================

cat("\n=== INCORPORATING MANUALLY IDENTIFIED SBETA CASES ===\n")

# Load manual cases
manual_cases <- read_csv("Manually Identified Sbeta Cases.csv", 
                        col_types = cols(patient_mrn = "c", .default = "c")) %>%
  # Convert numeric columns
  mutate(across(c(hematocrit, hgb_a2_pct, hgb_c_pct, hgb_f_pct, hgb_s_pct,
                  mchc, mcv, rbc_count, rdw, total_hgb, hgb_a_pct), 
                ~as.numeric(.))) %>%
  # Remove duplicates
  distinct(patient_mrn, .keep_all = TRUE)

cat("Loaded manual cases:", nrow(manual_cases), "unique patients\n")

# Score test quality for manual cases
manual_cases_scored <- manual_cases %>%
  mutate(
    # Count complete measurements
    has_hgb_s = !is.na(hgb_s_pct),
    has_hgb_a = !is.na(hgb_a_pct),
    has_hgb_c = !is.na(hgb_c_pct),
    has_hgb_f = !is.na(hgb_f_pct),
    has_hgb_a2 = !is.na(hgb_a2_pct),
    has_mcv = !is.na(mcv),
    
    # Calculate completeness score
    completeness_score = has_hgb_s + has_hgb_a + has_hgb_c + 
                        has_hgb_f + has_hgb_a2 + has_mcv,
    
    # Calculate clarity score based on clinical patterns
    clarity_score = case_when(
      # Clear SS pattern
      !is.na(hgb_s_pct) & hgb_s_pct > 85 & 
      (is.na(hgb_a_pct) | hgb_a_pct == 0) & 
      (is.na(hgb_c_pct) | hgb_c_pct == 0) ~ 10,
      
      # Clear SC pattern
      !is.na(hgb_s_pct) & !is.na(hgb_c_pct) & 
      hgb_s_pct >= 40 & hgb_s_pct <= 60 & 
      hgb_c_pct >= 40 & hgb_c_pct <= 60 ~ 10,
      
      # Clear beta-thal pattern
      !is.na(hgb_a2_pct) & hgb_a2_pct > 3.5 & 
      !is.na(hgb_a_pct) & hgb_a_pct > 0 ~ 9,
      
      # Moderately clear SS
      !is.na(hgb_s_pct) & hgb_s_pct > 75 & 
      (is.na(hgb_a_pct) | hgb_a_pct == 0) ~ 8,
      
      # Some clarity: high HbS
      !is.na(hgb_s_pct) & hgb_s_pct > 50 ~ 6,
      
      # Default: moderate clarity
      TRUE ~ 5
    ),
    
    # Calculate total quality score
    total_quality_score = completeness_score * 2 + clarity_score,
    
    # Add confidence level
    test_confidence_level = case_when(
      completeness_score >= 5 & clarity_score >= 8 ~ "High",
      completeness_score >= 4 & clarity_score >= 6 ~ "Medium", 
      TRUE ~ "Low"
    ),
    
    # Map genotypes to our standardized categories
    genotype = case_when(
      str_detect(`Actual Genotype`, "Sbeta\\+") ~ "Sbeta+",
      str_detect(`Actual Genotype`, "Sbeta0") ~ "Sbeta0",
      str_detect(`Actual Genotype`, "SS/") ~ "SS/HPFH",
      str_detect(`Actual Genotype`, "^SS") ~ "SS",
      str_detect(`Actual Genotype`, "AS") ~ "AS",
      TRUE ~ `Actual Genotype`
    ),
    
    # Add high quality indicator
    high_quality_test = (completeness_score >= 5 & clarity_score >= 8)
  )

cat("\nManual cases quality distribution:\n")
print(manual_cases_scored %>% count(test_confidence_level, sort = TRUE))

cat("\nManual cases genotype distribution:\n")
print(manual_cases_scored %>% count(genotype, sort = TRUE))

# Prepare manual cases for merging
manual_cases_final <- manual_cases_scored %>%
  select(patient_mrn, genotype, phenotype_interpretation,
         matches("hgb_.*"), matches("total_.*"), 
         hematocrit, mcv, mchc, rbc_count, rdw,
         completeness_score, clarity_score, total_quality_score,
         test_confidence_level, high_quality_test)

# Merge with existing data, preferring manual cases where there's overlap
quality_filtered_patients <- bind_rows(
  # Keep existing patients not in manual cases
  quality_filtered_patients %>%
    filter(!patient_mrn %in% manual_cases_final$patient_mrn),
  # Add manual cases
  manual_cases_final
)

cat("\nUpdated dataset after incorporating manual cases:\n")
cat("Total patients:", nrow(quality_filtered_patients), "\n")

# After incorporating manual cases (after complete_patients creation)
cat('After incorporating manual cases: Rows =', nrow(complete_patients), ' Unique patient_mrn =', length(unique(complete_patients$patient_mrn)), ' Missing patient_mrn =', sum(is.na(complete_patients$patient_mrn)), '\n')

# Update complete_patients reference
complete_patients <- quality_filtered_patients

# ===============================================================================
# 9. CREATE DERIVED FEATURES (WITH CORRECTED GENOTYPE CONSOLIDATION)
# ===============================================================================

cat("\n=== CREATING DERIVED FEATURES ===\n")

final_data <- complete_patients %>%
  left_join(hba_low_patients, by = "patient_mrn") %>%
  mutate(
    # Presence/absence indicators
    hgb_a_present = if_else(!is.na(hgb_a_pct) & hgb_a_pct > 0, 1, 0, missing = 0),
    hgb_s_present = if_else(!is.na(hgb_s_pct) & hgb_s_pct > 0, 1, 0, missing = 0),
    hgb_c_present = if_else(!is.na(hgb_c_pct) & hgb_c_pct > 0, 1, 0, missing = 0),
    
    # Data-driven thresholds based on your dataset
    hgb_f_elevated = if_else(!is.na(hgb_f_pct) & hgb_f_pct > 8.6, 1, 0, missing = 0),  # Median from your data
    hgb_a2_elevated = if_else(!is.na(hgb_a2_pct) & hgb_a2_pct > 4.0, 1, 0, missing = 0),  # Based on your β-thal patients
    
    # Ratios and derived measures
    hgb_s_to_total = if_else(!is.na(hgb_s_pct), hgb_s_pct / 100, NA_real_),
    hgb_a_to_s_ratio = if_else(!is.na(hgb_a_pct) & !is.na(hgb_s_pct) & hgb_s_pct > 0, 
                               hgb_a_pct / hgb_s_pct, NA_real_),
    
    # Total variant hemoglobin percentage
    total_variant_pct = rowSums(select(., hgb_s_pct, hgb_c_pct), na.rm = TRUE),
    
    # Replace missing percentages with 0 for key hemoglobins
    hgb_s_pct = if_else(is.na(hgb_s_pct), 0, hgb_s_pct),
    hgb_a_pct = if_else(is.na(hgb_a_pct), 0, hgb_a_pct),
    hgb_c_pct = if_else(is.na(hgb_c_pct), 0, hgb_c_pct),
    hgb_f_pct = if_else(is.na(hgb_f_pct), 0, hgb_f_pct),
    hgb_a2_pct = if_else(is.na(hgb_a2_pct), 3.2, hgb_a2_pct),  # Normal A2 ~3.2%
    
    # Add transfusion status feature for ML
    is_transfused = if_else(is.na(is_transfused), FALSE, is_transfused),
    
    # FIXED: Create simplified genotype categories with proper Sbeta0 consolidation
    genotype_simplified = case_when(
      str_detect(genotype, "^SS") ~ "SS",  # Catches SS, SS/HPFH, SS- Lepore
      str_detect(genotype, "^SC") ~ "SC", 
      str_detect(genotype, "^AS") ~ "AS",
      str_detect(genotype, "^AA") ~ "AA",
      # FIXED: Consolidated pattern to catch all Sbeta0 variants (including fixed Sbet0)
      str_detect(genotype, "^Sbeta0$") ~ "S_beta_0",
      str_detect(genotype, "^Sbeta\\+$") ~ "S_beta_plus",
      str_detect(genotype, "^CC") ~ "CC",
      str_detect(genotype, "^AC") ~ "AC",
      TRUE ~ "Other"
    ),
    
    # NEW: Add quality-based confidence indicator
    high_quality_test = (completeness_score >= 5 & clarity_score >= 8),
    test_confidence_level = case_when(
      completeness_score >= 5 & clarity_score >= 8 ~ "High",
      completeness_score >= 4 & clarity_score >= 6 ~ "Medium", 
      TRUE ~ "Low"
    ),
    had_hba_below_1_ever = if_else(is.na(had_hba_below_1_ever), FALSE, had_hba_below_1_ever)
  )

# After creating final_data (after final_data creation)
cat('After creating final_data: Rows =', nrow(final_data), ' Unique patient_mrn =', length(unique(final_data$patient_mrn)), ' Missing patient_mrn =', sum(is.na(final_data$patient_mrn)), '\n')

# ===============================================================================
# 10. FINAL VALIDATION AND SUMMARY
# ===============================================================================

cat("\n=== FINAL DATASET SUMMARY (WITH QUALITY IMPROVEMENTS) ===\n")
cat("Final dataset dimensions:", nrow(final_data), "rows x", ncol(final_data), "columns\n")

cat("\nCORRECTED Simplified genotype distribution:\n")
final_genotype_summary <- final_data %>% count(genotype_simplified, sort = TRUE)
print(final_genotype_summary)

# Test quality distribution summary
cat("\nTest quality distribution:\n")
quality_dist <- final_data %>% count(test_confidence_level, sort = TRUE)
print(quality_dist)

# Cross-check: show original vs simplified mapping
cat("\nOriginal → Simplified genotype mapping:\n")
genotype_mapping <- final_data %>%
  count(genotype, genotype_simplified, sort = TRUE) %>%
  filter(genotype_simplified %in% c("S_beta_0", "S_beta_plus"))

print(genotype_mapping)

cat("\nKey electrophoresis value summaries:\n")
electro_summary <- final_data %>%
  select(hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct) %>%
  summary()
print(electro_summary)

cat("\nQuality score summaries:\n")
quality_summary <- final_data %>%
  select(completeness_score, clarity_score, total_quality_score) %>%
  summary()
print(quality_summary)

# ===============================================================================
# 11. SAVE PROCESSED DATA (UPDATED)
# ===============================================================================

cat("\n=== SAVING IMPROVED PROCESSED DATA ===\n")

# Save the list of unique patient_mrn from final_data for debugging
expected_mrns <- unique(final_data$patient_mrn)
cat(sprintf("Saving %d unique expected MRNs to expected_mrns.rds\n", length(expected_mrns)))
saveRDS(expected_mrns, "expected_mrns.rds")

# Before writing the final CSV, ensure patient_mrn is character (it should be by now)
final_data <- final_data %>% mutate(patient_mrn = as.character(patient_mrn))
write_csv(final_data, "processed_hemoglobin_data_improved.csv")
cat("IMPROVED processed data saved to: processed_hemoglobin_data_improved.csv\n")

# Save transfused patients list for reference (NOW INCLUDED IN DATASET)
transfused_in_final <- final_data %>% filter(is_transfused) %>% select(patient_mrn)
write_csv(transfused_in_final, "included_transfused_patients.csv")
cat("Included transfused patients saved to: included_transfused_patients.csv\n")

# Save quality analysis report
quality_report <- list(
  raw_records = nrow(raw_data),
  transfused_patients_included = sum(final_data$is_transfused, na.rm = TRUE),
  non_transfused_patients = sum(!final_data$is_transfused, na.rm = TRUE),
  tests_before_selection = nrow(test_scores),
  unique_patients = length(unique(test_scores$patient_mrn)), 
  final_patients = nrow(final_data),
  avg_test_quality = round(mean(final_data$total_quality_score), 1),
  high_quality_tests = sum(final_data$high_quality_test),
  genotype_distribution = final_data %>% count(genotype_simplified, sort = TRUE),
  quality_distribution = final_data %>% count(test_confidence_level, sort = TRUE),
  transfusion_distribution = final_data %>% count(is_transfused, sort = TRUE)
)

writeLines(capture.output(quality_report), "data_processing_quality_report.txt")
cat("Quality analysis report saved to: data_processing_quality_report.txt\n")

cat("\n✅ IMPROVED DATA PREPARATION COMPLETE!\n")
cat("🔹 IMPROVEMENTS MADE:\n")
cat("   • Intelligent test selection per patient MRN\n")
cat("   • Quality scoring based on completeness and diagnostic clarity\n")
cat("   • Preference for definitive over ambiguous results\n")
cat("   • Enhanced data quality filtering\n")
cat("   • TRANSFUSION STATUS INCLUDED AS ML FEATURE (not excluded)\n")
cat("Sbet0 typos have been consolidated with Sbeta0.\n")
cat("Ready for enhanced machine learning pipeline with transfusion awareness.\n") 