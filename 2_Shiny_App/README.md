# Hemoglobin Genotype Prediction Shiny App (Enhanced 5-Class)

A Shiny web application for predicting hemoglobin genotypes using clinical decision rules and machine learning with SHAP explanations. **Now supports 5-class classification including Sickle Cell Trait (AS) with clinical constraints**.

## Features

- **Enhanced 5-class prediction system**: Now includes AS (Sickle Cell Trait) classification
- **Clinical constraints framework**: Prevents biochemically impossible predictions (e.g., SC without HgbC)
- **Two-tier prediction approach**: Clinical rules for high-confidence cases, machine learning for ambiguous cases
- **Transfusion awareness**: Accounts for recent transfusion status and its confounding effects
- **SHAP explanations**: Interpretable AI explanations for machine learning predictions
- **Real-time validation**: Input validation with helpful error messages
- **Professional UI**: Modern, responsive interface with Bootstrap styling
- **Downloadable reports**: Generate detailed prediction reports

## Installation

### Prerequisites

- R (version 4.0 or higher)
- Required R packages (install with the commands below)

### R Package Installation

```r
install.packages(c(
  "shiny",
  "shinydashboard", 
  "shinyWidgets",
  "DT",
  "plotly",
  "tidyverse",
  "tidymodels",
  "shapviz",
  "vip",
  "bslib",
  "dplyr",
  "conflicted"
))
```

### Setup

1. Clone or download this application directory
2. Ensure the trained model file `tidymodels_genotype_workflow_with_transfusion.rds` is in the `data/` folder
3. Open R and set your working directory to the app folder

## Usage

### Running the App

```r
# From R console in the app directory
shiny::runApp()

# Or specify the app directory path
shiny::runApp("/path/to/genotype_prediction_app")
```

### Using the Application

1. **Enter Patient Data**:
   - Required: HgbS%, HgbA%, HgbF%, HgbA2%, and transfusion status
   - Optional: CBC parameters and test quality information

2. **Get Predictions**:
   - Click "Predict Genotype" to run analysis
   - System will first try clinical rules, then ML if needed

3. **Interpret Results**:
   - View prediction confidence and rationale
   - For ML predictions, examine SHAP explanations
   - Review clinical guidelines for context

4. **Download Reports**:
   - Generate and download detailed prediction reports

## System Architecture

### Prediction Pipeline

1. **Input Validation**: Real-time validation of user inputs
2. **Enhanced Clinical Rules**: High-confidence pattern matching for 5 genotype classes:
   - SS (Sickle Cell Disease)
   - SC (Sickle-Hemoglobin C Disease)  
   - S-β⁰ (Sickle Beta-Zero Thalassemia)
   - S-β⁺ (Sickle Beta-Plus Thalassemia)
   - AS (Sickle Cell Trait)

3. **Enhanced Machine Learning**: XGBoost model for 5-class predictions with:
   - **Clinical constraints**: Biochemical possibility enforcement
   - Transfusion-aware features and interactions
   - SHAP explanations for interpretability
   - Confidence-based expert review flagging

### Key Features

- **Clinical Constraints**: Prevents impossible genotype-lab combinations (e.g., SC disease with HgbC = 0%)
- **Transfusion Handling**: Automatically accounts for confounding effects of recent transfusion
- **Quality Awareness**: Incorporates test quality information when available
- **Expert Review**: Flags low-confidence cases for manual review
- **Clinical Context**: Provides relevant clinical guidelines and interpretations

## Model Information

- **Algorithm**: Enhanced XGBoost via tidymodels framework with clinical constraints
- **Classes**: 5-class prediction (SS, SC, S-β⁰, S-β⁺, AS)
- **Training Data**: Enhanced dataset with AS patients
- **Performance**: ~97% cross-validation accuracy on 5-class problem with constraints
- **Features**: 25+ engineered features including transfusion interactions
- **Clinical Safety**: Post-processing constraints ensure biochemically plausible predictions
- **Special Handling**: Transfusion-aware rules for AS and S-β⁺ patterns

## File Structure

```
genotype_prediction_app/
├── app.R                    # Main Shiny application
├── utils/
│   ├── data_validation.R    # Input validation functions
│   ├── clinical_rules.R     # Clinical decision rules
│   ├── shap_generator.R     # SHAP explanation generation
│   └── helpers.R           # UI rendering and utility functions
├── data/
│   └── tidymodels_genotype_workflow_with_transfusion.rds
├── www/
│   └── styles.css          # Custom CSS styling
└── README.md               # This file
```

## Important Disclaimers

- **Research Use Only**: This tool is for research and educational purposes
- **Clinical Validation Required**: All predictions should be validated by clinical experts
- **Transfusion Effects**: Recent transfusion significantly affects interpretation
- **Atypical Cases**: Unusual patterns may require additional testing

## Technical Considerations

### Performance
- SHAP calculations may take 10-30 seconds for complex cases
- Model loading occurs once at startup
- Responsive design for mobile and desktop use

### Error Handling
- Comprehensive input validation
- Graceful handling of SHAP calculation failures
- Informative error messages for users

### Security
- No patient data persistence between sessions
- Input sanitization for all user inputs
- Session isolation to prevent data leakage

## Support

For technical issues or questions about the application:

1. Check input validation messages for common issues
2. Ensure all required fields are completed
3. Verify model file is properly located in `data/` folder
4. Check R console for detailed error messages

## License

This application is provided for research and educational use. Please ensure appropriate institutional review and compliance with data protection regulations when using with patient data. 