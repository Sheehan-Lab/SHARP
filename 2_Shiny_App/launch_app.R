#!/usr/bin/env Rscript

# Launch script for Hemoglobin Genotype Prediction Shiny App
# Run this script to start the application

cat("Starting Hemoglobin Genotype Prediction App...\n")

# Check if all required packages are installed
required_packages <- c(
  "shiny", "shinydashboard", "shinyWidgets", "DT", "plotly",
  "tidyverse", "tidymodels", "shapviz", "vip", "bslib", "dplyr", "conflicted"
)

missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}

# Load shiny
library(shiny)

# Check if model file exists
model_file <- "data/tidymodels_genotype_workflow_with_transfusion.rds"
if (!file.exists(model_file)) {
  stop("Model file not found: ", model_file, 
       "\nPlease ensure the trained model is in the data/ directory.")
}

cat("All dependencies satisfied. Launching app...\n")
cat("The app will open in your default web browser.\n")
cat("Press Ctrl+C to stop the app.\n\n")

# Run the app
runApp(launch.browser = TRUE) 