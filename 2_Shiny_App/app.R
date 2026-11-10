# ==============================================================================
# SHARP: Hemoglobin Genotype Prediction Shiny App (v2.0 - Single & Batch + SHAP)
# ==============================================================================
# Interactive web application for predicting hemoglobin genotypes using
# clinical rules and machine learning (Enhanced 5-Class XGBoost) with
# exact multi-class TreeSHAP explanations for both single patients and
# uploaded batch spreadsheets.
# ==============================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(tidyverse)
library(tidymodels)
library(xgboost)
library(shapviz)
library(vip)
library(bslib)

# Resolve namespace conflicts
library(conflicted)
conflicts_prefer(shiny::observe)
conflicts_prefer(shiny::dataTableOutput)
conflicts_prefer(shiny::renderDataTable)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)

# Source utility functions
source("utils/data_validation.R")
source("utils/clinical_rules.R")
source("utils/shap_generator.R")
source("utils/helpers.R")

# ==============================================================================
# USER INTERFACE
# ==============================================================================

ui <- page_sidebar(
  title = "SHARP 🩸→🧬",
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#3498db",
    secondary = "#95a5a6"
  ),

  tags$head(
    tags$style(HTML("
      /* Target the navbar container for the background color */
      .navbar {
        background-color: #d9534f !important; /* Bootstrap 'danger' red */
      }

      /* Target the title text 'SHARP' which is a 'navbar-brand' */
      .navbar-brand {
        color: white !important;
        font-size: 2.2rem !important;
        font-weight: 600 !important;
      }

      .batch-kpi-box {
        background: #f8f9fa;
        border-left: 4px solid #d9534f;
        border-radius: 6px;
        padding: 12px 15px;
        margin-bottom: 12px;
      }
      .batch-kpi-val {
        font-size: 1.5rem;
        font-weight: 700;
        color: #2c3e50;
      }
      .batch-kpi-label {
        font-size: 0.82rem;
        color: #6c757d;
        text-transform: uppercase;
        letter-spacing: 0.04em;
      }
    "))
  ),

  # Include custom CSS
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),

  # Sidebar with mode selector + inputs
  sidebar = sidebar(
    width = 400,

    # Mode Selector Card
    card(
      card_header("Analysis Mode", class = "bg-danger text-white"),
      radioButtons(
        "analysis_mode",
        label = NULL,
        choices = c(
          "Single Patient Prediction" = "single",
          "Batch Spreadsheet Analysis (CSV / Excel)" = "batch"
        ),
        selected = "single"
      )
    ),

    # ---------------- SINGLE PATIENT SIDEBAR CONTROLS ----------------
    conditionalPanel(
      condition = "input.analysis_mode == 'single'",

      # Patient Information
      card(
        card_header("Patient Information"),
        textInput("patient_id", "Patient ID (optional)",
                  placeholder = "Enter patient identifier")
      ),

      # Electrophoresis Data
      card(
        card_header(
          "Electrophoresis Results",
          class = "bg-danger text-white"
        ),
        helpText("Required fields marked with *"),

        numericInput("hgb_s_pct", "HgbS %*",
                     value = NULL, min = 0, max = 100, step = 0.1),
        numericInput("hgb_a_pct", "HgbA %*",
                     value = NULL, min = 0, max = 100, step = 0.1),
        numericInput("hgb_c_pct", "HgbC %",
                     value = 0, min = 0, max = 100, step = 0.1),
        numericInput("hgb_f_pct", "HgbF %*",
                     value = NULL, min = 0, max = 100, step = 0.1),
        numericInput("hgb_a2_pct", "HgbA2 %*",
                     value = NULL, min = 0, max = 100, step = 0.1),

        # Real-time validation feedback
        uiOutput("electrophoresis_validation")
      ),

      # Transfusion Status
      card(
        card_header(
          "Transfusion History",
          class = "bg-danger text-white"
        ),
        radioButtons("is_transfused", "Recent Transfusion (within 3 months)*",
                     choices = list("No" = 0, "Yes" = 1),
                     selected = character(0)),
        helpText("Critical for accurate interpretation")
      ),

      # CBC Data (Optional)
      card(
        card_header("CBC Parameters"),
        checkboxInput("include_cbc", "Include CBC data", TRUE),

        conditionalPanel(
          condition = "input.include_cbc",
          numericInput("total_hgb", "Total Hemoglobin (g/dL)",
                       value = NULL, min = 3, max = 20, step = 0.1),
          numericInput("mcv", "MCV (fL)",
                       value = NULL, min = 50, max = 120, step = 0.1),
          numericInput("mchc", "MCHC (g/dL)",
                       value = NULL, min = 25, max = 40, step = 0.1),
          numericInput("rdw", "RDW (%)",
                       value = NULL, min = 10, max = 30, step = 0.1),
          numericInput("hematocrit", "Hematocrit (%)",
                       value = NULL, min = 15, max = 60, step = 0.1),
          numericInput("rbc_count", "RBC Count (M/μL)",
                       value = NULL, min = 2, max = 8, step = 0.1)
        )
      ),

      # Test Quality
      card(
        card_header("Test Quality Assessment (Optional)"),
        selectInput("test_confidence_level", "Test Quality",
                    choices = list("Not specified" = "", "High" = "High",
                                   "Medium" = "Medium", "Low" = "Low")),
        numericInput("clarity_score", "Clarity Score (1-10)",
                     value = NULL, min = 1, max = 10, step = 1)
      ),

      # Action Buttons
      div(
        class = "d-grid gap-2",
        actionButton("predict_btn", "Predict Genotype",
                     class = "btn-primary btn-lg", icon = icon("dna")),
        actionButton("reset_btn", "Reset Form",
                     class = "btn-secondary"),
        br(),
        downloadButton("download_report", "Download Report",
                       class = "btn-outline-primary")
      )
    ),

    # ---------------- BATCH ANALYSIS SIDEBAR CONTROLS ----------------
    conditionalPanel(
      condition = "input.analysis_mode == 'batch'",
      card(
        card_header("1. Format Template", class = "bg-danger text-white"),
        p(class = "small text-muted mb-2",
          "Download the formatted template CSV containing required electrophoresis headers, optional CBC headers, and 7 sample patient profiles."),
        downloadButton("download_batch_template", "Download Template CSV",
                       class = "btn-outline-danger w-100")
      ),
      card(
        card_header("2. Upload Patient Spreadsheet", class = "bg-danger text-white"),
        fileInput(
          "batch_file",
          "Select CSV or Excel (.xlsx / .xls)",
          accept = c(".csv", ".xlsx", ".xls")
        ),
        helpText("Required columns: hgb_s_pct, hgb_a_pct, hgb_c_pct, hgb_f_pct, hgb_a2_pct, is_transfused. Optional: patient_id, total_hgb, mcv, mchc, rdw, hematocrit, rbc_count."),
        div(
          class = "d-grid gap-2 mt-2",
          actionButton("run_batch_btn", "Run Batch Analysis + SHAP",
                       class = "btn-primary btn-lg", icon = icon("gears")),
          actionButton("load_demo_batch_btn", "Load Demo Batch (7 Cases)",
                       class = "btn-outline-secondary btn-sm")
        )
      ),
      card(
        card_header("3. Export Results + SHAP"),
        p(class = "small text-muted mb-2",
          "Output spreadsheet includes all original columns, predicted genotypes, constrained 5-class probabilities, expert-review flags, and exact multi-class TreeSHAP columns."),
        downloadButton("download_batch_results", "Download Results CSV (with SHAP)",
                       class = "btn-success w-100")
      )
    )
  ),

  # Main content area
  layout_columns(

    # Left Column: Single Prediction Results OR Batch Summary & Format Guide
    card(
      full_screen = TRUE,
      card_header(
        uiOutput("left_card_title_ui"),
        class = "bg-danger text-white"
      ),
      uiOutput("left_panel_content_ui")
    ),

    # Right Column: Explanation & Interpretation + Batch Analysis Tab
    card(
      full_screen = TRUE,
      card_header("Explanation, Batch Results & Clinical Reference", class = "bg-danger text-white"),

      navset_card_tab(
        id = "sharp_main_tabs",
        nav_panel(
          "Explanation",
          value = "tab_explanation",
          uiOutput("shap_explanation")
        ),
        nav_panel(
          "Batch Analysis",
          value = "tab_batch",
          uiOutput("batch_results_right_ui")
        ),
        nav_panel(
          "Clinical Guidelines",
          value = "tab_guidelines",
          uiOutput("clinical_guidelines")
        ),
        nav_panel(
          "Feature Importance",
          value = "tab_importance",
          plotlyOutput("feature_importance_plot")
        ),
        nav_panel(
          "About",
          value = "tab_about",
          uiOutput("about_content")
        )
      )
    ),

    col_widths = c(6, 6)
  )
)

# ==============================================================================
# SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {

  # Load the latest trained 5-class model
  genotype_model <- reactive({
    readRDS("data/tidymodels_genotype_workflow_with_transfusion.rds")
  })

  # Reactive values to store single and batch results
  values <- reactiveValues(
    prediction_result = NULL,
    shap_data = NULL,
    input_data = NULL,
    validation_status = NULL,
    batch_results_df = NULL,
    batch_summary = NULL,
    batch_error = NULL
  )

  # Automatically switch tabs when user toggles analysis_mode
  observeEvent(input$analysis_mode, {
    if (identical(input$analysis_mode, "batch")) {
      nav_select("sharp_main_tabs", selected = "tab_batch", session = session)
    } else if (identical(input$analysis_mode, "single")) {
      nav_select("sharp_main_tabs", selected = "tab_explanation", session = session)
    }
  }, ignoreInit = TRUE)

  # Dynamic Left Card Header
  output$left_card_title_ui <- renderUI({
    if (identical(input$analysis_mode, "batch")) {
      span(icon("table"), " Batch Analysis Summary & Spreadsheet Specification")
    } else {
      span(icon("stethoscope"), " Prediction Results")
    }
  })

  # Real-time input validation
  observe({
    validation <- validate_inputs(input)
    values$validation_status <- validation
  })

  # Validation feedback UI
  output$electrophoresis_validation <- renderUI({
    if (!is.null(values$validation_status) && !values$validation_status$valid) {
      div(
        class = "alert alert-warning mt-2",
        HTML(paste("⚠️", paste(values$validation_status$errors, collapse = "<br>")))
      )
    }
  })

  # Main single-patient prediction logic
  observeEvent(input$predict_btn, {
    validation <- validate_inputs(input)
    if (!validation$valid) {
      showNotification(
        paste("Please fix the following errors:", paste(validation$errors, collapse = ", ")),
        type = "error",
        duration = 5
      )
      return()
    }

    withProgress(message = "Predicting genotype...", value = 0, {
      incProgress(0.1, detail = "Preparing input data")
      input_data <- prepare_input_data(input)
      values$input_data <- input_data

      incProgress(0.3, detail = "Applying Tier-1 clinical rules")
      rule_result <- apply_clinical_rules(input_data)

      incProgress(0.5, detail = "Running Tier-2 5-class XGBoost + biochemical constraints")
      ml_result <- predict_with_ml(input_data, genotype_model())

      target_cls <- if (isTRUE(rule_result$high_confidence) &&
                        rule_result$prediction %in% c("SS", "SC", "S_beta_0", "S_beta_plus", "AS")) {
        rule_result$prediction
      } else {
        ml_result$prediction
      }

      incProgress(0.7, detail = "Generating exact multi-class TreeSHAP explanation")
      shap_result <- safely_generate_shap(input_data, genotype_model(), target_class = target_cls)

      if (!is.null(shap_result$error)) {
        showNotification("SHAP analysis failed, but prediction is available", type = "warning")
      } else {
        values$shap_data <- shap_result$result
      }

      incProgress(0.9, detail = "Compiling results")
      final_result <- compile_final_result(rule_result, ml_result, input_data)
      # Always attach constrained ML probabilities for transparency
      if (is.null(final_result$probabilities) && !is.null(ml_result$probabilities)) {
        final_result$probabilities <- ml_result$probabilities
      }
      values$prediction_result <- final_result

      incProgress(1.0, detail = "Complete!")
    })

    nav_select("sharp_main_tabs", selected = "tab_explanation", session = session)
    showNotification("Prediction complete!", type = "message")
  })

  # Reset form
  observeEvent(input$reset_btn, {
    updateNumericInput(session, "hgb_s_pct", value = NA)
    updateNumericInput(session, "hgb_a_pct", value = NA)
    updateNumericInput(session, "hgb_c_pct", value = 0)
    updateNumericInput(session, "hgb_f_pct", value = NA)
    updateNumericInput(session, "hgb_a2_pct", value = NA)
    updateRadioButtons(session, "is_transfused", selected = character(0))
    updateCheckboxInput(session, "include_cbc", value = FALSE)
    updateTextInput(session, "patient_id", value = "")
    updateSelectInput(session, "test_confidence_level", selected = "")
    updateNumericInput(session, "clarity_score", value = NA)

    values$prediction_result <- NULL
    values$shap_data <- NULL
    values$input_data <- NULL

    showNotification("Form reset", type = "message")
  })

  # ---------------- BATCH ANALYSIS HANDLERS ----------------

  # Template download
  output$download_batch_template <- downloadHandler(
    filename = function() {
      "SHARP_Batch_Template.csv"
    },
    content = function(file) {
      readr::write_csv(get_sharp_batch_template(), file)
    }
  )

  # Helper to run batch analysis on a dataframe
  execute_sharp_batch <- function(raw_df) {
    withProgress(message = "Running SHARP Batch Analysis + TreeSHAP...", value = 0, {
      res <- run_sharp_batch_analysis(
        raw_df = raw_df,
        model = genotype_model(),
        progress_callback = function(val, det) {
          setProgress(value = val, detail = det)
        }
      )
      values$batch_results_df <- res$results_df
      values$batch_summary    <- res$summary_stats
      values$batch_error      <- res$error
    })

    if (!is.null(values$batch_error)) {
      showNotification(values$batch_error, type = "error", duration = 8)
    } else {
      nav_select("sharp_main_tabs", selected = "tab_batch", session = session)
      showNotification(
        sprintf("Successfully analyzed %d patients with full SHAP columns!", nrow(values$batch_results_df)),
        type = "message",
        duration = 5
      )
    }
  }

  # Run batch on uploaded file
  observeEvent(input$run_batch_btn, {
    req(input$batch_file)
    ext <- tolower(tools::file_ext(input$batch_file$name))
    raw_df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package 'readxl' is required to read Excel files. Please upload a .csv file.")
        }
        readxl::read_excel(input$batch_file$datapath)
      } else {
        readr::read_csv(input$batch_file$datapath, show_col_types = FALSE)
      }
    }, error = function(e) {
      showNotification(paste("Failed to read file:", e$message), type = "error", duration = 6)
      return(NULL)
    })
    req(raw_df)
    execute_sharp_batch(raw_df)
  })

  # Load demo batch
  observeEvent(input$load_demo_batch_btn, {
    demo_df <- get_sharp_batch_template()
    execute_sharp_batch(demo_df)
  })

  # Download batch results CSV
  output$download_batch_results <- downloadHandler(
    filename = function() {
      paste0("SHARP_Batch_Predictions_with_SHAP_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$batch_results_df)
      readr::write_csv(values$batch_results_df, file)
    }
  )

  # ---------------- UI RENDERING ----------------

  # Left Panel Content: switches between Single Patient Results and Batch Summary/Spec
  output$left_panel_content_ui <- renderUI({
    if (identical(input$analysis_mode, "batch")) {
      render_sharp_batch_left_panel()
    } else {
      render_prediction_results(values$prediction_result, values$input_data)
    }
  })

  # Render Batch Left Panel (KPI Summary + Format Specification Table)
  render_sharp_batch_left_panel <- function() {
    format_spec_df <- tibble::tibble(
      `Column Name` = c(
        "hgb_s_pct *", "hgb_a_pct *", "hgb_c_pct", "hgb_f_pct *", "hgb_a2_pct *", "is_transfused *",
        "patient_id", "total_hgb", "mcv", "mchc", "rdw", "hematocrit", "rbc_count", "clarity_score"
      ),
      `Requirement` = c(
        "Required", "Required", "Required (0 if absent)", "Required", "Required", "Required (0/1 or No/Yes)",
        "Optional", "Optional", "Optional", "Optional", "Optional", "Optional", "Optional", "Optional"
      ),
      `Description / Units` = c(
        "Hemoglobin S (%)", "Hemoglobin A (%)", "Hemoglobin C (%)", "Hemoglobin F (%)", "Hemoglobin A2 (%)",
        "Recent transfusion <90 days (0 = No, 1 = Yes)", "Patient or Sample ID",
        "Total Hemoglobin (g/dL)", "Mean Corpuscular Volume (fL)", "MCHC (g/dL)",
        "Red Cell Distribution Width (%)", "Hematocrit (%)", "RBC Count (M/μL)", "Electrophoresis Clarity (1-10)"
      )
    )

    tagList(
      if (!is.null(values$batch_error)) {
        div(class = "alert alert-danger", icon("triangle-exclamation"), " ", strong("Batch Error: "), values$batch_error)
      },

      if (!is.null(values$batch_summary)) {
        s <- values$batch_summary
        div(
          class = "mb-4",
          h5(class = "fw-bold mb-3", icon("chart-pie"), " Cohort Batch Analysis Summary"),
          div(
            class = "row g-2",
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", s$total_rows),
              div(class = "batch-kpi-label", "Total Patients")
            )),
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", paste0(s$coverage_pct, "%")),
              div(class = "batch-kpi-label", "Automated Coverage (>=90% Conf)")
            )),
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", s$expert_review_calls),
              div(class = "batch-kpi-label", "Flagged for Expert Review")
            )),
            div(class = "col-md-6", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", paste0(s$rule_classified, " / ", s$ml_classified)),
              div(class = "batch-kpi-label", "Tier-1 Rule Calls / Tier-2 ML Calls")
            )),
            div(class = "col-md-6", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", paste0(s$mean_confidence_pct, "%")),
              div(class = "batch-kpi-label", "Mean Prediction Confidence")
            ))
          ),
          div(
            class = "card mt-2",
            div(class = "card-header py-2", strong("Cohort Mean |SHAP| Feature Importance")),
            div(class = "card-body p-2", plotOutput("batch_shap_summary_plot", height = "280px"))
          )
        )
      } else {
        div(
          class = "alert alert-info mb-3",
          icon("circle-info"),
          strong(" Ready for Batch Analysis: "),
          "Upload a CSV or Excel file in the sidebar (or click 'Load Demo Batch') to generate predictions and full multi-class SHAP columns."
        )
      },

      div(
        class = "card",
        div(class = "card-header bg-light", strong(icon("list-check"), " Required & Optional Input Spreadsheet Columns")),
        div(
          class = "card-body p-2",
          p(class = "small text-muted mb-2",
            "Column names are matched case-insensitively and common clinical aliases (e.g., 'Hgb S %', 'HCT', 'Transfused') are automatically recognized."),
          DT::renderDataTable({
            datatable(
              format_spec_df,
              options = list(dom = "t", paging = FALSE, ordering = FALSE),
              rownames = FALSE
            )
          })
        )
      )
    )
  }

  # Render Batch Right Tab (Interactive Results + SHAP Table)
  output$batch_results_right_ui <- renderUI({
    if (is.null(values$batch_results_df)) {
      return(div(
        class = "alert alert-info",
        icon("table"),
        h4("Batch Spreadsheet Results & SHAP Output"),
        p("Switch to 'Batch Spreadsheet Analysis' in the sidebar and upload a CSV/Excel file (or click 'Load Demo Batch') to inspect row-by-row predictions, constrained 5-class probabilities, and per-feature SHAP contributions.")
      ))
    }

    div(
      div(
        class = "d-flex justify-content-between align-items-center mb-3",
        div(
          h4(class = "mb-0", "Batch Predictions & SHAP Transparency Table"),
          span(class = "text-muted small",
               sprintf("Showing %d rows and %d total columns (including grouped & raw feature TreeSHAP values)",
                       nrow(values$batch_results_df), ncol(values$batch_results_df)))
        ),
        downloadButton("download_batch_results_inline", "Download Full CSV (with SHAP)",
                       class = "btn-success btn-sm")
      ),
      DT::dataTableOutput("batch_results_table")
    )
  })

  output$download_batch_results_inline <- downloadHandler(
    filename = function() {
      paste0("SHARP_Batch_Predictions_with_SHAP_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$batch_results_df)
      readr::write_csv(values$batch_results_df, file)
    }
  )

  output$batch_results_table <- DT::renderDataTable({
    req(values$batch_results_df)
    datatable(
      values$batch_results_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })

  output$batch_shap_summary_plot <- renderPlot({
    req(values$batch_results_df)
    df_res <- values$batch_results_df
    group_cols <- c(
      "shap_hgb_s_group"       = "HbS Group",
      "shap_hgb_a_group"       = "HbA Group",
      "shap_hgb_c_group"       = "HbC Group",
      "shap_hgb_f_group"       = "HbF Group",
      "shap_hgb_a2_group"      = "HbA2 Group",
      "shap_transfusion_group" = "Transfusion Terms",
      "shap_mcv_group"         = "MCV Terms",
      "shap_total_hgb"         = "Total Hemoglobin",
      "shap_rdw"               = "RDW",
      "shap_rbc_count"         = "RBC Count",
      "shap_mchc"              = "MCHC"
    )
    present_cols <- intersect(names(group_cols), names(df_res))
    req(length(present_cols) > 0)

    imp_tbl <- tibble::tibble(
      Feature = unname(group_cols[present_cols]),
      MeanAbsSHAP = sapply(present_cols, function(col) mean(abs(as.numeric(df_res[[col]])), na.rm = TRUE))
    ) %>%
      dplyr::arrange(MeanAbsSHAP) %>%
      dplyr::mutate(Feature = factor(Feature, levels = Feature))

    ggplot(imp_tbl, aes(x = Feature, y = MeanAbsSHAP)) +
      geom_col(fill = "#d9534f", width = 0.7) +
      coord_flip() +
      labs(
        x = NULL,
        y = "Mean |SHAP| Contribution (log-odds)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(face = "bold"),
        panel.grid.major.y = element_blank()
      )
  })

  # Render SHAP explanation
  output$shap_explanation <- renderUI({
    render_shap_explanation(values$shap_data, values$prediction_result)
  })

  # SHAP waterfall plot
  output$shap_waterfall <- renderPlot({
    if (!is.null(values$shap_data) && !is.null(values$shap_data$waterfall_plot)) {
      values$shap_data$waterfall_plot
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No SHAP data available") +
        theme_void()
    }
  })

  # SHAP contributions table
  output$shap_contributions_table <- DT::renderDataTable({
    if (!is.null(values$shap_data) && !is.null(values$shap_data$feature_contributions)) {
      contrib_display <- values$shap_data$feature_contributions %>%
        mutate(
          Feature_Display = if ("feature_label_clean" %in% names(.)) feature_label_clean else feature
        ) %>%
        select(
          Feature = Feature_Display,
          `SHAP Value` = shap_value,
          Direction = direction,
          `Clinical Interpretation` = clinical_interpretation
        ) %>%
        mutate(
          `SHAP Value` = round(`SHAP Value`, 4)
        )

      datatable(
        contrib_display,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          columnDefs = list(
            list(width = "200px", targets = 3)
          )
        ),
        rownames = FALSE
      )
    } else {
      datatable(data.frame(Message = "No SHAP contributions available"))
    }
  })

  # Feature importance plot
  output$feature_importance_plot <- renderPlotly({
    if (!is.null(values$shap_data) && !is.null(values$shap_data$importance_plot)) {
      ggplotly(values$shap_data$importance_plot, tooltip = c("x", "y"))
    } else if (!is.null(values$input_data) && !is.null(values$prediction_result)) {
      fallback_plot <- create_simple_feature_plot(values$input_data)
      ggplotly(fallback_plot, tooltip = c("x", "y"))
    } else {
      plotly_empty() %>%
        layout(
          title = "Feature importance will appear here after prediction",
          annotations = list(
            text = "Run a prediction to see feature importance",
            showarrow = FALSE,
            font = list(size = 16)
          )
        )
    }
  })

  # Clinical guidelines
  output$clinical_guidelines <- renderUI({
    render_clinical_guidelines()
  })

  # About content
  output$about_content <- renderUI({
    render_about_content()
  })

  # Download report
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("genotype_prediction_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      generate_prediction_report(values$prediction_result, values$shap_data,
                                 values$input_data, file)
    }
  )
}

# ==============================================================================
# RUN APPLICATION
# ==============================================================================

shinyApp(ui = ui, server = server)
