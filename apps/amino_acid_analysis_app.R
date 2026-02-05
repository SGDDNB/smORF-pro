## Shiny wizard for Amino Acid Analysis
##
## This app guides the user step-by-step:
## 0) Upload Annotation CSV
## 1) Generate FASTA
## 2) DeepTMHMM
## 3) TargetP 2.0
## 4) Deeploc 2.0
## 5) Interproscan
## 6) Integration & results

library(shiny)
library(DT)

# Source files relative to apps/ directory
source("../config/amino_acid_config.R", local = TRUE)
source("../scripts/process_amino_acid_results.R", local = TRUE)

steps <- c(
  "Welcome & Input",
  "FASTA generation",
  "DeepTMHMM",
  "TargetP 2.0",
  "Deeploc 2.0",
  "Interproscan",
  "Integration & Results"
)

ui <- fluidPage(
  titlePanel("smORF-pro: Amino Acid Analysis Wizard"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Workflow progress"),
      uiOutput("progress_ui"),
      hr(),
      actionButton("prev_step", "Previous"),
      actionButton("next_step", "Next", class = "btn-primary"),
      br(), br(),
      helpText("You can only proceed once the current step is complete.")
    ),
    mainPanel(
      width = 9,
      uiOutput("step_ui")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    step = 0L,
    step_complete = rep(FALSE, length(steps)),
    annotation = NULL,
    fasta_path = NULL,
    deeptmhmm_path = NULL,
    targetp_path = NULL,
    deeploc_path = NULL,
    interpro_path = NULL,
    files_saved = FALSE
  )

  observeEvent(input$prev_step, {
    if (rv$step > 0L) {
      rv$step <- rv$step - 1L
    }
  })

  observeEvent(input$next_step, {
    if (rv$step < length(steps) && rv$step_complete[rv$step + 1L]) {
      rv$step <- rv$step + 1L
    } else {
      showNotification("Please complete the current step before continuing.", type = "warning")
    }
  })

  output$progress_ui <- renderUI({
    current <- rv$step
    tags$ol(
      lapply(seq_along(steps), function(i) {
        style <- if (i - 1L == current) {
          "font-weight:bold;"
        } else if (rv$step_complete[i]) {
          "color:darkgreen;"
        } else {
          ""
        }
        tags$li(style = style, steps[i])
      })
    )
  })

  output$step_ui <- renderUI({
    switch(
      as.character(rv$step),
      "0" = step0_ui(),
      "1" = step1_ui(rv),
      "2" = step2_ui(rv),
      "3" = step3_ui(rv),
      "4" = step4_ui(rv),
      "5" = step5_ui(rv),
      "6" = step6_ui(rv),
      div("Unknown step")
    )
  })

  observeEvent(input$upload_annotation, {
    req(input$annotation_file)
    path <- input$annotation_file$datapath
    df <- tryCatch(
      {
        read.csv(path, stringsAsFactors = FALSE)
      },
      error = function(e) {
        showNotification(paste("Could not read CSV:", e$message), type = "error")
        NULL
      }
    )
    req(df)
    required_cols <- c("ORF_id", "gene_id", "gene_name", "iORF_type",
                       "gene_biotype", "len", "Peptide.seq", "starts")
    missing <- setdiff(required_cols, colnames(df))
    if (length(missing) > 0) {
      showNotification(
        paste("Annotation file is missing required columns:", paste(missing, collapse = ", ")),
        type = "error"
      )
      return()
    }

    annotation <- data.frame(
      iORFID      = df$ORF_id,
      Gene_id     = df$gene_id,
      Gene_name   = df$gene_name,
      ORF_type    = df$iORF_type,
      Gene_type   = df$gene_biotype,
      Source      = "sORFs.org",
      Length      = df$len / 3,
      Peptide_seq = df$Peptide.seq,
      Start_codon = df$starts,
      stringsAsFactors = FALSE
    )

    rv$annotation <- annotation
    rv$step_complete[1] <- TRUE
    showNotification("Annotation file loaded successfully.", type = "message")
  })

  observeEvent(input$generate_fasta, {
    req(rv$annotation)
    fasta_rel <- amino_acid_config$paths$fasta_output
    # Adjust path to be relative to apps/ directory
    fasta_path_from_apps <- file.path("..", fasta_rel)
    fasta_abs <- normalizePath(fasta_path_from_apps, winslash = "/", mustWork = FALSE)
    dir.create(dirname(fasta_abs), recursive = TRUE, showWarnings = FALSE)

    withProgress(message = "Generating FASTA file...", value = 0.5, {
      create_fasta_from_annotation(rv$annotation, fasta_abs)
    })
    rv$fasta_path <- fasta_abs
    rv$step_complete[2] <- TRUE
    showNotification(paste("FASTA generated at", fasta_abs), type = "message")
  })

  observeEvent(input$deeptmhmm_file, {
    req(input$deeptmhmm_file)
    path <- input$deeptmhmm_file$datapath
    rv$deeptmhmm_path <- path
    try({
      updated <- process_deeptmhmm(path, rv$annotation)
      rv$annotation <- updated
      rv$step_complete[3] <- TRUE
      showNotification("DeepTMHMM file processed successfully.", type = "message")
    }, silent = FALSE)
  })

  observeEvent(input$targetp_file, {
    req(input$targetp_file)
    path <- input$targetp_file$datapath
    rv$targetp_path <- path
    try({
      updated <- process_targetp(path, rv$annotation)
      rv$annotation <- updated
      rv$step_complete[4] <- TRUE
      showNotification("TargetP file processed successfully.", type = "message")
    }, silent = FALSE)
  })

  observeEvent(input$deeploc_file, {
    req(input$deeploc_file)
    path <- input$deeploc_file$datapath
    rv$deeploc_path <- path
    try({
      updated <- process_deeploc(path, rv$annotation)
      rv$annotation <- updated
      rv$step_complete[5] <- TRUE
      showNotification("Deeploc file processed successfully.", type = "message")
    }, silent = FALSE)
  })

  observeEvent(input$interpro_file, {
    req(input$interpro_file)
    path <- input$interpro_file$datapath
    rv$interpro_path <- path
    try({
      updated <- process_interproscan(path, rv$annotation)
      rv$annotation <- updated
      rv$step_complete[6] <- TRUE
      showNotification("Interproscan file processed successfully.", type = "message")
    }, silent = FALSE)
  })

  observe({
    if (!is.null(rv$annotation)) {
      rv$step_complete[7] <- TRUE
    }
  })


  # Save all files to project folder
  observeEvent(input$save_project, {
    req(rv$annotation, input$project_name)
    
    if (input$project_name == "" || nchar(trimws(input$project_name)) == 0) {
      showNotification("Please enter a project name.", type = "warning")
      return()
    }
    
    # Sanitize project name (remove invalid characters for file paths)
    project_name <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$project_name))
    if (project_name == "") {
      showNotification("Project name contains only invalid characters. Please use letters, numbers, hyphens, or underscores.", type = "error")
      return()
    }
    
    # Create projects folder if it doesn't exist (relative to apps/, so ../projects)
    projects_base <- file.path("..", "projects")
    project_path <- file.path(projects_base, project_name)
    
    # Convert to absolute path for operations
    project_path_abs <- normalizePath(project_path, mustWork = FALSE)
    
    tryCatch({
      withProgress(message = "Saving files to project folder...", value = 0, {
        # Create folder structure
        setProgress(0.1, detail = "Creating folder structure...")
        dir.create(file.path(project_path_abs, "data"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis"), 
                   recursive = TRUE, showWarnings = FALSE)
        
        # Save Annotation.rds
        setProgress(0.2, detail = "Saving Annotation.rds...")
        saveRDS(rv$annotation, file.path(project_path_abs, "data", "Annotation.rds"))
        
        # Create Annotation_DT.rds (with factors for DT display)
        setProgress(0.3, detail = "Creating Annotation_DT.rds...")
        annotation_dt <- rv$annotation
        # Apply factors to appropriate columns (as in README)
        factor_cols <- c("ORF_type", "Gene_type", "Start_codon")
        # Add tool result columns if they exist
        tool_cols <- c("deepTMHMM", "TargetP", "Deeploc", "interproscan")
        for (col in c(factor_cols, tool_cols)) {
          if (col %in% colnames(annotation_dt) && !all(is.na(annotation_dt[[col]]))) {
            annotation_dt[[col]] <- as.factor(annotation_dt[[col]])
          }
        }
        saveRDS(annotation_dt, file.path(project_path_abs, "data", "Annotation_DT.rds"))
        
        # Save Annotation.csv
        setProgress(0.4, detail = "Saving Annotation.csv...")
        utils::write.csv(rv$annotation, 
                         file.path(project_path_abs, "data", "Annotation.csv"), 
                         row.names = FALSE)
        
        # Copy FASTA file
        setProgress(0.5, detail = "Copying FASTA file...")
        if (!is.null(rv$fasta_path) && file.exists(rv$fasta_path)) {
          fasta_dest <- file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis", "smORFs.fa")
          file.copy(rv$fasta_path, fasta_dest, overwrite = TRUE)
        }
        
        # Copy tool output files with standardized names
        setProgress(0.6, detail = "Copying tool output files...")
        
        if (!is.null(rv$deeptmhmm_path) && file.exists(rv$deeptmhmm_path)) {
          file.copy(rv$deeptmhmm_path,
                    file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis", 
                             "predicted_topologies.3line"),
                    overwrite = TRUE)
        }
        
        if (!is.null(rv$targetp_path) && file.exists(rv$targetp_path)) {
          file.copy(rv$targetp_path,
                    file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis", 
                             "output_protein_type.txt"),
                    overwrite = TRUE)
        }
        
        if (!is.null(rv$deeploc_path) && file.exists(rv$deeploc_path)) {
          file.copy(rv$deeploc_path,
                    file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis", 
                             "Deeploc_output.csv"),
                    overwrite = TRUE)
        }
        
        if (!is.null(rv$interpro_path) && file.exists(rv$interpro_path)) {
          file.copy(rv$interpro_path,
                    file.path(project_path_abs, "data_preparation", "Amino_Acid_Analysis", 
                             "interproscan_output.tsv"),
                    overwrite = TRUE)
        }
        
        # Create metadata file
        setProgress(0.9, detail = "Creating metadata file...")
        metadata <- list(
          project_name = project_name,
          project_path = project_path_abs,
          created_date = as.character(Sys.Date()),
          created_time = as.character(Sys.time()),
          tools_processed = list(
            deeptmhmm = !is.null(rv$deeptmhmm_path) && file.exists(rv$deeptmhmm_path),
            targetp = !is.null(rv$targetp_path) && file.exists(rv$targetp_path),
            deeploc = !is.null(rv$deeploc_path) && file.exists(rv$deeploc_path),
            interproscan = !is.null(rv$interpro_path) && file.exists(rv$interpro_path)
          ),
          num_smorfs = nrow(rv$annotation),
          annotation_columns = colnames(rv$annotation)
        )
        saveRDS(metadata, file.path(project_path_abs, "data", "amino_acid_metadata.rds"))
      
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$files_saved <- TRUE
      showNotification(
        paste("All files saved successfully to:", project_path_abs),
        type = "message",
        duration = 10
      )
    }, error = function(e) {
      showNotification(
        paste("Error saving files:", e$message),
        type = "error",
        duration = 10
      )
      rv$files_saved <- FALSE
    })
  })


  # Output to check which tools were processed
  output$deeptmhmm_processed <- reactive({
    !is.null(rv$deeptmhmm_path) && file.exists(rv$deeptmhmm_path)
  })
  outputOptions(output, "deeptmhmm_processed", suspendWhenHidden = FALSE)

  output$targetp_processed <- reactive({
    !is.null(rv$targetp_path) && file.exists(rv$targetp_path)
  })
  outputOptions(output, "targetp_processed", suspendWhenHidden = FALSE)

  output$deeploc_processed <- reactive({
    !is.null(rv$deeploc_path) && file.exists(rv$deeploc_path)
  })
  outputOptions(output, "deeploc_processed", suspendWhenHidden = FALSE)

  output$interpro_processed <- reactive({
    !is.null(rv$interpro_path) && file.exists(rv$interpro_path)
  })
  outputOptions(output, "interpro_processed", suspendWhenHidden = FALSE)

  output$annotation_preview <- DT::renderDataTable({
    req(rv$annotation)
    DT::datatable(rv$annotation, options = list(pageLength = 10))
  })

  output$download_annotation_csv <- downloadHandler(
    filename = function() {
      "Annotation_amino_acid_processed.csv"
    },
    content = function(file) {
      req(rv$annotation)
      utils::write.csv(rv$annotation, file, row.names = FALSE)
    }
  )

  output$download_annotation_rds <- downloadHandler(
    filename = function() {
      "Annotation_amino_acid_processed.rds"
    },
    content = function(file) {
      req(rv$annotation)
      saveRDS(rv$annotation, file)
    }
  )

  # Output to check if files were saved
  output$files_saved <- reactive({
    rv$files_saved
  })
  outputOptions(output, "files_saved", suspendWhenHidden = FALSE)
  
  # Display saved path message
  output$saved_path_message <- renderText({
    if (rv$files_saved && !is.null(input$project_name) && input$project_name != "") {
      project_name <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$project_name))
      saved_path <- file.path("projects", project_name)
      paste("Files saved to:", saved_path)
    } else {
      ""
    }
  })
}

step0_ui <- function() {
  tagList(
    h3("Step 0: Welcome & Input Annotation"),
    p("This wizard will guide you through the amino acid analysis workflow used in smORF-pro."),
    p("You will:"),
    tags$ul(
      tags$li("Upload your smORF annotation CSV"),
      tags$li("Let the app create a FASTA file"),
      tags$li("Run external web tools (DeepTMHMM, TargetP, Deeploc, Interproscan)"),
      tags$li("Upload their result files and let the app process them"),
      tags$li("Download an updated Annotation table with all results integrated")
    ),
    fileInput("annotation_file", "Upload smORF annotation CSV", accept = c(".csv")),
    actionButton("upload_annotation", "Load annotation"),
    br(), br(),
    DTOutput("annotation_preview")
  )
}

step1_ui <- function(rv) {
  fasta_rel <- amino_acid_config$paths$fasta_output
  tagList(
    h3("Step 1: Generate FASTA"),
    p("We will create a FASTA file from your annotation for use with external web tools."),
    p("Output path (relative to project root):"),
    code(fasta_rel),
    br(), br(),
    actionButton("generate_fasta", "Generate FASTA file")
  )
}

step2_ui <- function(rv) {
  fasta_rel <- amino_acid_config$paths$fasta_output
  tagList(
    h3("Step 2: DeepTMHMM Analysis"),
    h4("What you need to do externally"),
    tags$ol(
      tags$li("Open DeepTMHMM in your browser."),
      tags$li("Upload the FASTA file created in Step 1."),
      tags$li("Run the analysis and wait for it to complete."),
      tags$li("Download the results in 3line format."),
      tags$li("Return here and upload the downloaded 3line file below.")
    ),
    p("DeepTMHMM website:"),
    tags$a(href = amino_acid_config$urls$deeptmhmm, target = "_blank", amino_acid_config$urls$deeptmhmm),
    br(),
    p("FASTA file path (relative): ", code(fasta_rel)),
    hr(),
    fileInput("deeptmhmm_file", "Upload DeepTMHMM 3line file", accept = c(".3line", ".txt"))
  )
}

step3_ui <- function(rv) {
  fasta_rel <- amino_acid_config$paths$fasta_output
  tagList(
    h3("Step 3: TargetP 2.0 Analysis"),
    h4("What you need to do externally"),
    tags$ol(
      tags$li("Open TargetP 2.0 in your browser."),
      tags$li("Upload the FASTA file created in Step 1."),
      tags$li("Select 'Non-plant' and 'Short output (no figures)'."),
      tags$li("Submit the job and wait for results."),
      tags$li("Download the 'Prediction summary' file."),
      tags$li("Return here and upload the summary file below.")
    ),
    p("TargetP 2.0 website:"),
    tags$a(href = amino_acid_config$urls$targetp, target = "_blank", amino_acid_config$urls$targetp),
    br(),
    p("FASTA file path (relative): ", code(fasta_rel)),
    hr(),
    fileInput("targetp_file", "Upload TargetP prediction summary (.txt)", accept = c(".txt"))
  )
}

step4_ui <- function(rv) {
  fasta_rel <- amino_acid_config$paths$fasta_output
  tagList(
    h3("Step 4: Deeploc 2.0 Analysis"),
    h4("What you need to do externally"),
    tags$ol(
      tags$li("Open Deeploc 2.0 in your browser."),
      tags$li("Upload the FASTA file created in Step 1."),
      tags$li("Select 'High-quality (Slow)' and 'Short output (no figures)'."),
      tags$li("Submit the job and wait for results."),
      tags$li("Download the 'CSV summary' file."),
      tags$li("Return here and upload the CSV summary file below.")
    ),
    p("Deeploc 2.0 website:"),
    tags$a(href = amino_acid_config$urls$deeploc, target = "_blank", amino_acid_config$urls$deeploc),
    br(),
    p("FASTA file path (relative): ", code(fasta_rel)),
    hr(),
    fileInput("deeploc_file", "Upload Deeploc CSV summary", accept = c(".csv"))
  )
}

step5_ui <- function(rv) {
  fasta_rel <- amino_acid_config$paths$fasta_output
  tagList(
    h3("Step 5: Interproscan Analysis"),
    h4("What you need to do externally"),
    tags$ol(
      tags$li("Open Interproscan sequence search in your browser."),
      tags$li("Upload the FASTA file created in Step 1."),
      tags$li("Click 'Search' to submit the job."),
      tags$li("Wait for completion."),
      tags$li("Use 'Group Actions' → 'Download TSV output'."),
      tags$li("Return here and upload the TSV output file below.")
    ),
    p("Interproscan website:"),
    tags$a(href = amino_acid_config$urls$interpro, target = "_blank", amino_acid_config$urls$interpro),
    br(),
    p("FASTA file path (relative): ", code(fasta_rel)),
    hr(),
    fileInput("interpro_file", "Upload Interproscan TSV file", accept = c(".tsv", ".txt"))
  )
}

step6_ui <- function(rv) {
  tagList(
    h3("Step 6: Integration & Results"),
    p("All available results have been integrated into the Annotation table below."),
    
    h4("Save to Project Folder"),
    p("Enter a project name. All files will be saved to projects/[project_name]/ with standardized names, ready for use in the main Shiny app:"),
    
    textInput("project_name", "Project name:", 
              placeholder = "e.g., my_smorf_analysis",
              width = "100%"),
    helpText("A folder will be created at: projects/[project_name]/"),
    br(),
    
    conditionalPanel(
      condition = "input.project_name != ''",
      wellPanel(
        h5("Files to be saved:"),
        tags$ul(
          tags$li(strong("data/Annotation.rds"), " - Main annotation with all results"),
          tags$li(strong("data/Annotation_DT.rds"), " - Annotation formatted for DT display"),
          tags$li(strong("data/Annotation.csv"), " - Annotation in CSV format"),
          tags$li(strong("data/amino_acid_metadata.rds"), " - Metadata about this analysis"),
          tags$li(strong("data_preparation/Amino_Acid_Analysis/smORFs.fa"), " - FASTA file"),
          conditionalPanel(
            condition = "output.deeptmhmm_processed",
            tags$li(strong("data_preparation/Amino_Acid_Analysis/predicted_topologies.3line"), 
                   " - DeepTMHMM results")
          ),
          conditionalPanel(
            condition = "output.targetp_processed",
            tags$li(strong("data_preparation/Amino_Acid_Analysis/output_protein_type.txt"), 
                   " - TargetP results")
          ),
          conditionalPanel(
            condition = "output.deeploc_processed",
            tags$li(strong("data_preparation/Amino_Acid_Analysis/Deeploc_output.csv"), 
                   " - Deeploc results")
          ),
          conditionalPanel(
            condition = "output.interpro_processed",
            tags$li(strong("data_preparation/Amino_Acid_Analysis/interproscan_output.tsv"), 
                   " - Interproscan results")
          )
        ),
        br(),
        actionButton("save_project", "Save All Files to Project Folder", 
                     class = "btn-success btn-lg", style = "width: 100%;"),
        br(), br(),
        conditionalPanel(
          condition = "output.files_saved",
          div(class = "alert alert-success",
              tags$strong("✓ Files saved successfully!"),
              br(),
              textOutput("saved_path_message", inline = TRUE),
              br(),
              "Your project folder is now ready for the next analysis steps or for use with the main Shiny app.")
        )
      )
    ),
    
    hr(),
    h4("Preview & Download"),
    p("You can also download individual files here:"),
    DTOutput("annotation_preview"),
    br(),
    downloadButton("download_annotation_csv", "Download Annotation CSV", class = "btn-primary"),
    downloadButton("download_annotation_rds", "Download Annotation RDS", class = "btn-primary")
  )
}

shinyApp(ui, server)

