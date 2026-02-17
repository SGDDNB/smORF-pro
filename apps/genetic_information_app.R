## Shiny wizard for Genetic Information Analysis
##
## This app guides the user step-by-step:
## 0) Upload Annotation CSV or GTF
## 1) GTF File Creation (if CSV uploaded)
## 2) Conservation Analysis Setup
## 3) Clustal Omega Alignment (Optional)
## 4) GWAS Overlap Analysis
## 5) Integration & Results

library(shiny)
library(DT)

# Source files relative to apps/ directory
source("../config/genetic_information_config.R", local = TRUE)
source("../scripts/process_genetic_information.R", local = TRUE)
source("../scripts/manage_assets.R", local = TRUE)

steps <- c(
  "Welcome & Input",
  "GTF Creation",
  "MAF Indexing",
  "Conservation Analysis",
  "Clustal Omega (Optional)",
  "GWAS Overlap",
  "Integration & Results"
)

ui <- fluidPage(
  titlePanel("smORF-pro: Genetic Information Analysis Wizard"),
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
      # Error display area
      conditionalPanel(
        condition = "output.error_message != ''",
        div(id = "error_display", class = "alert alert-danger",
            tags$strong("Error: "),
            textOutput("error_message", inline = TRUE),
            style = "margin-bottom: 20px;")
      ),
      uiOutput("step_ui")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    step = 0L,
    step_complete = rep(FALSE, length(steps)),
    annotation = NULL,
    gtf_path = NULL,
    exon_file_path = NULL,
    intermediate_gtf_path = NULL,
    conservation_output_dir = NULL,
    gwas_overlaps = NULL,
    files_saved = FALSE,
    error_message = "",
    # Assets file lists (updated on refresh)
    maf_files_list = data.frame(),
    gwas_files_list = data.frame(),
    clustalo_available = FALSE,
    # MAF indexing status
    maf_indexed = FALSE,
    maf_index_output = "",
    maf_index_errors = FALSE
  )
  
  is_maf_indexed <- function(maf_files_list) {
    nrow(maf_files_list) > 0 &&
      "has_index" %in% colnames(maf_files_list) &&
      all(maf_files_list$has_index, na.rm = TRUE)
  }

  update_maf_index_status <- function(maf_files_list) {
    indexed <- is_maf_indexed(maf_files_list)
    rv$maf_indexed <- indexed
    if (indexed) {
      rv$step_complete[3] <- TRUE
    }
  }
  
  # Initialize: Check for Clustal Omega and scan assets folders
  observe({
    rv$clustalo_available <- check_clustalo_available()
    rv$maf_files_list <- check_maf_files(from_apps = TRUE)
    rv$gwas_files_list <- check_gwas_catalog(from_apps = TRUE)
    update_maf_index_status(rv$maf_files_list)
  })
  
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
      "2" = step2_ui(rv, input),
      "3" = step3_ui(rv, input),
      "4" = step4_ui(rv),
      "5" = step5_ui(rv, input),
      "6" = step6_ui(rv),
      step0_ui()
    )
  })
  
  # Step 0: Upload annotation or GTF
  observeEvent(input$upload_annotation, {
    req(input$annotation_file)
    path <- input$annotation_file$datapath
    rv$error_message <- ""
    
    tryCatch({
      annotation <- utils::read.csv(path, stringsAsFactors = FALSE)
      
      # Validate required columns
      required_cols <- c("Chr", "S_exon1", "E_exon1", "strand", "gene_id", "ORF_id", 
                         "gene_name", "iORF_type", "gene_biotype", "Peptide.seq")
      missing_cols <- setdiff(required_cols, colnames(annotation))
      
      if (length(missing_cols) > 0) {
        error_msg <- paste0(
          "Missing required columns in annotation file: ", 
          paste(missing_cols, collapse = ", "), 
          ". Required columns are: ", 
          paste(required_cols, collapse = ", ")
        )
        rv$error_message <- error_msg
        showNotification(error_msg, type = "error", duration = 10)
        return()
      }
      
      rv$annotation <- annotation
      rv$step_complete[1] <- TRUE
      rv$error_message <- ""
      showNotification("Annotation loaded successfully.", type = "message")
    }, error = function(e) {
      error_msg <- paste("Error loading annotation file:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  observeEvent(input$upload_gtf, {
    req(input$gtf_file)
    path <- input$gtf_file$datapath
    rv$error_message <- ""
    
    tryCatch({
      validate_gtf(path)
      rv$gtf_path <- path
      rv$step_complete[1] <- TRUE
      rv$step_complete[2] <- TRUE  # Skip GTF creation if GTF already provided
      rv$error_message <- ""
      showNotification("GTF file loaded successfully.", type = "message")
    }, error = function(e) {
      error_msg <- paste("Error loading GTF file:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 1: Create GTF from annotation
  observeEvent(input$create_gtf, {
    req(rv$annotation)
    
    # Clear previous error
    rv$error_message <- ""
    
    # Validate columns first
    required_cols <- c("Chr", "S_exon1", "E_exon1", "strand", "gene_id", "ORF_id", 
                       "gene_name", "iORF_type", "gene_biotype", "Peptide.seq")
    missing_cols <- setdiff(required_cols, colnames(rv$annotation))
    
    if (length(missing_cols) > 0) {
      error_msg <- paste0(
        "Missing required columns in annotation file: ", 
        paste(missing_cols, collapse = ", "), 
        ". Required columns are: ", 
        paste(required_cols, collapse = ", ")
      )
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
      return()
    }
    
    tryCatch({
      withProgress(message = "Creating GTF file...", value = 0, {
        setProgress(0.2, detail = "Creating exon file...")
        
        # Create exon file
        exon_rel <- genetic_information_config$paths$exon_file_output
        exon_path_from_apps <- file.path("..", exon_rel)
        exon_abs <- normalizePath(exon_path_from_apps, winslash = "/", mustWork = FALSE)
        dir.create(dirname(exon_abs), recursive = TRUE, showWarnings = FALSE)
        
        create_exon_file_from_annotation(rv$annotation, exon_abs)
        rv$exon_file_path <- exon_abs
        
        setProgress(0.4, detail = "Running Perl scripts...")
        
        # Run Perl scripts
        intermediate_rel <- genetic_information_config$paths$intermediate_gtf
        intermediate_path_from_apps <- file.path("..", intermediate_rel)
        intermediate_abs <- normalizePath(intermediate_path_from_apps, winslash = "/", mustWork = FALSE)
        
        final_rel <- genetic_information_config$paths$gtf_output
        final_path_from_apps <- file.path("..", final_rel)
        final_abs <- normalizePath(final_path_from_apps, winslash = "/", mustWork = FALSE)
        
        run_perl_scripts(exon_abs, intermediate_abs, final_abs, from_apps = TRUE)
        
        rv$gtf_path <- final_abs
        rv$intermediate_gtf_path <- intermediate_abs
        rv$step_complete[2] <- TRUE
        
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$error_message <- ""  # Clear error on success
      showNotification("GTF file created successfully.", type = "message")
    }, error = function(e) {
      error_msg <- paste("Error creating GTF:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 2: MAF Indexing - refresh MAF files list
  observeEvent(input$refresh_maf_files, {
    rv$maf_files_list <- check_maf_files(from_apps = TRUE)
    update_maf_index_status(rv$maf_files_list)
    showNotification("MAF files list refreshed.", type = "message")
  })
  
  # Index MAF files
  observeEvent(input$index_maf_files, {
    req(rv$gtf_path)
    
    # Get MAF directory (same logic as conservation analysis)
    maf_files <- NULL
    selected_rows <- input$maf_files_table_rows_selected
    if (!is.null(selected_rows) && length(selected_rows) > 0) {
      maf_files <- rv$maf_files_list$path[selected_rows]
    } else {
      assets_base <- get_assets_base(from_apps = TRUE)
      default_maf_dir <- file.path(assets_base, "maf_files")
      if (dir.exists(default_maf_dir)) {
        maf_files <- default_maf_dir
      } else {
        showNotification("No MAF files found. Please download files first.", type = "warning")
        return()
      }
    }
    
    rv$error_message <- ""
    
    tryCatch({
      withProgress(message = "Indexing MAF files...", value = 0, {
        setProgress(0.1, detail = "Starting indexing (this may take hours)...")
        
        result <- index_maf_files(maf_files, from_apps = TRUE)
        
        rv$maf_index_output <- result$output
        rv$maf_index_errors <- !result$success
        
        setProgress(0.95, detail = "Refreshing file list...")
        
        # Refresh the MAF files list to show renamed files
        rv$maf_files_list <- check_maf_files(from_apps = TRUE)
        update_maf_index_status(rv$maf_files_list)
        
        setProgress(1.0, detail = "Complete!")
      })
      
      if (result$success) {
        rv$step_complete[3] <- TRUE
        showNotification("MAF files indexed successfully!", type = "message", duration = 10)
      } else {
        showNotification("Some MAF files had errors. Check the output for details.", 
                         type = "warning", duration = 15)
      }
    }, error = function(e) {
      error_msg <- paste("Error indexing MAF files:", e$message)
      rv$error_message <- error_msg
      rv$maf_indexed <- FALSE
      rv$maf_index_errors <- TRUE
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 3: Conservation analysis
  
  # Download MAF files
  observeEvent(input$download_maf_files, {
    rv$error_message <- ""
    
    tryCatch({
      withProgress(message = "Downloading MAF files...", value = 0, {
        # Determine which chromosomes to download
        chromosomes <- NULL
        if (!is.null(rv$gtf_path) && file.exists(rv$gtf_path)) {
          setProgress(0.05, detail = "Determining chromosomes from GTF...")
          chromosomes <- get_chromosomes_from_gtf(rv$gtf_path)
        }
        
        # If no GTF or couldn't determine, download all standard chromosomes
        if (is.null(chromosomes) || length(chromosomes) == 0) {
          chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
          setProgress(0.05, detail = "Preparing to download all standard chromosomes...")
        } else {
          setProgress(0.05, detail = paste("Preparing to download", length(chromosomes), "chromosomes..."))
        }
        
        # Create progress callback function
        progress_callback <- function(file_num, total_files, filename) {
          progress_value <- 0.1 + (0.8 * (file_num / total_files))  # 10% to 90% for downloads
          detail_text <- paste0("Downloading ", filename, " (", file_num, " of ", total_files, ")")
          setProgress(progress_value, detail = detail_text)
        }
        
        setProgress(0.1, detail = "Starting downloads...")
        
        # Download files with progress callback
        result <- download_maf_files(chromosomes = chromosomes, from_apps = TRUE, 
                                    progress_callback = progress_callback)
        
        setProgress(0.95, detail = "Refreshing file list...")
        
        # Refresh the MAF files list
        rv$maf_files_list <- check_maf_files(from_apps = TRUE)
        update_maf_index_status(rv$maf_files_list)
        
        setProgress(1.0, detail = "Complete!")
      })
      
      # Show results
      if (length(result$downloaded) > 0) {
        showNotification(
          paste("Downloaded", result$total_downloaded, "MAF file(s).", 
                if (length(result$skipped) > 0) paste(length(result$skipped), "already existed.") else "",
                if (length(result$errors) > 0) paste(length(result$errors), "errors.") else ""),
          type = "message",
          duration = 10
        )
      } else if (result$total_skipped > 0) {
        showNotification(
          paste("All requested MAF files already exist (", result$total_skipped, " files)."),
          type = "message",
          duration = 5
        )
      } else {
        showNotification(
          "No files were downloaded.",
          type = "message",
          duration = 5
        )
      }
      
      if (length(result$errors) > 0) {
        error_msg <- paste("Some downloads failed:", paste(result$errors, collapse = "; "))
        rv$error_message <- error_msg
        showNotification(error_msg, type = "warning", duration = 10)
      }
      
    }, error = function(e) {
      error_msg <- paste("Error downloading MAF files:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  observeEvent(input$run_conservation, {
    req(rv$gtf_path)
    
    # Check if MAF files are indexed
    if (!rv$maf_indexed) {
      showNotification("Please index MAF files first (Step 2).", type = "warning")
      return()
    }
    
    # Get selected MAF files or use default directory
    maf_files <- NULL
    selected_rows <- input$maf_files_table_rows_selected
    if (!is.null(selected_rows) && length(selected_rows) > 0) {
      # Get paths of selected files
      maf_files <- rv$maf_files_list$path[selected_rows]
    } else if (!is.null(input$maf_directory) && nzchar(input$maf_directory)) {
      maf_files <- input$maf_directory
    } else {
      # Use default MAF directory from assets
      assets_base <- get_assets_base(from_apps = TRUE)
      default_maf_dir <- file.path(assets_base, "maf_files")
      if (dir.exists(default_maf_dir)) {
        maf_files <- default_maf_dir
      } else {
        showNotification("No MAF files selected and default directory (assets/maf_files/) not found. Please select files or download them first.", type = "warning")
        return()
      }
    }
    
    if (is.null(maf_files) || length(maf_files) == 0) {
      showNotification("Please select MAF files from the table or specify a directory.", type = "warning")
      return()
    }
    
    rv$error_message <- ""
    
    tryCatch({
      withProgress(message = "Running conservation analysis...", value = 0, {
        setProgress(0.1, detail = "Setting up output directories...")
        
        # Create output directory
        output_rel <- "data_preparation/Genetic_information/Conservation"
        output_path_from_apps <- file.path("..", output_rel)
        output_abs <- normalizePath(output_path_from_apps, winslash = "/", mustWork = FALSE)
        dir.create(output_abs, recursive = TRUE, showWarnings = FALSE)
        
        rv$conservation_output_dir <- output_abs
        
        setProgress(0.3, detail = "Running Python script (this may take hours for first run - MAF indexing)...")
        
        # Run conservation analysis
        run_conservation_analysis(rv$gtf_path, maf_files, output_abs, from_apps = TRUE)
        
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$error_message <- ""
      showNotification("Conservation analysis completed successfully!", type = "message", duration = 10)
      rv$step_complete[4] <- TRUE
    }, error = function(e) {
      error_msg <- paste("Error running conservation analysis:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 3: Clustal Omega - refresh file list
  observeEvent(input$refresh_clustalo_files, {
    # Could scan for existing MSA files
    showNotification("File list refreshed.", type = "message")
  })
  
  observeEvent(input$run_clustalo, {
    req(rv$conservation_output_dir)
    
    if (!rv$clustalo_available) {
      showNotification("Clustal Omega not found. Please activate the conda environment: conda activate smorfpro", 
                       type = "error")
      return()
    }
    
    rv$error_message <- ""
    
    tryCatch({
      withProgress(message = "Running Clustal Omega...", value = 0, {
        setProgress(0.1, detail = "Preparing output folders...")
        
        msa_dir <- file.path(rv$conservation_output_dir, "fasta", "protein", "msa")
        dist_dir <- file.path(rv$conservation_output_dir, "fasta", "protein", "dist_mat")
        dir.create(msa_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)
        
        script_rel <- "data_preparation/Genetic_information/run_clustalo.sh"
        script_path_from_apps <- file.path("..", script_rel)
        script_abs <- normalizePath(script_path_from_apps, winslash = "/", mustWork = FALSE)
        
        if (!file.exists(script_abs)) {
          stop("Clustal Omega script not found: ", script_abs)
        }
        
        setProgress(0.3, detail = "Running Clustal Omega (this may take a while)...")
        
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(rv$conservation_output_dir)
        
        result <- system(paste("bash", shQuote(script_abs)), intern = TRUE)
        
        exit_code <- attr(result, "status")
        if (is.null(exit_code)) exit_code <- 0
        
        if (exit_code != 0) {
          stop("Clustal Omega failed with exit code: ", exit_code)
        }
        
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$error_message <- ""
      rv$step_complete[5] <- TRUE
    }, error = function(e) {
      error_msg <- paste("Error running Clustal Omega:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 4: GWAS overlap - refresh GWAS files list
  observeEvent(input$refresh_gwas_files, {
    rv$gwas_files_list <- check_gwas_catalog(from_apps = TRUE)
    showNotification("GWAS files list refreshed.", type = "message")
  })
  
  observeEvent(input$find_gwas_overlaps, {
    req(rv$gtf_path)
    
    gwas_file <- NULL
    selected_rows <- input$gwas_files_table_rows_selected
    if (!is.null(selected_rows) && length(selected_rows) > 0) {
      # Get path of selected file
      gwas_file <- rv$gwas_files_list$path[selected_rows[1]]
    } else if (!is.null(input$gwas_file_upload)) {
      gwas_file <- input$gwas_file_upload$datapath
    }
    
    if (is.null(gwas_file) || !file.exists(gwas_file)) {
      showNotification("Please select a file from the table or upload a GWAS catalog file.", type = "warning")
      return()
    }
    
    rv$error_message <- ""
    
    tryCatch({
      withProgress(message = "Finding GWAS overlaps...", value = 0, {
        setProgress(0.3, detail = "Reading GTF and GWAS files...")
        
        output_rel <- "data_preparation/Genetic_information/GWAS/GWAS_smORF_overlap.csv"
        output_path_from_apps <- file.path("..", output_rel)
        output_abs <- normalizePath(output_path_from_apps, winslash = "/", mustWork = FALSE)
        
        setProgress(0.6, detail = "Finding overlaps...")
        
        overlaps <- find_gwas_overlaps(rv$gtf_path, gwas_file, output_abs)
        rv$gwas_overlaps <- overlaps
        
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$error_message <- ""
      showNotification(paste("Found", nrow(overlaps), "overlaps."), type = "message")
      rv$step_complete[6] <- TRUE
    }, error = function(e) {
      error_msg <- paste("Error finding GWAS overlaps:", e$message)
      rv$error_message <- error_msg
      showNotification(error_msg, type = "error", duration = 10)
    })
  })
  
  # Step 5: Save to project folder
  observeEvent(input$save_project, {
    req(input$project_name)
    
    if (input$project_name == "" || nchar(trimws(input$project_name)) == 0) {
      showNotification("Please enter a project name.", type = "warning")
      return()
    }
    
    project_name <- gsub("[^A-Za-z0-9_-]", "_", trimws(input$project_name))
    if (project_name == "") {
      showNotification("Project name contains only invalid characters.", type = "error")
      return()
    }
    
    projects_base <- file.path("..", "projects")
    project_path <- file.path(projects_base, project_name)
    project_path_abs <- normalizePath(project_path, mustWork = FALSE)
    
    tryCatch({
      withProgress(message = "Saving files to project folder...", value = 0, {
        setProgress(0.1, detail = "Creating folder structure...")
        dir.create(file.path(project_path_abs, "data"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(project_path_abs, "data_preparation", "Genetic_information"), 
                   recursive = TRUE, showWarnings = FALSE)
        
        # Save GTF file
        if (!is.null(rv$gtf_path) && file.exists(rv$gtf_path)) {
          setProgress(0.2, detail = "Copying GTF file...")
          file.copy(rv$gtf_path,
                    file.path(project_path_abs, "data_preparation", "Genetic_information", "smORFs.gtf"),
                    overwrite = TRUE)
        }
        
        # Copy conservation results if available
        if (!is.null(rv$conservation_output_dir) && dir.exists(rv$conservation_output_dir)) {
          setProgress(0.4, detail = "Copying conservation results...")
          conservation_dest <- file.path(project_path_abs, "data_preparation", "Genetic_information", "Conservation")
          if (dir.exists(conservation_dest)) {
            unlink(conservation_dest, recursive = TRUE)
          }
          file.copy(rv$conservation_output_dir, 
                    file.path(project_path_abs, "data_preparation", "Genetic_information"),
                    recursive = TRUE, overwrite = TRUE)
        }
        
        # Save GWAS overlaps if available
        if (!is.null(rv$gwas_overlaps)) {
          setProgress(0.6, detail = "Saving GWAS overlaps...")
          gwas_dest <- file.path(project_path_abs, "data_preparation", "Genetic_information", "GWAS")
          dir.create(gwas_dest, recursive = TRUE, showWarnings = FALSE)
          write.csv(rv$gwas_overlaps,
                    file.path(gwas_dest, "GWAS_smORF_overlap.csv"),
                    row.names = FALSE)
        }
        
        # Create metadata
        setProgress(0.9, detail = "Creating metadata...")
        metadata <- list(
          project_name = project_name,
          project_path = project_path_abs,
          created_date = as.character(Sys.Date()),
          created_time = as.character(Sys.time()),
          analyses_completed = list(
            gtf_created = !is.null(rv$gtf_path) && file.exists(rv$gtf_path),
            conservation_run = !is.null(rv$conservation_output_dir) && dir.exists(rv$conservation_output_dir),
            gwas_overlaps = !is.null(rv$gwas_overlaps)
          ),
          gtf_path = if (!is.null(rv$gtf_path)) rv$gtf_path else NA_character_,
          conservation_dir = if (!is.null(rv$conservation_output_dir)) rv$conservation_output_dir else NA_character_,
          num_gwas_overlaps = if (!is.null(rv$gwas_overlaps)) nrow(rv$gwas_overlaps) else 0
        )
        saveRDS(metadata, file.path(project_path_abs, "data", "genetic_information_metadata.rds"))
        
        setProgress(1.0, detail = "Complete!")
      })
      
      rv$files_saved <- TRUE
      showNotification(
        paste("All files saved successfully to:", project_path_abs),
        type = "message",
        duration = 10
      )
    }, error = function(e) {
      showNotification(paste("Error saving files:", e$message), type = "error", duration = 10)
      rv$files_saved <- FALSE
    })
  })
  
  # Outputs
  output$gtf_created <- reactive({
    !is.null(rv$gtf_path) && file.exists(rv$gtf_path)
  })
  outputOptions(output, "gtf_created", suspendWhenHidden = FALSE)
  
  output$conservation_completed <- reactive({
    !is.null(rv$conservation_output_dir) && dir.exists(rv$conservation_output_dir)
  })
  outputOptions(output, "conservation_completed", suspendWhenHidden = FALSE)
  
  get_conservation_summary <- function(output_dir) {
    expected_files <- c(
      "start_codon_conservation.txt",
      "percentage_ids.txt",
      "stop_codon_conservation.txt",
      "inframe.txt",
      "protein_percentage_ids.txt"
    )
    paths <- file.path(output_dir, expected_files)
    exists <- file.exists(paths)
    size_mb <- ifelse(exists, round(file.info(paths)$size / 1024 / 1024, 2), NA_real_)
    data.frame(
      file = expected_files,
      exists = exists,
      size_mb = size_mb,
      stringsAsFactors = FALSE
    )
  }
  
  output$conservation_summary <- renderTable({
    req(rv$conservation_output_dir)
    get_conservation_summary(rv$conservation_output_dir)
  })
  
  
  output$gwas_completed <- reactive({
    !is.null(rv$gwas_overlaps)
  })
  outputOptions(output, "gwas_completed", suspendWhenHidden = FALSE)
  
  output$clustalo_available <- reactive({
    rv$clustalo_available
  })
  outputOptions(output, "clustalo_available", suspendWhenHidden = FALSE)
  
  output$files_saved <- reactive({
    rv$files_saved
  })
  outputOptions(output, "files_saved", suspendWhenHidden = FALSE)
  
  output$annotation_preview <- DT::renderDataTable({
    req(rv$annotation)
    DT::datatable(rv$annotation, options = list(pageLength = 10))
  })
  
  output$gwas_overlaps_preview <- DT::renderDataTable({
    req(rv$gwas_overlaps)
    DT::datatable(rv$gwas_overlaps, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$maf_files_table <- DT::renderDataTable({
    req(rv$maf_files_list)
    if (nrow(rv$maf_files_list) == 0) {
      return(data.frame(Message = "No MAF files found in assets/maf_files/"))
    }
    
    # Create display table with gzip indicator
    display_df <- rv$maf_files_list[, c("file", "size_mb", "modified")]
    
    # Add gzip indicator to filename if present
    if ("is_gzipped" %in% colnames(rv$maf_files_list)) {
      display_df$file <- paste0(
        display_df$file,
        ifelse(rv$maf_files_list$is_gzipped, " (gzipped)", "")
      )
    }
    
    DT::datatable(display_df, 
                  options = list(pageLength = 10),
                  selection = "multiple",
                  colnames = c("File", "Size (MB)", "Modified"))
  })
  
  output$gwas_files_table <- DT::renderDataTable({
    req(rv$gwas_files_list)
    if (nrow(rv$gwas_files_list) == 0) {
      return(data.frame(Message = "No GWAS catalog files found in assets/gwas_catalog/"))
    }
    DT::datatable(rv$gwas_files_list[, c("file", "size_mb", "modified")], 
                  options = list(pageLength = 5),
                  selection = "single")
  })
  
  output$num_overlaps <- renderText({
    if (!is.null(rv$gwas_overlaps)) {
      as.character(nrow(rv$gwas_overlaps))
    } else {
      "0"
    }
  })
  
  output$annotation_loaded <- reactive({
    !is.null(rv$annotation)
  })
  outputOptions(output, "annotation_loaded", suspendWhenHidden = FALSE)
  
  output$error_message <- renderText({
    rv$error_message
  })
  outputOptions(output, "error_message", suspendWhenHidden = FALSE)
  
  # MAF indexing outputs
  output$maf_indexed <- reactive({
    rv$maf_indexed
  })
  outputOptions(output, "maf_indexed", suspendWhenHidden = FALSE)
  
  output$maf_index_errors <- reactive({
    rv$maf_index_errors
  })
  outputOptions(output, "maf_index_errors", suspendWhenHidden = FALSE)
  
  output$maf_index_output <- renderText({
    rv$maf_index_output
  })
  outputOptions(output, "maf_index_output", suspendWhenHidden = FALSE)
}

# UI Helper Functions

fileSelectionUI <- function(id, file_list, file_type, multi_select = FALSE) {
  ns <- NS(id)
  
  tagList(
    h5(paste("Available", file_type, "files:")),
    if (nrow(file_list) == 0) {
      div(class = "alert alert-info",
          paste("No", file_type, "files found. Please download files and save them to assets/", 
                tolower(gsub(" ", "_", file_type)), "/", sep = ""))
    } else {
      tagList(
        DT::dataTableOutput(ns("files_table")),
        br(),
        actionButton(ns("refresh"), "Refresh file list", class = "btn-sm")
      )
    }
  )
}

step0_ui <- function() {
  tagList(
    h3("Step 0: Welcome & Input"),
    p("This wizard will guide you through the genetic information analysis workflow used in smORF-pro."),
    p("You will:"),
    tags$ul(
      tags$li("Upload your smORF annotation CSV (to create GTF) OR upload an existing GTF file"),
      tags$li("Run conservation analysis using UCSC MAF files"),
      tags$li("Optionally run Clustal Omega for multiple sequence alignments"),
      tags$li("Find overlaps with GWAS catalog mutations"),
      tags$li("Save all results to a project folder")
    ),
    hr(),
    h4("Option 1: Upload Annotation CSV"),
    p("If you have an annotation CSV with exon coordinates, we can create a GTF file."),
    p(strong("Required columns:"), 
      "Chr, S_exon1, E_exon1, strand, gene_id, ORF_id, gene_name, iORF_type, gene_biotype, Peptide.seq"),
    p(strong("Optional columns:"), 
      "S_exon2, E_exon2, S_exon3, E_exon3 (for additional exons)"),
    p("See Example/Example_smORFs.csv for the expected format."),
    fileInput("annotation_file", "Upload smORF annotation CSV", accept = c(".csv")),
    actionButton("upload_annotation", "Load annotation"),
    br(), br(),
    conditionalPanel(
      condition = "output.annotation_loaded",
      DTOutput("annotation_preview")
    ),
    hr(),
    h4("Option 2: Upload Existing GTF File"),
    p("If you already have a GTF file, you can upload it directly:"),
    fileInput("gtf_file", "Upload GTF file", accept = c(".gtf")),
    actionButton("upload_gtf", "Load GTF")
  )
}

step1_ui <- function(rv) {
  tagList(
    h3("Step 1: GTF File Creation"),
    p("We will create a GTF file from your annotation for conservation and GWAS analysis."),
    p("This step will:"),
    tags$ul(
      tags$li("Create an exon file from your annotation"),
      tags$li("Run Perl scripts to convert to GTF format"),
      tags$li("Save the GTF file for use in subsequent steps")
    ),
    br(),
    conditionalPanel(
      condition = "!output.gtf_created",
      actionButton("create_gtf", "Create GTF File", class = "btn-primary btn-lg")
    ),
    conditionalPanel(
      condition = "output.gtf_created",
      div(class = "alert alert-success",
          tags$strong("✓ GTF file created successfully!"),
          br(),
          "GTF file path: ",
          code(if (!is.null(rv$gtf_path)) rv$gtf_path else "N/A")
      )
    )
  )
}

step2_ui <- function(rv, input) {
  tagList(
    h3("Step 2: MAF File Indexing"),
    p("Before running conservation analysis, MAF files must be indexed. This step will:"),
    tags$ul(
      tags$li("Detect if MAF files are gzipped (compressed)"),
      tags$li("Automatically decompress gzipped files if needed"),
      tags$li("Check each MAF file for issues (corruption, encoding problems)"),
      tags$li("Create index files (.mafindex) for fast access"),
      tags$li("Report any problems found"),
      tags$li("Note: This may take hours/days on first run")
    ),
    hr(),
    h4("MAF Files"),
    p("Files in assets/maf_files/ (both .maf and .maf.gz files are supported):"),
    DT::dataTableOutput("maf_files_table"),
    br(),
    div(
      actionButton("refresh_maf_files", "Refresh file list", class = "btn-sm"),
      actionButton("download_maf_files", "Download Missing MAF Files", 
                   class = "btn-primary btn-sm", style = "margin-left: 10px;"),
      br(), br(),
      helpText("Downloads all standard chromosomes (chr1-chr22, chrX, chrY, chrM) or only those needed based on your GTF file."),
      helpText("Note: Downloaded files may be gzipped - they will be automatically decompressed during indexing.")
    ),
    hr(),
    conditionalPanel(
      condition = "output.gtf_created",
      actionButton("index_maf_files", "Index MAF Files", 
                   class = "btn-primary btn-lg"),
      br(), br(),
      conditionalPanel(
        condition = "output.maf_indexed",
        div(class = "alert alert-success",
            tags$strong("✓ MAF files indexed successfully!"),
            br(),
            "You can now proceed to conservation analysis.")
      ),
      conditionalPanel(
        condition = "output.maf_index_errors",
        div(class = "alert alert-danger",
            tags$strong("⚠ Some MAF files had errors during indexing."),
            br(),
            "Check the output below for details. You may need to re-download problematic files.")
      )
    ),
    conditionalPanel(
      condition = "!output.gtf_created",
      div(class = "alert alert-warning",
          "Please complete GTF creation first.")
    ),
    hr(),
    conditionalPanel(
      condition = "output.maf_index_output != ''",
      h4("Indexing Output"),
      verbatimTextOutput("maf_index_output", placeholder = TRUE)
    )
  )
}

step3_ui <- function(rv, input) {
  tagList(
    h3("Step 3: Conservation Analysis"),
    p("Conservation analysis aligns smORF CDS with 100 vertebrates to find conserved regions."),
    p(strong("Note:"), "MAF files must be indexed first (Step 2)."),
    hr(),
    h4("MAF Files"),
    p("Files in assets/maf_files/:"),
    DT::dataTableOutput("maf_files_table"),
    br(),
    div(
      actionButton("refresh_maf_files", "Refresh file list", class = "btn-sm"),
      br(), br()
    ),
    hr(),
    h4("Alternative: Specify MAF Directory"),
    textInput("maf_directory", "Or specify custom MAF directory path:", 
              placeholder = "/path/to/maf/files",
              width = "100%"),
    hr(),
    conditionalPanel(
      condition = "output.maf_indexed",
      actionButton("run_conservation", "Run Conservation Analysis", 
                   class = "btn-primary btn-lg")
    ),
    conditionalPanel(
      condition = "!output.maf_indexed",
      div(class = "alert alert-warning",
          "Please complete MAF indexing first (Step 2).")
    ),
    conditionalPanel(
      condition = "output.conservation_completed",
      div(class = "alert alert-success",
          tags$strong("✓ Conservation analysis completed!"),
          br(),
          "Results saved to: ",
          code(if (!is.null(rv$conservation_output_dir)) rv$conservation_output_dir else "N/A")
      )
    )
    ,
    conditionalPanel(
      condition = "output.conservation_completed",
      h4("Conservation output summary"),
      p("Expected output files and their status:"),
      tableOutput("conservation_summary")
    )
  )
}

step4_ui <- function(rv) {
  tagList(
    h3("Step 4: Clustal Omega Alignment (Optional)"),
    p("Clustal Omega creates multiple sequence alignments and distance matrices from protein alignments."),
    div(class = "alert alert-warning",
        tags$strong("Warning: "),
        "Running this step will overwrite existing MSA and distance matrix outputs."),
    conditionalPanel(
      condition = "output.clustalo_available",
      div(class = "alert alert-success",
          tags$strong("✓ Clustal Omega is available!"),
          " (installed via conda)")
    ),
    conditionalPanel(
      condition = "!output.clustalo_available",
      div(class = "alert alert-warning",
          tags$strong("Clustal Omega not found."),
          br(),
          "Make sure you've activated the conda environment: ",
          code("conda activate smorfpro"),
          br(),
          "Clustal Omega is already included in environment.yml and will be installed automatically."
      )
    ),
    hr(),
    conditionalPanel(
      condition = "output.conservation_completed",
      p("Conservation results are available. You can run Clustal Omega on the protein alignments."),
      actionButton("run_clustalo", "Run Clustal Omega", class = "btn-primary"),
      br(), br(),
      actionButton("refresh_clustalo_files", "Refresh file list", class = "btn-sm")
    ),
    conditionalPanel(
      condition = "!output.conservation_completed",
      div(class = "alert alert-info",
          "Please complete conservation analysis first.")
    )
  )
}

step5_ui <- function(rv, input) {
  tagList(
    h3("Step 5: GWAS Overlap Analysis"),
    p("Find overlaps between smORF CDS and known GWAS mutations."),
    
    h4("Instructions"),
    tags$ol(
      tags$li("Download GWAS catalog from:"),
      tags$li(tags$a(href = genetic_information_config$urls$gwas_catalog, 
                     target = "_blank", 
                     genetic_information_config$urls$gwas_catalog)),
      tags$li("Download 'All Associations v1.0' (TSV format, ~300+ MB)"),
      tags$li("Save to assets/gwas_catalog/ folder"),
      tags$li("Or select from existing files below or upload a new file")
    ),
    hr(),
    h4("GWAS Catalog Files"),
    p("Files in assets/gwas_catalog/:"),
    DT::dataTableOutput("gwas_files_table"),
    br(),
    actionButton("refresh_gwas_files", "Refresh file list", class = "btn-sm"),
    hr(),
    h4("Or Upload New File"),
    fileInput("gwas_file_upload", "Upload GWAS catalog TSV file", accept = c(".tsv", ".txt")),
    hr(),
    conditionalPanel(
      condition = "output.gtf_created",
      actionButton("find_gwas_overlaps", "Find Overlaps", class = "btn-primary btn-lg")
    ),
    conditionalPanel(
      condition = "!output.gtf_created",
      div(class = "alert alert-warning",
          "Please complete GTF creation first.")
    ),
    conditionalPanel(
      condition = "output.gwas_completed",
      div(class = "alert alert-success",
          tags$strong("✓ GWAS overlap analysis completed!"),
          br(),
          "Found ", 
          textOutput("num_overlaps", inline = TRUE),
          " overlaps."
      ),
      hr(),
      h4("Overlap Results Preview"),
      DTOutput("gwas_overlaps_preview")
    )
  )
}

step6_ui <- function(rv) {
  tagList(
    h3("Step 6: Integration & Results"),
    p("All available results have been processed. Save everything to a project folder."),
    
    h4("Save to Project Folder"),
    p("Enter a project name. All files will be saved to projects/[project_name]/ with standardized names:"),
    
    textInput("project_name", "Project name:", 
              placeholder = "e.g., my_genetic_analysis",
              width = "100%"),
    helpText("A folder will be created at: projects/[project_name]/"),
    br(),
    
    conditionalPanel(
      condition = "input.project_name != ''",
      wellPanel(
        h5("Files to be saved:"),
        tags$ul(
          conditionalPanel(
            condition = "output.gtf_created",
            tags$li(strong("data_preparation/Genetic_information/smORFs.gtf"), 
                   " - GTF file")
          ),
          conditionalPanel(
            condition = "output.conservation_completed",
            tags$li(strong("data_preparation/Genetic_information/Conservation/"), 
                   " - Conservation analysis results")
          ),
          conditionalPanel(
            condition = "output.gwas_completed",
            tags$li(strong("data_preparation/Genetic_information/GWAS/GWAS_smORF_overlap.csv"), 
                   " - GWAS overlap results")
          ),
          tags$li(strong("data/genetic_information_metadata.rds"), 
                 " - Metadata about this analysis")
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
              "Your project folder is now ready for use with the main Shiny app.")
        )
      )
    ),
    
    hr(),
    conditionalPanel(
      condition = "output.gwas_completed",
      h4("GWAS Overlap Results"),
      DTOutput("gwas_overlaps_preview"),
      br()
    )
  )
}

shinyApp(ui, server)
