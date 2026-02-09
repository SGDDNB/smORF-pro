## Helper functions for genetic information analysis processing
## These functions encapsulate the code currently shown in README / README.Rmd
## so they can be reused by scripts and Shiny apps.

suppressed_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Please install it first.", pkg), call. = FALSE)
  }
}

#' Create exon file format from annotation data.frame
#'
#' @param annotation data.frame with exon coordinate columns (S_exon1, E_exon1, etc.)
#' @param output_path path to output exon file
#' @return invisible output_path
create_exon_file_from_annotation <- function(annotation, output_path) {
  suppressed_require("data.table")
  
  required_cols <- c("Chr", "S_exon1", "E_exon1", "strand", "gene_id", "ORF_id", 
                     "gene_name", "iORF_type", "gene_biotype", "Peptide.seq")
  missing <- setdiff(required_cols, colnames(annotation))
  if (length(missing) > 0) {
    stop("Annotation is missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Create 1 line per exon (handle any number of exons)
  # Find all exon columns (S_exon1, E_exon1, S_exon2, E_exon2, etc.)
  exon_start_cols <- grep("^S_exon[0-9]+$", colnames(annotation), value = TRUE)
  exon_end_cols <- grep("^E_exon[0-9]+$", colnames(annotation), value = TRUE)
  
  # Extract exon numbers and find maximum
  exon_numbers <- unique(as.numeric(gsub("^S_exon|^E_exon", "", c(exon_start_cols, exon_end_cols))))
  exon_numbers <- exon_numbers[!is.na(exon_numbers)]
  max_exon <- if (length(exon_numbers) > 0) max(exon_numbers) else 1
  
  # Base columns needed for all exons (select by name, not position)
  base_cols <- c("Chr", "strand", "gene_id", "ORF_id", "gene_name", 
                 "iORF_type", "gene_biotype", "Peptide.seq")
  
  # Create a data frame for each exon
  exon_list <- list()
  
  for (exon_num in 1:max_exon) {
    s_col <- paste0("S_exon", exon_num)
    e_col <- paste0("E_exon", exon_num)
    
    # Check if both start and end columns exist for this exon
    if (s_col %in% colnames(annotation) && e_col %in% colnames(annotation)) {
      # Select base columns plus this exon's start/end columns
      exon_cols <- c(base_cols, s_col, e_col)
      exon_data <- annotation[, exon_cols, drop = FALSE]
      
      # Rename exon columns to S_exon1/E_exon1 for consistency when combining
      colnames(exon_data)[colnames(exon_data) == s_col] <- "S_exon1"
      colnames(exon_data)[colnames(exon_data) == e_col] <- "E_exon1"
      
      exon_list <- c(exon_list, list(exon_data))
    }
  }
  
  # If no exons found, use exon1 (required column)
  if (length(exon_list) == 0) {
    exon_cols <- c(base_cols, "S_exon1", "E_exon1")
    exon_data <- annotation[, exon_cols, drop = FALSE]
    exon_list <- list(exon_data)
  }
  
  # Combine all exons
  exon_data <- data.table::rbindlist(exon_list, fill = TRUE)
  exon_data$iORF_id <- exon_data$ORF_id
  
  # Remove lines with no exon information
  exon_data <- exon_data[!is.na(exon_data$S_exon1), ]
  exon_data <- exon_data[order(exon_data$iORF_id), ]
  
  # Create exon file in format for Perl script
  exon_file <- data.frame(
    chr = exon_data$Chr,
    start = exon_data$S_exon1,
    end = exon_data$E_exon1,
    strand = exon_data$strand,
    gene = exon_data$gene_id,
    ORF = exon_data$iORF_id,
    exon_id = paste0(exon_data$Chr, exon_data$S_exon1, exon_data$E_exon1),
    gene_name = exon_data$gene_name,
    orf_type = exon_data$iORF_type,
    gene_type = exon_data$gene_biotype,
    pept = exon_data$Peptide.seq,
    Source = "Ho et al. 2020"
  )
  
  # Ensure output directory exists
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  write.table(exon_file, output_path, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t")
  
  invisible(output_path)
}

#' Run Perl scripts to convert exon file to GTF
#'
#' @param exon_file_path path to exon file
#' @param intermediate_gtf_path path for intermediate GTF (from exonToGTF.pl)
#' @param final_gtf_path path for final GTF (from collapseGTF.pl)
#' @param perl_script_dir directory containing Perl scripts (default: data_preparation/Genetic_information)
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return invisible final_gtf_path
run_perl_scripts <- function(exon_file_path, intermediate_gtf_path, final_gtf_path,
                            perl_script_dir = NULL, from_apps = FALSE) {
  if (is.null(perl_script_dir)) {
    if (from_apps) {
      perl_script_dir <- normalizePath("../data_preparation/Genetic_information", 
                                     winslash = "/", mustWork = FALSE)
    } else {
      perl_script_dir <- normalizePath("data_preparation/Genetic_information", 
                                     winslash = "/", mustWork = FALSE)
    }
  }
  
  exon_to_gtf_script <- file.path(perl_script_dir, "exonToGTF.pl")
  collapse_gtf_script <- file.path(perl_script_dir, "collapseGTF.pl")
  
  if (!file.exists(exon_to_gtf_script)) {
    stop("Perl script not found: ", exon_to_gtf_script, call. = FALSE)
  }
  if (!file.exists(collapse_gtf_script)) {
    stop("Perl script not found: ", collapse_gtf_script, call. = FALSE)
  }
  
  # Ensure output directories exist
  dir.create(dirname(intermediate_gtf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(final_gtf_path), recursive = TRUE, showWarnings = FALSE)
  
  # Run first Perl script: exonToGTF.pl
  source <- "Ho_et_al_2020"
  cmd1 <- sprintf('perl "%s" "%s" "%s" "%s"',
                  exon_to_gtf_script, intermediate_gtf_path, exon_file_path, source)
  result1 <- system(cmd1, intern = TRUE)
  
  if (attr(result1, "status") != 0 && !is.null(attr(result1, "status"))) {
    stop("Error running exonToGTF.pl. Command: ", cmd1, call. = FALSE)
  }
  
  # Run second Perl script: collapseGTF.pl
  cmd2 <- sprintf('perl "%s" "%s" "%s"',
                  collapse_gtf_script, final_gtf_path, intermediate_gtf_path)
  result2 <- system(cmd2, intern = TRUE)
  
  if (attr(result2, "status") != 0 && !is.null(attr(result2, "status"))) {
    stop("Error running collapseGTF.pl. Command: ", cmd2, call. = FALSE)
  }
  
  if (!file.exists(final_gtf_path)) {
    stop("GTF file was not created: ", final_gtf_path, call. = FALSE)
  }
  
  invisible(final_gtf_path)
}

#' Basic validation of GTF file format
#'
#' @param gtf_path path to GTF file
#' @return TRUE if passes basic checks, otherwise stops with error
validate_gtf <- function(gtf_path) {
  if (is.null(gtf_path) || !nzchar(gtf_path)) {
    stop("No GTF file specified.", call. = FALSE)
  }
  if (!file.exists(gtf_path)) {
    stop("GTF file does not exist: ", gtf_path, call. = FALSE)
  }
  
  # Check file extension
  ext <- tolower(tools::file_ext(gtf_path))
  if (ext != "gtf") {
    warning("File does not have .gtf extension: ", gtf_path)
  }
  
  # Check if file is readable and has content
  info <- file.info(gtf_path)
  if (info$size == 0) {
    stop("GTF file is empty: ", gtf_path, call. = FALSE)
  }
  
  # Try to read first few lines to check format
  con <- file(gtf_path, "r")
  first_lines <- readLines(con, n = 5)
  close(con)
  
  if (length(first_lines) == 0) {
    stop("GTF file appears to be empty or unreadable: ", gtf_path, call. = FALSE)
  }
  
  # GTF files should have 9 tab-separated columns
  # Skip comment lines starting with #
  non_comment_lines <- first_lines[!grepl("^#", first_lines)]
  if (length(non_comment_lines) > 0) {
    fields <- strsplit(non_comment_lines[1], "\t")[[1]]
    if (length(fields) < 9) {
      warning("GTF file may not be in correct format. Expected 9 tab-separated columns, found ", 
              length(fields), " in first line.")
    }
  }
  
  return(TRUE)
}

#' Index MAF files with diagnostics
#'
#' @param maf_files_list list of MAF file paths (or directory containing MAF files)
#' @param chromosomes vector of specific chromosomes to index (NULL = all found)
#' @param python_script_path path to Python indexing script
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return list with indexing results
index_maf_files <- function(maf_files_list, chromosomes = NULL,
                            python_script_path = NULL, from_apps = FALSE) {
  if (is.null(python_script_path)) {
    if (from_apps) {
      python_script_path <- normalizePath("../data_preparation/Genetic_information/index_maf_files.py",
                                         winslash = "/", mustWork = FALSE)
    } else {
      python_script_path <- normalizePath("data_preparation/Genetic_information/index_maf_files.py",
                                         winslash = "/", mustWork = FALSE)
    }
  }
  
  if (!file.exists(python_script_path)) {
    stop("Python indexing script not found: ", python_script_path, call. = FALSE)
  }
  
  # Determine MAF directory
  maf_dir <- NULL
  if (length(maf_files_list) == 1 && dir.exists(maf_files_list)) {
    maf_dir <- normalizePath(maf_files_list, winslash = "/", mustWork = TRUE)
  } else if (length(maf_files_list) > 0 && all(file.exists(maf_files_list))) {
    maf_dir <- unique(dirname(normalizePath(maf_files_list, winslash = "/", mustWork = TRUE)))
    if (length(maf_dir) > 1) {
      stop("MAF files are in multiple directories. Please provide a single directory or files from the same directory.", call. = FALSE)
    }
    maf_dir <- maf_dir[1]
  } else {
    stop("Invalid MAF files specification. Provide either a directory path or a list of MAF file paths.", call. = FALSE)
  }
  
  # Build Python command
  python_cmd <- sprintf('python "%s" --maf-dir "%s"',
                        python_script_path, maf_dir)
  
  if (!is.null(chromosomes) && length(chromosomes) > 0) {
    chr_args <- paste(chromosomes, collapse = " ")
    python_cmd <- paste(python_cmd, "--chromosomes", chr_args)
  }
  
  message("Indexing MAF files in: ", maf_dir)
  message("This may take several hours for the first run...")
  
  # Execute Python script and capture output
  result <- system(python_cmd, intern = TRUE)
  
  # Check exit code
  exit_code <- attr(result, "status")
  if (is.null(exit_code)) exit_code <- 0
  
  # Parse output
  output_text <- paste(result, collapse = "\n")
  
  # Return results
  list(
    exit_code = exit_code,
    output = output_text,
    maf_dir = maf_dir,
    success = exit_code == 0
  )
}

#' Run conservation analysis using Python script
#'
#' @param gtf_path path to GTF file
#' @param maf_files_list list of MAF file paths (or directory containing MAF files)
#' @param output_dir directory for output (fasta and protein subdirectories)
#' @param python_script_path path to Python script (default: bio_parse_orf_gtf_Hubner_no_length_restriction_maf.py)
#' @param species_file path to species_names.txt file
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return invisible output_dir
run_conservation_analysis <- function(gtf_path, maf_files_list, output_dir,
                                     python_script_path = NULL, species_file = NULL, from_apps = FALSE) {
  if (is.null(python_script_path)) {
    if (from_apps) {
      python_script_path <- normalizePath("../data_preparation/Genetic_information/bio_parse_orf_gtf_Hubner_no_length_restriction_maf.py",
                                         winslash = "/", mustWork = FALSE)
    } else {
      python_script_path <- normalizePath("data_preparation/Genetic_information/bio_parse_orf_gtf_Hubner_no_length_restriction_maf.py",
                                         winslash = "/", mustWork = FALSE)
    }
  }
  
  if (!file.exists(python_script_path)) {
    stop("Python script not found: ", python_script_path, call. = FALSE)
  }
  
  # Get species file path
  if (is.null(species_file)) {
    # Use manage_assets function if available, otherwise construct path directly
    if (exists("get_species_names_path")) {
      species_file <- get_species_names_path(from_apps = from_apps)
    } else {
      # Construct path directly
      if (from_apps) {
        assets_base <- normalizePath("../assets", winslash = "/", mustWork = FALSE)
        repo_path <- normalizePath("../data_preparation/Genetic_information/species_names.txt", 
                                  winslash = "/", mustWork = FALSE)
      } else {
        assets_base <- normalizePath("assets", winslash = "/", mustWork = FALSE)
        repo_path <- normalizePath("data_preparation/Genetic_information/species_names.txt", 
                                  winslash = "/", mustWork = FALSE)
      }
      default_path <- file.path(assets_base, "species_names", "species_names.txt")
      
      if (file.exists(default_path)) {
        species_file <- default_path
      } else if (file.exists(repo_path)) {
        species_file <- repo_path
      } else {
        species_file <- default_path
      }
    }
  }
  
  if (!file.exists(species_file)) {
    stop("Species names file not found: ", species_file, call. = FALSE)
  }
  
  # Determine MAF directory
  maf_dir <- NULL
  if (length(maf_files_list) == 1 && dir.exists(maf_files_list)) {
    # Single directory path
    maf_dir <- normalizePath(maf_files_list, winslash = "/", mustWork = TRUE)
  } else if (length(maf_files_list) > 0 && all(file.exists(maf_files_list))) {
    # List of MAF files - use their common directory
    maf_dir <- unique(dirname(normalizePath(maf_files_list, winslash = "/", mustWork = TRUE)))
    if (length(maf_dir) > 1) {
      stop("MAF files are in multiple directories. Please provide a single directory or files from the same directory.", call. = FALSE)
    }
    maf_dir <- maf_dir[1]
  } else {
    stop("Invalid MAF files specification. Provide either a directory path or a list of MAF file paths.", call. = FALSE)
  }
  
  # Validate GTF file
  if (!file.exists(gtf_path)) {
    stop("GTF file not found: ", gtf_path, call. = FALSE)
  }
  
  # Normalize paths
  gtf_path <- normalizePath(gtf_path, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  species_file <- normalizePath(species_file, winslash = "/", mustWork = TRUE)
  
  # Create output directories
  fasta_dir <- file.path(output_dir, "fasta")
  protein_dir <- file.path(output_dir, "fasta", "protein")
  dir.create(fasta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(protein_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build Python command
  python_cmd <- sprintf('python "%s" --gtf "%s" --maf-dir "%s" --output-dir "%s" --species "%s"',
                        python_script_path, gtf_path, maf_dir, output_dir, species_file)
  
  # Execute Python script
  message("Running conservation analysis...")
  message("GTF file: ", gtf_path)
  message("MAF directory: ", maf_dir)
  message("Output directory: ", output_dir)
  message("Species file: ", species_file)
  message("This may take several hours for the first run (MAF indexing)...")
  
  result <- system(python_cmd, intern = FALSE)
  
  if (result != 0) {
    stop("Python script failed with exit code: ", result, call. = FALSE)
  }
  
  # Check if output files were created
  expected_files <- c("start_codon_conservation.txt", "percentage_ids.txt", 
                     "stop_codon_conservation.txt", "inframe.txt", 
                     "protein_percentage_ids.txt")
  created_files <- file.path(output_dir, expected_files)
  missing_files <- expected_files[!file.exists(created_files)]
  
  if (length(missing_files) > 0) {
    warning("Some expected output files were not created: ", paste(missing_files, collapse = ", "))
  }
  
  message("Conservation analysis completed. Results saved to: ", output_dir)
  
  invisible(output_dir)
}

#' Clean conservation output files (remove quotes from filenames)
#'
#' @param fasta_dir directory containing fasta files
#' @param protein_dir directory containing protein files
#' @return invisible list of renamed files
clean_conservation_files <- function(fasta_dir, protein_dir) {
  renamed <- list(fasta = character(0), protein = character(0))
  
  # Clean nucleotide fasta files
  if (dir.exists(fasta_dir)) {
    fasta_files <- list.files(fasta_dir, pattern = "\\.fa$", full.names = TRUE)
    for (file in fasta_files) {
      new_name <- gsub('"', "", basename(file))
      if (new_name != basename(file)) {
        new_path <- file.path(dirname(file), new_name)
        file.rename(file, new_path)
        renamed$fasta <- c(renamed$fasta, new_path)
      }
    }
  }
  
  # Clean protein files
  if (dir.exists(protein_dir)) {
    protein_files <- list.files(protein_dir, pattern = "\\.txt$", full.names = TRUE)
    for (file in protein_files) {
      new_name <- gsub('"', "", basename(file))
      if (new_name != basename(file)) {
        new_path <- file.path(dirname(file), new_name)
        file.rename(file, new_path)
        renamed$protein <- c(renamed$protein, new_path)
      }
    }
  }
  
  invisible(renamed)
}

#' Convert species names in protein files from genome names to common names
#'
#' @param protein_dir directory containing protein alignment files
#' @param species_file path to species_names.txt file
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return invisible number of files processed
convert_species_names <- function(protein_dir, species_file = NULL, from_apps = FALSE) {
  suppressed_require("data.table")
  suppressed_require("stringr")
  
  if (is.null(species_file)) {
    # Use function from manage_assets.R
    source("../scripts/manage_assets.R", local = TRUE)
    species_file <- get_species_names_path(from_apps)
  }
  
  if (!file.exists(species_file)) {
    stop("Species names file not found: ", species_file, call. = FALSE)
  }
  
  species <- data.table::fread(species_file, sep = "\t")
  
  if (!dir.exists(protein_dir)) {
    stop("Protein directory does not exist: ", protein_dir, call. = FALSE)
  }
  
  protein_files <- list.files(protein_dir, pattern = "\\.txt$", full.names = TRUE)
  processed <- 0
  
  for (file_path in protein_files) {
    file_i <- data.table::fread(file_path, header = FALSE)
    
    # Remove empty lines
    idx_to_rm <- which(file_i[, 1] == "")
    if (length(idx_to_rm) > 0) {
      idx_to_rm <- c(idx_to_rm, idx_to_rm - 1)
      file_i <- file_i[-idx_to_rm, ]
    }
    
    order_to_move <- c()
    for (sp in 1:nrow(species)) {
      genome_name <- species$GENOME[sp]
      common_name <- species$COMMON[sp]
      
      if (length(grep(genome_name, unlist(file_i))) > 0) {
        matches <- grep(genome_name, unlist(file_i))
        order_to_move <- c(order_to_move, matches, matches + 1)
        
        # Replace genome name with common name
        for (col_idx in 1:ncol(file_i)) {
          file_i[[col_idx]] <- stringr::str_replace(file_i[[col_idx]], 
                                                   pattern = genome_name,
                                                   replacement = common_name)
        }
      }
    }
    
    # Reorder if needed
    if (length(order_to_move) > 0) {
      file_i <- as.data.frame(file_i[order_to_move, ])
    }
    
    write.table(file_i, file_path, row.names = FALSE, col.names = FALSE, 
                quote = FALSE)
    processed <- processed + 1
  }
  
  invisible(processed)
}

#' Find overlaps between smORF CDS and GWAS mutations
#'
#' @param gtf_path path to GTF file
#' @param gwas_file path to GWAS catalog TSV file
#' @param output_path path to output CSV file
#' @return data.frame with overlap results
find_gwas_overlaps <- function(gtf_path, gwas_file, output_path = NULL) {
  suppressed_require("rtracklayer")
  suppressed_require("GenomicRanges")
  suppressed_require("data.table")
  
  # Read GWAS catalog
  GWAS <- data.table::fread(gwas_file)
  GWAS$CHR_POS <- as.numeric(GWAS$CHR_POS)
  GWAS <- GWAS[!is.na(GWAS$CHR_POS), ]
  
  # Read GTF and filter to orfCDS
  GTF <- rtracklayer::import(gtf_path)
  GTF <- GTF[GTF$type == "orfCDS", ]
  
  # Create GRanges objects
  smORF_granges <- GenomicRanges::GRanges(
    IRanges::IRanges(start = GTF$start, end = GTF$end),
    seqnames = GTF$seqnames,
    iORF_ID = GTF$iORF_id
  )
  
  GWAS_granges <- GenomicRanges::GRanges(
    IRanges::IRanges(start = GWAS$CHR_POS, end = GWAS$CHR_POS),
    seqnames = GWAS$CHR_ID,
    SNP_ID = GWAS$SNPS
  )
  
  # Find overlaps
  smORF_granges_list <- GenomicRanges::GNCList(smORF_granges)
  overlap <- GenomicRanges::findOverlaps(GWAS_granges, smORF_granges_list)
  
  # Combine results
  GWAS_results <- cbind(
    as.data.frame(GTF[overlap@to, ]),
    as.data.frame(GWAS[overlap@from, ])
  )
  
  # Save if output path provided
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(GWAS_results, output_path, row.names = FALSE)
  }
  
  return(GWAS_results)
}

#' Process GWAS mutations to check for amino acid changes
#' 
#' Note: This requires smORF objects with nucleotide sequences, which may not be available.
#' This function is optional and may not be used in all workflows.
#'
#' @param gwas_overlaps data.frame from find_gwas_overlaps
#' @param smorf_objects_dir directory containing smORF RDS objects
#' @param output_path path to output CSV file
#' @return data.frame with mutation analysis
process_gwas_mutations <- function(gwas_overlaps, smorf_objects_dir, output_path = NULL) {
  suppressed_require("seqinr")
  suppressed_require("stringr")
  
  if (!dir.exists(smorf_objects_dir)) {
    warning("smORF objects directory not found: ", smorf_objects_dir, 
            ". Skipping mutation processing.")
    return(gwas_overlaps)
  }
  
  gwas_overlaps$peptide_seq <- NA
  gwas_overlaps$seq_after_mut <- NA
  gwas_overlaps$seq_change <- NA
  
  for (i in 1:nrow(gwas_overlaps)) {
    iorf_id <- gwas_overlaps$iORF_id[i]
    smorf_file <- file.path(smorf_objects_dir, paste0(iorf_id, ".rds"))
    
    if (!file.exists(smorf_file)) {
      next
    }
    
    smORF_object <- readRDS(smorf_file)
    nt_seq <- smORF_object$nt_seq
    nt_seq <- gsub("-", "", nt_seq)
    nt_seq <- toupper(nt_seq)
    original_AA_seq <- seqinr::c2s(seqinr::translate(seqinr::s2c(nt_seq)))
    gwas_overlaps$peptide_seq[i] <- original_AA_seq
    
    smORF_i_gtf <- smORF_object$gtf
    GWAS_i <- gwas_overlaps[i, ]
    change <- "No"
    
    snp <- stringr::str_sub(GWAS_i$STRONGEST.SNP.RISK.ALLELE,
                            start = nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE),
                            end = nchar(GWAS_i$STRONGEST.SNP.RISK.ALLELE))
    
    if (snp %in% c("A", "C", "G", "T")) {
      pos_snp <- GWAS_i$CHR_POS
      
      if (nrow(smORF_i_gtf) == 1) {
        if (smORF_i_gtf$strand == "+") {
          pos_relative_snp <- pos_snp - smORF_i_gtf$start + 1
          substr(nt_seq, pos_relative_snp, pos_relative_snp) <- snp
        } else if (smORF_i_gtf$strand == "-") {
          pos_relative_snp <- smORF_i_gtf$end - pos_snp + 1
          substr(nt_seq, pos_relative_snp, pos_relative_snp) <- chartr("ATGC", "TACG", snp)
        }
      }
      
      new_seq <- seqinr::c2s(seqinr::translate(seqinr::s2c(nt_seq)))
      gwas_overlaps$seq_after_mut[i] <- new_seq
      
      if (new_seq != original_AA_seq) {
        change <- "Yes"
      }
    }
    
    gwas_overlaps$seq_change[i] <- change
  }
  
  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(gwas_overlaps, output_path, row.names = FALSE)
  }
  
  return(gwas_overlaps)
}
