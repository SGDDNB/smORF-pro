## Helper functions for managing assets folder (downloaded data files)
## These functions scan, validate, and provide information about files in assets/

#' Get the base path to assets folder
#' 
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return path to assets folder
get_assets_base <- function(from_apps = FALSE) {
  if (from_apps) {
    return(normalizePath("../assets", winslash = "/", mustWork = FALSE))
  } else {
    return(normalizePath("assets", winslash = "/", mustWork = FALSE))
  }
}

#' Scan assets subfolder and return list of files with metadata
#'
#' @param asset_type one of "maf_files", "gwas_catalog", "species_names", "clustalo"
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return data.frame with columns: file, path, size_bytes, size_mb, modified, exists
scan_assets_folder <- function(asset_type, from_apps = FALSE) {
  valid_types <- c("maf_files", "gwas_catalog", "species_names", "clustalo")
  if (!asset_type %in% valid_types) {
    stop("asset_type must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  assets_base <- get_assets_base(from_apps)
  folder_path <- file.path(assets_base, asset_type)
  
  if (!dir.exists(folder_path)) {
    result_df <- data.frame(
      file = character(0),
      path = character(0),
      size_bytes = numeric(0),
      size_mb = numeric(0),
      modified = character(0),
      exists = logical(0),
      stringsAsFactors = FALSE
    )
    if (asset_type == "maf_files") {
      result_df$is_gzipped <- logical(0)
    }
    return(result_df)
  }
  
  # For MAF files, look for both .maf and .maf.gz
  if (asset_type == "maf_files") {
    files <- list.files(folder_path, pattern = "\\.maf(\\.gz)?$", full.names = TRUE)
  } else {
    files <- list.files(folder_path, full.names = TRUE, recursive = FALSE)
  }
  
  if (length(files) == 0) {
    result_df <- data.frame(
      file = character(0),
      path = character(0),
      size_bytes = numeric(0),
      size_mb = numeric(0),
      modified = character(0),
      exists = logical(0),
      stringsAsFactors = FALSE
    )
    if (asset_type == "maf_files") {
      result_df$is_gzipped <- logical(0)
    }
    return(result_df)
  }
  
  # Filter to files only (not directories)
  files <- files[file.info(files)$isdir == FALSE]
  
  if (length(files) == 0) {
    result_df <- data.frame(
      file = character(0),
      path = character(0),
      size_bytes = numeric(0),
      size_mb = numeric(0),
      modified = character(0),
      exists = logical(0),
      stringsAsFactors = FALSE
    )
    if (asset_type == "maf_files") {
      result_df$is_gzipped <- logical(0)
    }
    return(result_df)
  }
  
  file_info <- file.info(files)
  result <- data.frame(
    file = basename(files),
    path = files,
    size_bytes = file_info$size,
    size_mb = round(file_info$size / 1024 / 1024, 2),
    modified = format(file_info$mtime, "%Y-%m-%d %H:%M:%S"),
    exists = file_info$size > 0,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Add is_gzipped column for MAF files
  if (asset_type == "maf_files") {
    result$is_gzipped <- grepl("\\.gz$", result$file)
  }
  
  # Sort by modified date (newest first)
  result <- result[order(result$modified, decreasing = TRUE), ]
  
  return(result)
}

#' Check for available MAF files
#'
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return data.frame with file information (from scan_assets_folder)
check_maf_files <- function(from_apps = FALSE) {
  scan_assets_folder("maf_files", from_apps)
}

#' Check for available GWAS catalog files
#'
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return data.frame with file information (from scan_assets_folder)
check_gwas_catalog <- function(from_apps = FALSE) {
  scan_assets_folder("gwas_catalog", from_apps)
}

#' Get file information (size, date, etc.)
#'
#' @param file_path path to file
#' @return list with size_bytes, size_mb, modified, exists
get_file_info <- function(file_path) {
  if (is.null(file_path) || !nzchar(file_path) || !file.exists(file_path)) {
    return(list(
      size_bytes = 0,
      size_mb = 0,
      modified = NA_character_,
      exists = FALSE
    ))
  }
  
  info <- file.info(file_path)
  list(
    size_bytes = info$size,
    size_mb = round(info$size / 1024 / 1024, 2),
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    exists = TRUE
  )
}

#' Basic validation of MAF file format
#'
#' @param file_path path to MAF file
#' @return TRUE if passes basic checks, otherwise stops with error
validate_maf_file <- function(file_path) {
  if (is.null(file_path) || !nzchar(file_path)) {
    stop("No MAF file specified.", call. = FALSE)
  }
  if (!file.exists(file_path)) {
    stop("MAF file does not exist: ", file_path, call. = FALSE)
  }
  
  # Check file extension
  ext <- tolower(tools::file_ext(file_path))
  if (ext != "maf") {
    warning("File does not have .maf extension: ", file_path)
  }
  
  # Check if file is readable and has content
  info <- file.info(file_path)
  if (info$size == 0) {
    stop("MAF file is empty: ", file_path, call. = FALSE)
  }
  
  # Try to read first few lines to check format
  con <- file(file_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  
  if (length(first_line) == 0) {
    stop("MAF file appears to be empty or unreadable: ", file_path, call. = FALSE)
  }
  
  # MAF files typically start with "##maf" or "a" (alignment block)
  if (!grepl("^##maf|^a\\s|^#", first_line, ignore.case = TRUE)) {
    warning("MAF file may not be in correct format. First line: ", substr(first_line, 1, 50))
  }
  
  return(TRUE)
}

#' Basic validation of GWAS catalog file format
#'
#' @param file_path path to GWAS catalog TSV file
#' @return TRUE if passes basic checks, otherwise stops with error
validate_gwas_file <- function(file_path) {
  if (is.null(file_path) || !nzchar(file_path)) {
    stop("No GWAS catalog file specified.", call. = FALSE)
  }
  if (!file.exists(file_path)) {
    stop("GWAS catalog file does not exist: ", file_path, call. = FALSE)
  }
  
  # Check file extension
  ext <- tolower(tools::file_ext(file_path))
  if (!ext %in% c("tsv", "txt", "csv")) {
    warning("GWAS catalog file does not have expected extension (.tsv, .txt, or .csv): ", file_path)
  }
  
  # Check if file is readable and has content
  info <- file.info(file_path)
  if (info$size == 0) {
    stop("GWAS catalog file is empty: ", file_path, call. = FALSE)
  }
  
  # GWAS catalog files should be reasonably large (300+ MB)
  if (info$size < 100 * 1024 * 1024) {  # Less than 100 MB
    warning("GWAS catalog file seems small (", round(info$size / 1024 / 1024, 2), 
            " MB). Expected ~300+ MB for full catalog.")
  }
  
  # Try to read first line to check if it's TSV
  con <- file(file_path, "r")
  first_line <- readLines(con, n = 1)
  close(con)
  
  if (length(first_line) == 0) {
    stop("GWAS catalog file appears to be empty or unreadable: ", file_path, call. = FALSE)
  }
  
  # Check if it looks like TSV (has tabs)
  if (!grepl("\t", first_line)) {
    warning("GWAS catalog file may not be TSV format (no tabs found in first line).")
  }
  
  return(TRUE)
}

#' Check if Clustal Omega tool is available
#'
#' @return logical, TRUE if clustalo is found in PATH
check_clustalo_available <- function() {
  result <- tryCatch({
    system2("clustalo", "--version", stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  
  if (is.logical(result)) {
    return(result)
  }
  
  # If we got output, tool is available
  return(TRUE)
}

#' Get species names file path
#'
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @return path to species_names.txt file
get_species_names_path <- function(from_apps = FALSE) {
  assets_base <- get_assets_base(from_apps)
  default_path <- file.path(assets_base, "species_names", "species_names.txt")
  
  # Also check repo version as fallback
  if (from_apps) {
    repo_path <- normalizePath("../data_preparation/Genetic_information/species_names.txt", 
                              winslash = "/", mustWork = FALSE)
  } else {
    repo_path <- normalizePath("data_preparation/Genetic_information/species_names.txt", 
                              winslash = "/", mustWork = FALSE)
  }
  
  # Prefer assets version, fallback to repo version
  if (file.exists(default_path)) {
    return(default_path)
  } else if (file.exists(repo_path)) {
    return(repo_path)
  } else {
    return(default_path)  # Return expected path even if doesn't exist
  }
}

#' Download MAF files from UCSC
#'
#' @param chromosomes vector of chromosome names to download (e.g., c("chr1", "chr2", ...))
#'                   If NULL, downloads all standard chromosomes (1-22, X, Y, M)
#' @param output_dir directory to save MAF files (default: assets/maf_files)
#' @param from_apps logical, if TRUE assumes called from apps/ directory
#' @param base_url base URL for UCSC MAF files
#' @param progress_callback optional function(file_num, total_files, filename) called for each file
#' @return list with downloaded files and any errors
download_maf_files <- function(chromosomes = NULL, output_dir = NULL, 
                               from_apps = FALSE, 
                               base_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/",
                               progress_callback = NULL) {
  
  if (is.null(output_dir)) {
    assets_base <- get_assets_base(from_apps)
    output_dir <- file.path(assets_base, "maf_files")
  }
  
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Default chromosomes if not specified
  if (is.null(chromosomes)) {
    chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  }
  
  downloaded <- character(0)
  errors <- character(0)
  skipped <- character(0)
  
  # Filter out chromosomes that already have files
  files_to_download <- character(0)
  for (chr in chromosomes) {
    filename <- paste0(chr, ".maf")
    filepath <- file.path(output_dir, filename)
    if (!file.exists(filepath)) {
      files_to_download <- c(files_to_download, chr)
    } else {
      skipped <- c(skipped, filename)
    }
  }
  
  total_files <- length(files_to_download)
  
  if (total_files == 0) {
    return(list(
      downloaded = character(0),
      errors = character(0),
      skipped = skipped,
      total_requested = length(chromosomes),
      total_downloaded = 0,
      total_skipped = length(skipped)
    ))
  }
  
  # Check which tool is available (curl or wget)
  has_curl <- nzchar(Sys.which("curl"))
  has_wget <- nzchar(Sys.which("wget"))
  
  if (!has_curl && !has_wget) {
    stop("Neither 'curl' nor 'wget' is available. Please install one of them to download files.")
  }
  
  # Use curl if available (more common on macOS), otherwise wget
  use_curl <- has_curl
  
  file_num <- 0
  for (chr in files_to_download) {
    file_num <- file_num + 1
    filename <- paste0(chr, ".maf")
    filepath <- file.path(output_dir, filename)
    url <- paste0(base_url, filename)
    
    # Call progress callback if provided
    if (!is.null(progress_callback)) {
      progress_callback(file_num, total_files, filename)
    }
    
    tryCatch({
      if (use_curl) {
        # Use curl with progress bar
        cmd <- sprintf('curl -L -o "%s" "%s"', filepath, url)
        result <- system(cmd, intern = TRUE)
      } else {
        # Use wget with progress bar
        cmd <- sprintf('wget -O "%s" "%s"', filepath, url)
        result <- system(cmd, intern = TRUE)
      }
      
      # Check if file was downloaded successfully
      if (file.exists(filepath) && file.info(filepath)$size > 0) {
        downloaded <- c(downloaded, filename)
      } else {
        errors <- c(errors, paste("Failed to download", filename))
      }
    }, error = function(e) {
      errors <- c(errors, paste("Error downloading", filename, ":", e$message))
    })
  }
  
  list(
    downloaded = downloaded,
    errors = errors,
    skipped = skipped,
    total_requested = length(chromosomes),
    total_downloaded = length(downloaded),
    total_skipped = length(skipped)
  )
}

#' Get chromosomes needed from GTF file
#'
#' @param gtf_path path to GTF file
#' @return vector of unique chromosome names found in GTF
get_chromosomes_from_gtf <- function(gtf_path) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    warning("rtracklayer package not available. Cannot determine chromosomes from GTF.")
    return(NULL)
  }
  
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    warning("GenomicRanges package not available. Cannot determine chromosomes from GTF.")
    return(NULL)
  }
  
  tryCatch({
    gtf <- rtracklayer::import(gtf_path)
    chromosomes <- unique(as.character(GenomicRanges::seqnames(gtf)))
    # Filter to standard chromosome format (chr1, chr2, etc.)
    chromosomes <- chromosomes[grepl("^chr", chromosomes)]
    return(chromosomes)
  }, error = function(e) {
    warning("Could not read GTF to determine chromosomes: ", e$message)
    return(NULL)
  })
}
