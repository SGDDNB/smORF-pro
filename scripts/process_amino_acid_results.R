## Helper functions for amino acid analysis processing
## These functions encapsulate the code currently shown in README / README.Rmd
## so they can be reused by scripts and Shiny apps.

suppressed_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Please install it first.", pkg), call. = FALSE)
  }
}

#' Create FASTA file from Annotation data.frame
#'
#' @param annotation data.frame with at least columns iORFID and Peptide_seq
#' @param fasta_path output FASTA file path
#' @return invisible fasta_path
create_fasta_from_annotation <- function(annotation, fasta_path) {
  suppressed_require("seqinr")

  required_cols <- c("iORFID", "Peptide_seq")
  missing <- setdiff(required_cols, colnames(annotation))
  if (length(missing) > 0) {
    stop("Annotation is missing required columns: ", paste(missing, collapse = ", "))
  }

  seqinr::write.fasta(
    sequences = as.list(annotation$Peptide_seq),
    names = annotation$iORFID,
    as.string = FALSE,
    file.out = fasta_path
  )

  invisible(fasta_path)
}


#' Basic file-format validation for uploaded tool outputs
#'
#' @param path path to uploaded file
#' @param tool_name one of \"deeptmhmm\", \"targetp\", \"deeploc\", \"interpro\"
#' @return TRUE if passes basic checks, otherwise stops with error
validate_file_format <- function(path, tool_name) {
  if (is.null(path) || !nzchar(path)) {
    stop("No file selected for ", tool_name, ".", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))

  expected_ext <- switch(
    tool_name,
    deeptmhmm = c("3line", "txt"),
    targetp   = c("txt"),
    deeploc   = c("csv"),
    interpro  = c("tsv", "txt"),
    character()
  )

  if (length(expected_ext) && !ext %in% expected_ext) {
    stop(
      sprintf(
        "Unexpected file extension '%s' for %s. Expected: %s",
        ext, tool_name, paste(expected_ext, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  TRUE
}


#' Process DeepTMHMM 3-line format output and add deepTMHMM column
#'
#' @param path path to DeepTMHMM 3line output
#' @param annotation Annotation data.frame
#' @return updated Annotation data.frame
process_deeptmhmm <- function(path, annotation) {
  suppressed_require("stringr")

  validate_file_format(path, "deeptmhmm")

  file_3line <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  file_3line <- file_3line[substr(file_3line$V1, 1, 1) == ">", , drop = FALSE]
  file_3line <- as.data.frame(sub(">", "", file_3line$V1), stringsAsFactors = FALSE)
  colnames(file_3line) <- "raw"

  split <- stringr::str_split_fixed(file_3line$raw, " \\| ", 2)
  file_3line$iORFID <- split[, 1]
  file_3line$TMHMM  <- split[, 2]

  annotation$deepTMHMM <- file_3line$TMHMM[match(annotation$iORFID, file_3line$iORFID)]
  annotation
}


#' Process TargetP 2.0 prediction summary and add TargetP column
#'
#' @param path path to TargetP prediction summary txt
#' @param annotation Annotation data.frame
#' @return updated Annotation data.frame
process_targetp <- function(path, annotation) {
  suppressed_require("data.table")

  validate_file_format(path, "targetp")

  targetP <- data.table::fread(path, data.table = FALSE)
  # Expect columns including \"# ID\" and \"Prediction\"
  if (!all(c("# ID", "Prediction") %in% colnames(targetP))) {
    stop("TargetP file does not contain expected columns '# ID' and 'Prediction'.", call. = FALSE)
  }

  annotation$TargetP <- targetP$Prediction[match(annotation$iORFID, targetP[["# ID"]])]
  annotation
}


#' Process Deeploc CSV summary and add Deeploc column
#'
#' @param path path to Deeploc CSV summary
#' @param annotation Annotation data.frame
#' @return updated Annotation data.frame
process_deeploc <- function(path, annotation) {
  validate_file_format(path, "deeploc")

  deeploc <- utils::read.csv(path, stringsAsFactors = FALSE)

  expected_cols <- c("Protein_ID", "Localizations")
  missing <- setdiff(expected_cols, colnames(deeploc))
  if (length(missing) > 0) {
    stop("Deeploc file is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  annotation$Deeploc <- deeploc$Localizations[match(annotation$iORFID, deeploc$Protein_ID)]
  annotation
}


#' Process Interproscan TSV output and add interproscan summary column
#'
#' @param path path to Interproscan TSV summary
#' @param annotation Annotation data.frame
#' @return updated Annotation data.frame
process_interproscan <- function(path, annotation) {
  suppressed_require("data.table")

  validate_file_format(path, "interpro")

  interpro <- data.table::fread(path, data.table = FALSE)

  if (ncol(interpro) < 15) {
    stop("Interproscan TSV file has fewer than 15 columns; please check you downloaded the correct format.", call. = FALSE)
  }

  colnames(interpro)[1:15] <- c(
    "smORF_ID", "md5", "length", "Analysis", "Signature_accession",
    "Signature_description", "start", "stop", "score", "Status", "Date",
    "interpro_accession", "interpro_description", "GO_Annotation", "Pathway_Annotation"
  )

  smORFs_with_signals <- unique(interpro$smORF_ID)

  annotation$interproscan <- "No domain found"

  for (smORF_i in smORFs_with_signals) {
    interpro_i <- interpro[interpro$smORF_ID == smORF_i, , drop = FALSE]
    signals_i <- unique(c(unlist(interpro_i[, c("Signature_description", "interpro_description")], use.names = FALSE)))
    signals_i <- signals_i[signals_i != "-"]
    if (length(signals_i) == 0) next
    annotation$interproscan[annotation$iORFID == smORF_i] <- paste(signals_i, collapse = "; ")
  }

  annotation
}


#' Convenience function to integrate all available amino acid analysis results
#'
#' Each argument is an optional file path; if NULL or empty, that step is skipped.
#'
#' @param annotation Annotation data.frame
#' @param deeptmhmm_path optional DeepTMHMM 3line file
#' @param targetp_path optional TargetP summary file
#' @param deeploc_path optional Deeploc CSV summary
#' @param interpro_path optional Interproscan TSV file
#' @return updated Annotation data.frame
integrate_amino_acid_results <- function(annotation,
                                         deeptmhmm_path = NULL,
                                         targetp_path = NULL,
                                         deeploc_path = NULL,
                                         interpro_path = NULL) {
  out <- annotation

  if (!is.null(deeptmhmm_path) && nzchar(deeptmhmm_path)) {
    out <- process_deeptmhmm(deeptmhmm_path, out)
  }
  if (!is.null(targetp_path) && nzchar(targetp_path)) {
    out <- process_targetp(targetp_path, out)
  }
  if (!is.null(deeploc_path) && nzchar(deeploc_path)) {
    out <- process_deeploc(deeploc_path, out)
  }
  if (!is.null(interpro_path) && nzchar(interpro_path)) {
    out <- process_interproscan(interpro_path, out)
  }

  out
}

