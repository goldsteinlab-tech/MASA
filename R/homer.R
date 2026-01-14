#' Run HOMER Motif Analysis
#'
#' Runs HOMER's motif discovery tool (`findMotifsGenome.pl`) on a set of peaks and optionally processes the results.
#'
#' @param input_file Character. Path to the input peak file (BED format).
#' @param output_dir Character. Directory to save HOMER results (will be created if needed).
#' @param genome Character. Genome identifier (e.g., "hg38", "mm10"). Default is "mm10".
#' @param homer_path Character. Path to HOMER's findMotifsGenome.pl script (default: "findMotifsGenome.pl" in PATH).
#' @param additional_args Character. Additional arguments for HOMER as a single string (e.g., "-len 10 -size 200").
#' @param process_results Logical. If TRUE, process and return results as R objects. Default is TRUE.
#'
#' @return If \code{process_results} is TRUE, returns a list of processed HOMER results; otherwise, returns NULL invisibly.
#' @export
#'
#' @examples
#' # run_homer("peaks.bed", "homer_output", genome = "hg38", additional_args = "-len 10 -size 200")
run_homer <- function(
  input_file,
  output_dir,
  genome = "mm10",
  homer_path = "findMotifsGenome.pl",
  additional_args = NULL,
  process_results = TRUE
) {
  # Check if the input file exists
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Build the HOMER command
  command <- paste(
    shQuote(homer_path),
    shQuote(input_file),
    genome,
    shQuote(output_dir)
  )
  if (!is.null(additional_args)) {
    command <- paste(command, additional_args)
  }
  message("Running HOMER with the following command:\n", command)
  system(command, intern = FALSE)
  # Optionally process results
  if (process_results) {
    return(process_homer_results(output_dir))
  } else {
    invisible(NULL)
  }
}

#' Process HOMER Results
#'
#' Loads and parses key output files from a HOMER motif analysis directory.
#'
#' @param output_dir Character. Directory containing HOMER output.
#'
#' @return List containing key results: known motif results, motif definitions, parameters, HTML summary, and motif logos.
#' @export
#'
#' @examples
#' # process_homer_results("homer_output")
process_homer_results <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  files <- list.files(output_dir, full.names = TRUE)
  homer_results <- list()
  for (file in files) {
    if (grepl("knownResults.txt$", file)) {
      homer_results$known_results <- read.delim(file, header = TRUE, sep = "\t")
    } else if (grepl("homerMotifs.all.motifs$", file)) {
      homer_results$motifs <- readLines(file)
    } else if (grepl("motifFindingParameters.txt$", file)) {
      homer_results$parameters <- readLines(file)
    } else if (grepl("homerResults.html$", file)) {
      homer_results$html_summary <- file
    } else if (grepl("logo", file)) {
      if (is.null(homer_results$logos)) homer_results$logos <- list()
      homer_results$logos[[basename(file)]] <- file
    }
  }
  return(homer_results)
}

