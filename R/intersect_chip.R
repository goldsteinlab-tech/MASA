#' Run bedtoolsr intersect on peak files with reference using bedtoolsr
#'
#' This function runs `bedtoolsr::bt.intersect` for each peak file in a specified directory
#' against a reference BED file (e.g., H3K27ac ChIP-seq peaks). Each result is saved in the
#' output directory using a customizable result name suffix.
#'
#' @param atac_dir Path to the directory containing peak files (e.g., ATAC-seq)
#' @param chip_file Path to the reference chip-seq BED file
#' @param results_dir Path to the directory where result `.bed` files will be stored
#' @param result_name Suffix used in the output file names
#' @param file_extension File extension to match (e.g., "bed", "narrowPeak"). Do not include the dot.
#'
#' @return This function does not return an R object. It writes `.bed` files to the specified output directory.
#' @return Use -wa and -u as intersect option to keep only the unique ATAc-seq peaks intersected with the chip-seq peaks
#' @export
intersect_atac_with_chip <- function(
  atac_dir,
  chip_file,
  results_dir,
  result_name,
  file_extension = "narrowPeak"
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  # Construct the file pattern from the extension
  pattern <- paste0("\\.", file_extension, "$")
  
  # List all matching peak files
  atac_files <- list.files(
    path = atac_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  # Load the reference ChIP-seq peak file
  chip_bed <- read.table(chip_file, header = FALSE, sep = "\t")
  
  # intersect each peak file with the reference
  for (atac_file in atac_files) {
    base <- tools::file_path_sans_ext(basename(atac_file))
    out_file <- file.path(results_dir, paste0(base, "_", result_name, ".bed"))
    
    atac_bed <- read.table(atac_file, header = FALSE, sep = "\t")
    
    intersected <- bedtoolsr::bt.intersect(a = atac_bed, b = chip_bed, wa = TRUE, u = TRUE)
    
    write.table(intersected, file = out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
