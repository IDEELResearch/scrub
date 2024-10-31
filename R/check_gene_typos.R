#' Check for Gene Name Typos
#'
#' This function checks a column of gene names for common typos and prints warnings 
#' when it finds any. It searches for probable misspellings of "k13", "crt", and "mdr1".
#'
#' @param gene_column A character vector containing gene names to check for typos.
#' @return Prints a warning message for each typo found, indicating the probable 
#'         correct gene name. No value is returned.
#' @examples
#' genes <- c("kletch13", "ctr", "mdr1", "mdr", "kelch 13")
#' check_gene_typos(genes)
#' # Output:
#' # Probable typo found for k13: did you mean k13?
#' # Probable typo found for crt: did you mean crt?
#' # Probable typo found for mdr1: did you mean mdr1?
#'
#' 
check_gene_typos <- function(gene_column) {
  # Define probable typos for k13, crt, and mdr1
  probable_typos <- list(
    "k13" = c("kelch13", "kelch 13", "kelch_13", "kelch", "kletch13", "klech 13"),
    "crt" = c("ctr"),
    "mdr1" = c("mrd1", "mdr")
  )
  
  for (correct_gene in names(probable_typos)) {
    typos_found <- gene_column[gene_column %in% probable_typos[[correct_gene]]]
    if (length(typos_found) > 0) {
      message(paste("Probable typo found for", correct_gene, ": did you mean", correct_gene, "?"))
    }
  }
}
