#' Convert wwarn k13 SNP formats for Stave
#'
#' @details
#' 
#' This function processes k13 SNP data into a format compatible with the `STAVE` tool. 
#'
#' NB FOR STAVE: "|" reflects phased 
#' 
#' NB FOR STAVE "/" reflects unphased 
#'
#' ## Input format conventions:
#' - `_`: Denotes separate, distinct haplotypes (different clones).
#' - `/`: Denotes additional non-synonymous (NS) changes in the same clone.
#'
#' ## Output format:
#' - `WT`: Wild type.
#' - Mutations are prefixed with "k13-" and standardised.
#'
#' @param marker Character. A string representing k13 SNP data in various formats.
#' @return Character. Processed k13 SNP data in a standardised format for `STAVE`.
#' Returns `NA` if the input is missing or errors if unexpected type
#' @export
#' @examples
#' wwarn_format_k13_for_stave("")  # Blank
#' wwarn_format_k13_for_stave("WT")  # WT
#' wwarn_format_k13_for_stave("P441L")  # Single SNP homozygous
#' wwarn_format_k13_for_stave("R528G/T")  # Mixed Infection
#' wwarn_format_k13_for_stave("T474I_V520I")  # Phased double mutant
wwarn_format_k13_for_stave <- function(marker) {
  
  # Case 1: NA or "-" - Missing
  if (is.na(marker) || marker == "-") {
    return(NA)
  }
  
  # Case 2: Empty string - Wild Type
  if (marker == "WT") {
    return("k13:WT")
  }
  
  # Case 3: Simple SNP format (no `/` or `,`)
  if (!grepl("\\/|\\,|\\_|WT", marker)) {
    return(clean_mutations(paste0("k13-", gsub("*", "", toupper(marker), fixed = TRUE))))
  }
  
  # Case 4: Unphased at same locus (just / and no _)
  if (grepl("\\/", marker) && !grepl("\\_", marker)) {
    marker <- gsub("(\\w)(.*)", "\\2", marker)
    loci <- gsub("(\\d*)(.*)", "\\1", marker)
    snps <- gsub("(\\d*)(.*)", "\\2", marker)
    return(paste0("k13:", loci, ":", snps))
  }
  
  # Case 5: phased and on same clone (_ and no /)
  if (grepl("\\_", marker) && !grepl("\\/", marker)) {
    spl <- toupper(strsplit(marker, "\\_")[[1]])
    all_loci <- gsub("(^\\w)(\\d*)(\\w)", "\\2", spl)
    all_snps <- gsub("(^\\w)(\\d*)(\\w)", "\\3", spl)
    res <- paste0("k13:",
                  paste0(sort(as.numeric(all_loci)), collapse = "_"),
                  ":",
                  paste0(all_snps[order(as.numeric(all_loci))], collapse = "_"))
    return(res)
  }
  
  stop("Unexpected mutation")
  
}
