#' Convert pf7k k13 SNP formats for Stave
#'
#' @details
#' 
#' This function processes k13 SNP data into a format compatible with the `STAVE` tool. 
#' It handles different scenarios of input data, such as wild types, phased/unphased mutations, 
#' and distinct haplotypes. The function accommodates variations in notation, including the use 
#' of `,`, `/`, and `*` in the input string.
#'
#' NB FOR STAVE: "|" reflects phased 
#' 
#' NB FOR STAVE "/" reflects unphased 
#'
#' ## Input format conventions:
#' - **Uppercase letters**: Homozygous mutations.
#' - **Lowercase letters**: Heterozygous mutations.
#' - **Lowercase without a second mutation**: Heterozygous at the SNP level but same amino acid.
#' - `,`: Denotes separate, distinct haplotypes (different clones).
#' - `/`: Denotes additional non-synonymous (NS) changes in the same clone.
#' - `*`: Indicates unphased data for the sample.
#'
#' ## Output format:
#' - `WT`: Wild type.
#' - Mutations are prefixed with "k13-" and standardised.
#' - Multiple haplotypes or NS changes are separated with `&&` for downstream splitting.
#'
#' @param marker Character. A string representing k13 SNP data in various formats.
#' @return Character. Processed k13 SNP data in a standardised format for `STAVE`.
#' Returns `NA` if the input is missing.
#' @export
#' @examples
#' pf7k_format_k13_for_stave("")  # Blank
#' pf7k_format_k13_for_stave("C580Y")  # Single SNP homozygous
#' pf7k_format_k13_for_stave("c580y")  # Single SNP heterozygous but same SNP
#' pf7k_format_k13_for_stave("C580Y,G449A")  # Distinct haplotypes
#' pf7k_format_k13_for_stave("C580Y/G449A")  # NS changes in same clone
#' pf7k_format_k13_for_stave("C580Y,G449A/I543T")  # Mixed haplotypes and NS changes
pf7k_format_k13_for_stave <- function(marker) {
  
  # Case 1: NA or "-" - Missing
  if (is.na(marker) || marker == "-") {
    return(NA)
  }
  
  # Case 2: Empty string - Wild Type
  if (marker == "") {
    return("WT")
  }
  
  # Case 3: Simple SNP format (no `/` or `,`)
  if (!grepl("\\/|\\,", marker)) {
    return(clean_mutations(paste0("k13-", gsub("*", "", toupper(marker), fixed = TRUE))))
  }
  
  # Case 4: Multiple distinct haplotypes separated by `,`
  if (grepl("\\,", marker) && !grepl("\\/", marker)) {
    marker <- gsub("*", "", marker, fixed = TRUE)
    spl <- toupper(strsplit(marker, "\\,|\\/")[[1]])
    res <- vapply(paste0("k13-", spl), clean_mutations, character(1))
    return(paste0(as.character(res), collapse = "&&"))
  }
  
  # Case 5: Non-synonymous changes in the same clone (separated by `/`)
  if (grepl("\\/", marker) && !grepl("\\,", marker)) {
    marker <- gsub("*", "", marker, fixed = TRUE)
    spl <- toupper(strsplit(marker, "\\,|\\/")[[1]])
    all_loci <- unique(gsub("(^\\w)(\\d*)(\\w)", "\\2", spl))
    all_snps <- unique(gsub("(^\\w)(\\d*)(\\w)", "\\3", spl))
    res <- paste0("k13:",
                  paste0(sort(as.numeric(all_loci)), collapse = "_"),
                  ":",
                  paste0(all_snps[order(as.numeric(all_loci))], collapse = "_"))
    return(res)
  }
  
  # Case 6: Mixed haplotypes and NS changes (`/` and `,`)
  if (grepl("\\,", marker) && grepl("\\/", marker)) {
    marker <- gsub("*", "", marker, fixed = TRUE)
    spl <- toupper(strsplit(marker, "\\,")[[1]])
    spl2 <- strsplit(spl, "\\/")
    all_loci <- lapply(spl2, function(x) unique(gsub("(^\\w)(\\d*)(\\w)", "\\2", x)))
    all_snps <- lapply(spl2, function(x) unique(gsub("(^\\w)(\\d*)(\\w)", "\\3", x)))
    res <- lapply(seq_along(all_loci), function(x) {
      paste0("k13:",
             paste0(sort(as.numeric(all_loci[[x]])), collapse = "_"),
             ":",
             paste0(all_snps[[x]][order(as.numeric(all_loci[[x]]))], collapse = "_"))
    })
    return(paste0(as.character(unlist(res)), collapse = "&&"))
  }
  
}


#' Convert Biallelic SNP Formats for Stave
#'
#' This function processes biallelic SNP data and converts it into a format compatible 
#' with the `stave` tool. It ensures compatibility by handling phased, unphased, 
#' and missing data, as well as removing redundant SNPs.
#'
#' ## Input format conventions:
#' - **Phased data**: Denoted by `,`, converted to `|` in the output.
#' - **Unphased data**: Marked with `*`, adjusted to `/` for compatibility.
#' - **Missing data**: Represented by `-`, converted to `NA`.
#' - **Redundant SNPs**: Identical alleles (e.g., `G|G`) are simplified.
#'
#' ## Output format:
#' - SNPs are prefixed with the gene locus in the form `gene:locus:SNP`.
#' - Phased SNPs are joined with `|`.
#' - Unphased SNPs are joined with `/`.
#'
#' @param mut Character vector. A vector of SNP mutations in a biallelic format.
#' @param gls Character. The gene locus and SNP identifier, inferred from the input variable name by default.
#' @return Character vector. Processed SNP data in a standardised format for `stave`. Missing values are returned as `NA`.
#' @export
#' @examples
#' 
#' # Example input with phased and unphased SNPs
#' crt_76 <- c("K,T", "K*,T", "T", "-") 
#' pf7k_format_bia_for_stave(crt_76)  
#' 
#' # Handles redundant and unphased SNPs
#' pf7k_format_bia_for_stave(c("G,G", "A,C", "A*|G"), gls = "dhps_581_a")  
pf7k_format_bia_for_stave <- function(mut, gls = deparse(substitute(mut))) {
  
  # Extract gene locus and SNP identifier
  gl <- gsub("(^[[:alnum:]]*_\\d*)(_[[:alnum:]]$)", "\\1", gls)
  gl <- gsub("_", ":", gl)
  
  # Convert to stave phased format
  res <- gsub("\\,", "\\|", mut)
  
  # Handle unphased mutations (marked by `*`)
  unphased <- grep("\\*", res)
  if (length(unphased) > 0) {
    res[unphased] <- gsub("\\|", "\\/", res[unphased])  # Convert `|` to `/`
    res[unphased] <- gsub("\\*", "", res[unphased])    # Remove `*`
  }
  
  # Handle missing data (marked by `-`)
  missing <- grep("\\-", res)
  if (length(missing) > 0) {
    res[missing] <- NA
  }
  
  # Remove redundant SNPs (e.g., `G|G`)
  rept <- strsplit(res, "\\|")
  if (any(lengths(rept) == 2)) {
    pos <- which(lengths(rept) == 2)
    fix <- unlist(lapply(rept[pos], function(x) x[1] == x[2]))
    if (any(fix)) {
      correction <- unlist(lapply(rept[pos[which(fix)]], function(x) x[1]))
      res[pos[which(fix)]] <- correction
    }
  }
  
  # Append gene locus and SNP identifier
  res[!is.na(res)] <- paste0(gl, ":", res[!is.na(res)])
  
  return(res)
}
