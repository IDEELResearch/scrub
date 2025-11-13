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
  
  # Case 6: Mixed phased/unphased (_ and / both present)
  if (grepl("\\_", marker) && grepl("\\/", marker)) {
    spl <- toupper(strsplit(marker, "\\_")[[1]])  # split by '_'
    
    all_loci <- c()
    all_snps <- c()
    
    for (s in spl) {
      # If there is a "/", handle it accordingly
      if (grepl("\\/", s)) {
        # Extract base locus, e.g., from "N197D/N" → 197
        loc <- gsub("(^\\w)(\\d+)(\\w)(\\/\\w)?", "\\2", s)
        snp <- gsub("(^\\w)(\\d+)(\\w\\/\\w)", "\\3", s)
      } else {
        loc <- gsub("(^\\w)(\\d+)(\\w)", "\\2", s)
        snp <- gsub("(^\\w)(\\d+)(\\w)", "\\3", s)
      }
      all_loci <- c(all_loci, loc)
      all_snps <- c(all_snps, snp)
    }
    
    res <- paste0("k13:",
                  paste0(sort(as.numeric(all_loci)), collapse = "_"),
                  ":",
                  paste0(all_snps[order(as.numeric(all_loci))], collapse = "_"))
    return(res)
  }
  
  stop("Unexpected mutation")
  
}

#' Convert WWARN crt SNP formats for STAVE
#'
#' @details
#'
#' This function processes crt SNP data into a format compatible with the `STAVE` tool.
#'
#' NB FOR STAVE: "|" reflects phased  
#'
#' NB FOR STAVE: "/" reflects unphased
#'
#' ## Input format conventions:
#' - `WT`: Denotes wild type.
#' - `76T`: Single mutation (codon + amino acid).
#' - `76K/T`: Unphased call (e.g. mixed infection at the same locus).
#' - `-` or `NA`: Treated as missing.
#'
#' ## Output format:
#' - `crt:WT`: Wild type.
#' - `crt:76:T`: Single mutation.
#' - `crt:76:K/T`: Mixed/unphased mutation.
#'
#' Convert WWARN crt SNP formats for STAVE (vectorized)
#'
#' @param marker Character vector. SNP calls (e.g. "76T", "76K/T", "WT", NA).
#' @return Character vector in STAVE-compatible format.
#' @export
#' @examples
#' wwarn_format_crt_for_stave("WT")        # crt:WT
#' wwarn_format_crt_for_stave("76T")       # crt:76:T
#' wwarn_format_crt_for_stave("76K/T")     # crt:76:K/T
#' wwarn_format_crt_for_stave("-")         # NA
#' wwarn_format_crt_for_stave(NA)          # NA
wwarn_format_crt_for_stave <- function(marker) {
  vapply(marker, function(m) {
    if (is.na(m) || m == "-") {
      return(NA_character_)
    }
    
    m <- toupper(m)
    
    if (m == "WT") {
      return("crt:WT")
    }
    
    if (grepl("^\\d+[A-Z]$", m)) {
      loc <- sub("^(\\d+)[A-Z]$", "\\1", m)
      aa <- sub("^\\d+([A-Z])$", "\\1", m)
      return(paste0("crt:", loc, ":", aa))
    }
    
    if (grepl("^\\d+[A-Z]/[A-Z]$", m)) {
      loc <- sub("^(\\d+)[A-Z]/[A-Z]$", "\\1", m)
      aa <- sub("^\\d+([A-Z]/[A-Z])$", "\\1", m)
      return(paste0("crt:", loc, ":", aa))
    }
    
    stop(paste("Unexpected mutation format for crt:", m))
  }, FUN.VALUE = character(1))
}



#' Convert WWARN mdr1 SNP formats for STAVE
#'
#' @details
#'
#' This function processes mdr1 SNP data into a format compatible with the `STAVE` tool.
#'
#' NB FOR STAVE: "|" reflects phased  
#'
#' NB FOR STAVE: "/" reflects unphased
#'
#' ## Input format conventions:
#' - `WT`: Denotes wild type.
#' - `86Y`: Single mutation (codon + amino acid).
#' - `86N/Y`: Unphased call (e.g. mixed infection at the same locus).
#' - `-` or `NA`: Treated as missing.
#'
#' ## Output format:
#' - `mdr1:WT`: Wild type.
#' - `mdr1:86:Y`: Single mutation.
#' - `mdr1:86:N/Y`: Mixed/unphased mutation.
#'
##' Convert WWARN mdr1 SNP formats for STAVE (vectorized)
#'
#' @param marker Character vector. SNP calls (e.g. "86Y", "86N/Y", "WT", NA).
#' @return Character vector in STAVE-compatible format.
#' @export
#' 
#' @examples
#' wwarn_format_mdr1_for_stave("WT")        # mdr1:WT
#' wwarn_format_mdr1_for_stave("86Y")       # mdr1:86:Y
#' wwarn_format_mdr1_for_stave("86N/Y")     # mdr1:86:N/Y
#' wwarn_format_mdr1_for_stave("-")         # NA
#' wwarn_format_mdr1_for_stave(NA)          # NA
wwarn_format_mdr1_for_stave <- function(marker) {
  vapply(marker, function(m) {
    if (is.na(m) || m == "-") {
      return(NA_character_)
    }
    
    m <- toupper(m)
    
    if (m == "WT") {
      return("mdr1:WT")
    }
    
    if (grepl("^\\d+[A-Z]$", m)) {
      loc <- sub("^(\\d+)[A-Z]$", "\\1", m)
      aa <- sub("^\\d+([A-Z])$", "\\1", m)
      return(paste0("mdr1:", loc, ":", aa))
    }
    
    if (grepl("^\\d+[A-Z]/[A-Z]$", m)) {
      loc <- sub("^(\\d+)[A-Z]/[A-Z]$", "\\1", m)
      aa <- sub("^\\d+([A-Z]/[A-Z])$", "\\1", m)
      return(paste0("mdr1:", loc, ":", aa))
    }
    
    stop(paste("Unexpected mutation format for mdr1:", m))
  }, FUN.VALUE = character(1))
}
