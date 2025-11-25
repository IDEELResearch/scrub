#' Standardise codon 86 of mdr1 gene
#' 
#' @details
#' 
#' This function is designed to standardise the mdr1 mutations in codon 86 so that they can be 
#' easily collapsed, combined, cleaned and the complement can be imputed
#' 
#' @name standardise_mdr1_86
#' @export
#' 
#' @param x mdr1 mutation character string
#' @return character string with standardises 86 codon
#' 
#' @examples
#' # test this function
#' standardise_mdr1_86("mdr1:86:N") == "mdr1:86:N"
#' standardise_mdr1_86("mdr1:86:Y/N") == "mdr1:86:N/Y" 
#' standardise_mdr1_86("mdr1:86_184_1246:N_Y_D") == "mdr1:86:N"
#' standardise_mdr1_86("mdr1:86_184:Y/F_N") == "mdr1:86:Y" 
#' 
standardise_mdr1_86 <- function(x) {
  # split into 3 parts: "mdr1", positions, alleles
  parts <- str_split_fixed(x, ":", 3)
  # alleles for *all positions*
  alleles_all <- parts[,3]
  # extract only the allele block for the first position (86)
  first_block <- str_extract(alleles_all, "^[A-Z/]+")
  # keep only A or A/B
  first_allele <- str_extract(first_block, "^[A-Z](?:/[A-Z])?")
  # ----- NEW RULE: drop F if mixed -----
  # N/F → N
  # Y/F → Y
  first_allele <- gsub("^N/F$", "N", first_allele)
  first_allele <- gsub("^Y/F$", "Y", first_allele)
  
  # normalise N/Y ordering
  first_allele <- ifelse(first_allele == "Y/N", "N/Y", first_allele)
  
  paste0("mdr1:86:", first_allele)
}

#' Extract the allele at codon 76 from CRT haplotype strings
#'
#' @description
#' Parses a CRT haplotype string of the form `"crt:<codons>:<alleles>"` and
#' extracts the allele corresponding to codon **76**, regardless of its position
#' within the codon list. The function returns a normalised string of the form
#' `"crt:76:<allele>"`.
#'
#' This function is robust to wide variation in CRT haplotype encodings, such as:
#' - `crt:72_73_74_75_76:CVIET`
#' - `crt:74_75_76_356:M_N_K_I`
#' - `crt:76:T`
#' - `crt:72_76:N_K`
#' - `crt:72_73_74_75_76:C_V_I_E_T`
#'
#' @param x A character vector containing CRT haplotype strings. Each element
#'   must include the gene name (`crt`), a codon position block separated by
#'   underscores, and an allele block with alleles separated by underscores.
#'
#' @return
#' A character vector of the same length as `x`, where each element is of the
#' form `"crt:76:<allele>"`.  
#' If codon 76 is not present in the input string, `NA` is returned for that
#' element.
#'
#' @details
#' The function works by:
#' 1. Splitting the haplotype string into gene, codon block, and allele block.
#' 2. Converting the codon block into a numeric vector (e.g., `"72_73_76"` → `c(72, 73, 76)`).
#' 3. Splitting the allele block into individual alleles (e.g., `"C_V_I_E_T"` → `c("C","V","I","E","T")`).
#' 4. Matching the index of codon **76** in the codon vector.
#' 5. Extracting the allele with the same index from the allele vector.
#'
#' This approach ensures correct matching even when codons are reordered, when
#' non-adjacent codons are present, or when the allele block uses compressed
#' or expanded notation.
#'
#' @examples
#' extract_crt_76("crt:74_75_76_356:M_N_K_I")
#' extract_crt_76("crt:76:T")
#' extract_crt_76("crt:72_76:N_K")
#'
#' @export
extract_crt_76 <- function(x) {
  parts <- stringr::str_split_fixed(x, ":", 3)
  codon_block <- parts[,2]
  allele_block <- parts[,3]
  
  codons <- strsplit(codon_block, "_")
  codons <- lapply(codons, as.integer)
  
  alleles <- strsplit(allele_block, "_")
  
  idx_76 <- lapply(codons, function(v) which(v == 76))
  
  allele_76 <- mapply(function(a, i) {
    if (length(i) == 0) return(NA_character_)
    a[i]
  }, alleles, idx_76)
  
  paste0("crt:76:", allele_76)
}
