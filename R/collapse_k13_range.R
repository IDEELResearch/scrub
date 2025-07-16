#' Collapse K13 Mutation Range
#'
#' This function takes a `gene_mutation` string in the format "k13:start-end:*" 
#' and collapses the codon positions into a standardized format, returning 
#' concatenated amino acid and codon information. The function filters the 
#' reference amino acids based on the specified start and end positions.
#'
#' @param gene_mutation A string representing the mutation in the format "k13:start-end:*".
#' @param mutation_key Dataframe for mutation codon lookup
#' @return A string in the format "K13:start_pos_end_pos:amino_acids", or the original 
#'         `gene_mutation` if no matching codons are found in the reference table.
#' @examples
#' mutation_key <- data.frame(PROTEIN = "k13", CODON = 440:450, REF = letters[1:11])
#' collapse_k13_range("k13:440-445:*", mutation_key)
#' # Returns: "K13:440_441_442_443_444_445:a_b_c_d_e_f"
#' 
#' @export
collapse_k13_range <- function(gene_mutation, mutation_key) {
  parts <- unlist(strsplit(gene_mutation, "[:-]"))
  start_pos <- as.numeric(parts[2])
  end_pos <- as.numeric(parts[3])
  
  # Filter the reference amino acids for the specified range
  ref_amino_acids <- mutation_key %>%
    dplyr::filter(.data$PROTEIN == "k13", .data$CODON >= start_pos, .data$CODON <= end_pos) %>%
    dplyr::arrange(.data$CODON)
  
  # If no matching codons are found, return the original value and log the issue
  if (nrow(ref_amino_acids) == 0) {
    cat("Warning: No matching codons found for gene_mutation:", gene_mutation, "\n")
    return(gene_mutation)
  }
  
  # Concatenate codon positions and amino acids
  codon_string <- paste(ref_amino_acids$CODON, collapse = "_")
  amino_acid_string <- paste(ref_amino_acids$REF, collapse = "_")
  collapsed_mutation <- paste("K13", codon_string, amino_acid_string, sep = ":")
  
  return(collapsed_mutation)
}

#' Convert K13 Asterisk Mutation
#'
#' This function converts a single-codon mutation string in the format `"k13:codon:*"` 
#' into a standardized format by looking up the reference amino acid from a provided mutation key dataframe.
#' If the specified codon position does not exist in the mutation key, it will return `"K13:codon:NA"`.
#'
#' @param gene_mutation A string representing the mutation in the format `"k13:codon:*"`.
#' @param mutation_key A dataframe containing mutation codon information with columns 
#'        `"PROTEIN"`, `"CODON"`, and `"REF"`.
#' @return A string in the format `"K13:codon:amino_acid"`. If no matching codon is found, returns 
#'         `"K13:codon:NA"`.
#' @examples
#' mutation_key <- data.frame(PROTEIN = "k13", CODON = 440:450, REF = letters[1:11])
#' convert_k13_asterisk("k13:442:*", mutation_key)
#' # Returns: "K13:442:c"
#'
#' @export
convert_k13_asterisk <- function(gene_mutation, mutation_key) {
  parts <- strsplit(gene_mutation, ":")
  codon <- as.numeric(unlist(lapply(parts, "[[", 2)))
  
  # Filter the reference amino acids for the specific SNP
  ref_amino_acids <- mutation_key[mutation_key$PROTEIN == "k13",]$REF[
    match(codon, mutation_key[mutation_key$PROTEIN == "k13",]$CODON)
  ]
  
  return(paste("K13", codon, ref_amino_acids, sep = ":"))
  
}
