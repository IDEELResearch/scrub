#' Collapse K13 Mutation Range
#'
#' This function takes a `gene_mutation` string in the format "k13:start-end:*" 
#' and collapses the codon positions into a standardized format, returning 
#' concatenated amino acid and codon information. The function filters the 
#' reference amino acids based on the specified start and end positions.
#'
#' @param gene_mutation A string representing the mutation in the format "k13:start-end:*".
#' @return A string in the format "K13:start_pos_end_pos:amino_acids", or the original 
#'         `gene_mutation` if no matching codons are found in the reference table.
#' @examples
#' mutation_key <- data.frame(PROTEIN = "k13", CODON = 440:450, REF = letters[1:11])
#' collapse_k13_range("k13:440-445:*")
#' # Returns: "K13:440_441_442_443_444_445:a_b_c_d_e_f"
#' 
#' @export
collapse_k13_range <- function(gene_mutation) {
  parts <- unlist(strsplit(gene_mutation, "[:-]"))
  start_pos <- as.numeric(parts[2])
  end_pos <- as.numeric(parts[3])
  
  # Filter the reference amino acids for the specified range
  ref_amino_acids <- mutation_key %>%
    filter(PROTEIN == "k13", CODON >= start_pos, CODON <= end_pos) %>%
    arrange(CODON)
  
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
