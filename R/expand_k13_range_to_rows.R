#' Expand K13 Mutation Range Into Individual Codon Rows
#'
#' This function takes a `gene_mutation` string in the format "k13:start-end:*" 
#' and returns multiple rows, one for each codon within a specified mutation list.
#'
#' @param gene_mutation A string in format "k13:start-end:*"
#' @param row_data The full data row corresponding to this gene_mutation (as a named list or data.frame row)
#' @param mutation_key Lookup table with columns: PROTEIN, CODON, REF
#' @param mutation_positions Vector of candidate/validated codon positions
#' 
#' @return A data.frame with expanded rows, or NULL if no matching codons
#' 
expand_k13_range_to_rows <- function(gene_mutation, row_data, mutation_key, mutation_positions) {
  gene_mutation <- tolower(gene_mutation)
  
  # Ensure valid range string
  if (!grepl("^k13:[0-9]+-[0-9]+:\\*$", gene_mutation)) {
    return(NULL)
  }
  
  parts <- unlist(strsplit(gene_mutation, "[:-]"))
  start_pos <- as.numeric(parts[2])
  end_pos <- as.numeric(parts[3])
  
  # Filter matching codons
  matching_codons <- mutation_key %>%
    dplyr::filter(
      PROTEIN == "k13",
      CODON >= start_pos,
      CODON <= end_pos,
      CODON %in% mutation_positions
    ) %>%
    dplyr::arrange(CODON)
  
  if (nrow(matching_codons) == 0) {
    return(NULL)
  }
  
  # Create expanded rows
  expanded_rows <- lapply(seq_len(nrow(matching_codons)), function(i) {
    new_row <- as.list(row_data)
    codon <- matching_codons$CODON[i]
    ref <- matching_codons$REF[i]
    
    new_row$gene_mutation <- paste0("K13:", codon, ":", ref)
    new_row$mutant_num <- new_row$total_num
    
    as.data.frame(new_row, stringsAsFactors = FALSE)
  })
  
  return(dplyr::bind_rows(expanded_rows))
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
