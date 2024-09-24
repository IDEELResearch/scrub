#' Clean mutation names from old formats to standardized format
#'
#' This function cleans mutation names, converting from the old format (e.g., those used by WWARN) 
#' into a standardized long string format. The function standardizes separators, extracts gene and 
#' mutation information, and properly formats codons and amino acid mutations.
#'
#' @param string A string representing a mutation name in the old format (e.g., "mdr1-86Y", "k13-WT").
#' 
#' @return A cleaned mutation name string in the format "gene:codon:amino" for mutations or 
#' "gene:WT" for wild-type (WT) or "gene:CNV" for copy number variants.
#' 
#' @export
#' 
#' @examples
#' clean_mutations("mdr1-86Y")  # returns "mdr1:86:Y"
#' clean_mutations("k13-WT")    # returns "k13:WT"
#' clean_mutations("k13-580Y")  # returns "k13:580:Y"
#'
clean_mutations <- function(string) {
  # Replace underscores with hyphens
  string <- stringr::str_replace_all(string, "_", "-")
  
  # Split the string into gene and mutation
  gene <- stringr::str_split(string, "-")[[1]][1]
  mut <- stringr::str_split(string, "-")[[1]][2]
  
  # Extract only letters and specific symbols from the mutation part (e.g., mixed infections)
  amino <- gsub("[^a-zA-Z/|]", "", mut)
  
  # Process based on the type of mutation
  if (amino == "WT") {
    gene_mut <- paste0(gene, ":", amino)
  } else if (amino == "CNV") {
    gene_mut <- paste0(gene, ":", amino)
  } else {
    # Extract the codon as a number and combine with the amino acid mutation
    codon <- readr::parse_number(mut)
    gene_mut <- paste0(gene, ":", codon, ":", amino)
  }
  
  # Return the cleaned mutation name
  return(gene_mut)
}
