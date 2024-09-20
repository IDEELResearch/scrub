# function to clean mutation names from the old format such as those used by WWARN into the long string
# format that we decided upon

clean_mutations <- function(string) {
  gene <- str_split(string,"-")[[1]][1]
  mut <- str_split(string,"-")[[1]][2]
  amino <- gsub("[^a-zA-Z/|]", "", mut) # add symbols indicating mixed infections
  
  if(amino == "WT") {
    gene_mut <- paste0(gene,":",amino)
  } else {
    codon <- parse_number(mut)
    gene_mut <- paste0(gene,":",codon,":",amino)
  }
  return(gene_mut)
}

## test function
# clean_mutations("mdr1-86Y")
# clean_mutations("k13-WT")
# clean_mutations("k13-580Y")
