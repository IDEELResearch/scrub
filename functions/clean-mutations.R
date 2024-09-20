# function to clean mutation names from the old format such as those used by WWARN into the long string
# format that we decided upon
# updated to handle the formats in the pf7k database and pull out the things we are interested in

clean_mutations <- function(string) {
  string <- str_replace_all(string, "_", "-")
  gene <- str_split(string,"-")[[1]][1]
  mut <- str_split(string,"-")[[1]][2]
  amino <- gsub("[^a-zA-Z/|]", "", mut) # add symbols indicating mixed infections
  
  if(amino == "WT") {
    gene_mut <- paste0(gene,":",amino)
  } else if(amino == "CNV"){
    gene_mut <- paste0(gene,":",amino)
  } 
  else {
    codon <- parse_number(mut)
    gene_mut <- paste0(gene,":",codon,":",amino)
  }
  return(gene_mut)
}

# FOR UNDERSTANDING
# Upper case = homozygous mutations
# lower case = heterozygous
# lower case without second mutation noted is heterozygous at SNP level but same AA
# , dictates separate and distinct haplotypes (different clones)
# / dictates additional NS changes in the same clone
# * in the context of mutations, these indicate that the sample could not be phased

# specific function for dealing with the features of the pf7k k13 encoding
clean_pf7k_k13 <- function(string) {
  string <- str_replace_all(string, "_", "-")
  gene <- str_split(string,"-")[[1]][1]
  mut <- str_split(string,"-")[[1]][2]
  amino <- gsub("[^a-zA-Z/|]", "", mut) # add symbols indicating mixed infections
  
  if(amino == "WT") {
    gene_mut <- paste0(gene,":",amino)
  } else if(amino == "CNV"){
    gene_mut <- paste0(gene,":",amino)
  } 
  else {
    codon <- parse_number(mut)
    gene_mut <- paste0(gene,":",codon,":",amino)
  }
  return(gene_mut)
}

## test function
# clean_mutations("mdr1-86Y")
# clean_mutations("k13-WT")
# clean_mutations("k13-580Y")
