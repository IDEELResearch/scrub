#' Process a single survey ID to impute missing k13 mutation records
#'
#' This function processes an individual survey ID from the master_table_formatted 
#' dataset to identify and impute missing k13 mutation records. It checks the 
#' coverage range from k13_min:k13_max in study_overview, extracts existing k13 
#' mutations, and adds missing variants based  on the reference mutation key 
#' defining variant_num and total_num from average of total_num of all k13 alleles 
#' in a survey between the specified k13 range from the study overview.
#' 
#' @param each_survey A character string representing a unique survey ID to be processed.
#' @return A tibble containing imputed records for missing k13 mutations, or NULL if no 
#'         imputation is needed or the survey has missing/inconsistent k13 range data.
#' @examples
#' # Example usage (not run):
#' # process_survey("S0001")
#' 
#' @export
process_survey <- function(each_survey) {
  #print(paste("Processing survey_ID:", each_survey))  # debug message
  setTxtProgressBar(pb, which(unique_surveys == each_survey))  # update progress bar
  
  master_table_each_survey <- master_table_formatted %>%
    filter(survey_ID == each_survey)
  
  existing_k13_ranges <- master_table_each_survey %>%
    filter(str_detect(variant_string, "^k13:")) %>%
    pull(variant_string)
  
  existing_codons <- str_extract_all(existing_k13_ranges, "(?<=k13:)([\\d_]+)") %>%
    unlist() %>%
    str_split("_") %>%
    unlist() %>%
    as.numeric() %>%
    unique()
  
  unique_k13_min <- unique(master_table_each_survey$k13_min)
  unique_k13_max <- unique(master_table_each_survey$k13_max)
  
  # Convert empty strings ("") to NA and ensure numeric
  unique_k13_min <- suppressWarnings(as.numeric(ifelse(unique_k13_min == "", NA, unique_k13_min)))
  unique_k13_max <- suppressWarnings(as.numeric(ifelse(unique_k13_max == "", NA, unique_k13_max)))
  
  # Skip surveys with missing or invalid k13_min/k13_max
  if (length(unique_k13_min) != 1 || length(unique_k13_max) != 1 || is.na(unique_k13_min) || is.na(unique_k13_max)) {
    message(sprintf("Skipping survey_ID: %s - Missing or invalid k13_min/k13_max values", each_survey))
    return(NULL)  # Skip this survey
  }
  
  k13_covered_range <- unique_k13_min:unique_k13_max  # create range safely
  
  missing_codons <- setdiff(k13_covered_range, existing_codons)
  
  # Only keep codons that exist in `k13_mutations`
  # Extract numeric codon numbers from k13_mutations$mut
  valid_k13_codons <- as.numeric(gsub("[A-Z]$", "", k13_mutations$mut))  # Remove trailing letter
  missing_codons <- missing_codons[missing_codons %in% valid_k13_codons]  # Keep only valid codons
  
  missing_variants <- mutation_key %>%
    filter(CODON %in% missing_codons) %>%
    mutate(variant_string = paste0("k13:", CODON, ":", REF),
           formatted_mut = paste0(CODON, ":", REF)) %>%
    select(CODON, REF, variant_string, formatted_mut)
  
  mean_total_num <- master_table_each_survey %>%
    filter(tolower(gene) == "k13") %>%
    pull(total_num) %>%
    mean(na.rm = TRUE) %>%
    round()
  
  if (nrow(missing_variants) > 0) {
    template_row <- master_table_each_survey[1, ] %>%
      select(-gene, -mut, -variant_string, -total_num, -variant_num, -survey_ID)
    
    new_records <- map_dfr(1:nrow(missing_variants), function(i) {
      new_row <- template_row
      new_row$gene <- "k13"
      new_row$mut <- missing_variants$formatted_mut[i]  
      new_row$variant_string <- missing_variants$variant_string[i]
      new_row$total_num <- mean_total_num
      new_row$variant_num <- mean_total_num  
      new_row$survey_ID <- each_survey
      new_row$prev <- 1
      return(new_row)
    })
    
    return(new_records)
  } else {
    return(NULL)
  }
}