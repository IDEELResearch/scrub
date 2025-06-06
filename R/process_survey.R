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
  # Update the progress bar (external objects: pb, unique_surveys assumed)
  setTxtProgressBar(pb, which(unique_surveys == each_survey))
  
  # === STEP 1: Filter to records for this survey only ===
  master_table_each_survey <- master_table_formatted %>%
    filter(survey_ID == each_survey)
  
  # === STEP 2: Extract all k13 variant records ===
  k13_variants_all <- master_table_each_survey %>%
    filter(str_detect(variant_string, "^k13:"))
  
  # === STEP 3: Calculate mean total_num from ALL k13 records (valid or not) ===
  # This ensures even coverage-reported ranges/wildcards contribute to the mean
  mean_total_num <- k13_variants_all %>%
    pull(total_num) %>%
    mean(na.rm = TRUE) %>%
    round()
  
  # === STEP 4: Separate valid from invalid k13 records ===
  # Valid: individual codons; Invalid: entries with dash (-) or asterisk (*)
  k13_valid <- k13_variants_all %>%
    filter(!str_detect(variant_string, "[-*]"))
  
  k13_removed <- k13_variants_all %>%
    filter(str_detect(variant_string, "[-*]"))
  
  # === STEP 5: Remove invalid k13 records from the survey dataset ===
  master_table_each_survey <- anti_join(master_table_each_survey, k13_removed, by = colnames(k13_removed))
  
  # === STEP 6: Identify codon numbers already present in valid k13 entries ===
  # Extract numeric codons from variant strings like "k13:622:I"
  existing_codons <- k13_valid %>%
    pull(variant_string) %>%
    str_extract("(?<=k13:)[0-9]+") %>%
    as.numeric() %>%
    unique()
  
  # === STEP 7: Clean and validate the k13_min and k13_max values ===
  unique_k13_min <- unique(master_table_each_survey$k13_min)
  unique_k13_max <- unique(master_table_each_survey$k13_max)
  
  # Convert blank strings to NA and cast to numeric
  unique_k13_min <- suppressWarnings(as.numeric(ifelse(unique_k13_min == "", NA, unique_k13_min)))
  unique_k13_max <- suppressWarnings(as.numeric(ifelse(unique_k13_max == "", NA, unique_k13_max)))
  
  # Skip if min/max are invalid or inconsistent
  if (length(unique_k13_min) != 1 || length(unique_k13_max) != 1 ||
      is.na(unique_k13_min) || is.na(unique_k13_max)) {
    message(sprintf("Skipping survey_ID: %s - Missing or invalid k13_min/k13_max values", each_survey))
    return(NULL)
  }
  
  # === STEP 8: Build full covered codon range based on min/max ===
  k13_covered_range <- unique_k13_min:unique_k13_max
  
  # === STEP 9: Identify codons missing from the valid entries ===
  # Reference valid codons from the mutation key
  valid_k13_codons <- as.numeric(gsub("[A-Z]$", "", k13_mutations$mut))
  
  missing_codons <- setdiff(k13_covered_range, existing_codons)
  missing_codons <- missing_codons[missing_codons %in% valid_k13_codons]
  
  # === STEP 10: Construct imputed variant entries for the missing codons ===
  missing_variants <- mutation_key %>%
    filter(CODON %in% missing_codons) %>%
    mutate(
      variant_string = paste0("k13:", CODON, ":", REF),
      formatted_mut = paste0(CODON, ":", REF)
    ) %>%
    select(CODON, REF, variant_string, formatted_mut)
  
  # === STEP 11: Generate imputed rows if missing variants exist ===
  if (nrow(missing_variants) > 0) {
    # Use the first row as a metadata template, excluding mutation-specific columns
    template_row <- master_table_each_survey[1, ] %>%
      select(-gene, -mut, -variant_string, -total_num, -variant_num, -survey_ID)
    
    # Impute new rows with average total coverage
    new_records <- map_dfr(1:nrow(missing_variants), function(i) {
      new_row <- template_row
      new_row$gene <- "k13"
      new_row$mut <- missing_variants$formatted_mut[i]
      new_row$variant_string <- missing_variants$variant_string[i]
      new_row$total_num <- mean_total_num
      new_row$variant_num <- mean_total_num  # Assume full allele representation
      new_row$survey_ID <- each_survey
      new_row$prev <- 1  # Full prevalence assumed for imputed data
      return(new_row)
    })
    
    return(new_records)
  } else {
    # Nothing to impute
    return(NULL)
  }
}
