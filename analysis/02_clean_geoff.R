# Load required libraries
library(dplyr)
library(lubridate)
library(here)

################################################################################
#
# Script: 02_clean_geoffs.R
# Purpose: Perform various data cleaning and reformatting functions to 
# facilitate rbind() with wwarn and other cleaned data before deduplication
# and stave append
# 
################################################################################

# Define functions

# Collapse k13 range into a standard format
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

# Adjust invalid dates (e.g., 2019-02-31) with year-only, year-month, and full date support
adjust_invalid_date <- function(date_str, is_start = TRUE) {
  date_fixed <- suppressWarnings(
    case_when(
      grepl("^[0-9]{4}$", date_str) ~ {
        if (is_start) ymd(paste0(date_str, "-01-01")) else ymd(paste0(date_str, "-12-31"))
      },
      grepl("^[0-9]{4}-[0-9]{2}$", date_str) ~ {
        if (is_start) ymd(paste0(date_str, "-01")) else ceiling_date(ymd(paste0(date_str, "-01")), "month") - days(1)
      },
      grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_str) ~ ymd(date_str),
      TRUE ~ as.Date(NA)
    )
  )
  return(date_fixed)
}

# Function to check for gene name typos and print warnings
check_gene_typos <- function(gene_column) {
  # Define probable typos for k13, crt, and mdr1
  probable_typos <- list(
    "k13" = c("kelch13", "kelch 13", "kelch_13", "kelch", "kletch13", "klech 13"),
    "crt" = c("ctr"),
    "mdr1" = c("mrd1", "mdr")
  )
  
  for (correct_gene in names(probable_typos)) {
    typos_found <- gene_column[gene_column %in% probable_typos[[correct_gene]]]
    if (length(typos_found) > 0) {
      message(paste("Probable typo found for", correct_gene, ": did you mean", correct_gene, "?"))
    }
  }
}

# Correcting known bad strings in the substudy column and counting instances
correct_substudy_entries <- function(substudy_column) {
  # Correct entries and count the corrections
  corrections <- list(
    day0extracted = sum(grepl("^day 0$", substudy_column)),  # Count "day 0"
    untreated_extracted = sum(grepl("^untreated_extracted$", substudy_column)),  # Count "untreated_extracted"
    untreadextracted = sum(grepl("^untreadextracted$", substudy_column))  # Count "untreadextracted"
  )
  
  # Perform the corrections
  substudy_column <- gsub("^day 0$", "day0extracted", substudy_column)
  substudy_column <- gsub("^untreated_extracted$", "untreatedextracted", substudy_column)
  substudy_column <- gsub("^untreadextracted$", "untreatedextracted", substudy_column)
  
  # Print the corrections made and their counts
  for (correction in names(corrections)) {
    if (corrections[[correction]] > 0) {
      message(paste("Corrected", corrections[[correction]], "instances of", correction))
    }
  }
  
  return(substudy_column)
}

# Check if all entries are valid for substudy
check_substudy_entries <- function(substudy_column) {
  # Correct bad strings first and track corrections
  substudy_column <- correct_substudy_entries(substudy_column)
  
  # Define allowed categories
  allowed_categories <- c(
    "untreatedextracted", "untreatedcalculated", "day0extracted", "day0calculated",
    "day3calculated", "day24calculated", "treatedextracted", "treatedcalculated"
  )
  
  # Regular expression for "day X extracted" or "day X calculated"
  day_extracted_calculated_pattern <- "^day[0-9]+(extracted|calculated)$"
  
  # Find any entries not matching the allowed categories or the day X extracted/calculated pattern
  invalid_entries <- substudy_column[!substudy_column %in% allowed_categories & 
                                       !grepl(day_extracted_calculated_pattern, substudy_column)]
  
  # If there are invalid entries, print them
  if (length(invalid_entries) > 0) {
    message("Invalid entries found in 'substudy':")
    print(invalid_entries)
  } else {
    message("All 'substudy' entries are valid after corrections.")
  }
  
  return(substudy_column)
}

# Function to create combined dataframe with all columns from the original master_table
create_combined_df <- function(df, mapping, default_database = "GEOFF") {
  result_df <- df  # Start with all columns from the original dataframe
  for (new_col in names(mapping)) {
    if (!is.na(mapping[[new_col]]) && mapping[[new_col]] %in% colnames(df)) {
      result_df[[new_col]] <- df[[mapping[[new_col]]]]
    } else if (new_col == "database") {
      result_df[[new_col]] <- default_database
    } else {
      result_df[[new_col]] <- NA
    }
  }
  return(result_df)
}

# Load data and perform operations
# Load the combined geoff data table created in the first script
master_table <- readRDS(here("analysis", "data-derived", "01_read_geoffs_output_table.rds"))

# Load mutation key for k13 reference ranges
mutation_key_path <- here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")
mutation_key <- read.csv(mutation_key_path)

# Filter the combined geoff data table for untreated data only
master_table <- master_table %>% filter(!substudy %in% c("treatedextracted", "treatedcalculated"))

# Expand gene mutation ranges for reference range syntax
indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table$gene_mutation)))
master_table$gene_mutation[indices_to_transform] <- sapply(
  master_table$gene_mutation[indices_to_transform],
  collapse_k13_range
)

# Adjust and add dates and survey IDs
master_table$collection_start <- sapply(master_table$date_start, adjust_invalid_date, is_start = TRUE)
master_table$collection_end <- sapply(master_table$date_end, adjust_invalid_date, is_start = FALSE)

# Compute collection day
master_table$collection_day <- sapply(1:nrow(master_table), function(i) {
  if (!is.na(master_table$collection_start[i]) & !is.na(master_table$collection_end[i])) {
    return(as.Date(mean(c(as.numeric(master_table$collection_start[i]), as.numeric(master_table$collection_end[i]))), origin = "1970-01-01"))
  } else if (!is.na(master_table$collection_start[i])) {
    return(master_table$collection_start[i])
  } else if (!is.na(master_table$collection_end[i])) {
    return(master_table$collection_end[i])
  } else {
    return(NA)
  }
})

# Convert numeric back to Date class after transformations
master_table$collection_start <- as.Date(master_table$collection_start, origin = "1970-01-01")
master_table$collection_end <- as.Date(master_table$collection_end, origin = "1970-01-01")
master_table$collection_day <- as.Date(master_table$collection_day, origin = "1970-01-01")

# Add survey_ID and other fields
master_table <- master_table %>%
  mutate(
    survey_ID = paste0(study_uid, "_", first_author_surname, "_", site_name, "_", publication_year, "_", sample(1:1000000, size = nrow(master_table), replace = FALSE)),
    survey_ID = gsub("[^a-zA-Z0-9_]", "", survey_ID),
    survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT"),
    gene = sub(":.*", "", gene_mutation),
    mut = sapply(strsplit(gene_mutation, ":"), function(x) paste(tail(x, 2), collapse = ":"))
  )

# Check for gene name typos and print warnings if any are found
check_gene_typos(master_table$gene)

# Apply substudy corrections and validation
master_table$substudy <- check_substudy_entries(master_table$substudy)

# Define a column mapping for the WWARN merge
column_mapping <- list(
  admin_0 = "country",
  admin_1 = NA,
  site = "site_name",
  lat = "lat_n",
  long = "lon_e",
  year = "collection_day",
  study_start_year = "date_start",
  study_end_year = "date_end",
  x = "mutant_num",
  n = "total_num",
  prev = NA,
  gene = "gene",
  mut = "mut",
  gene_mut = "gene_mutation",
  annotation = NA,
  database = "database",
  url = "study_url",
  source = "publication_status",
  iso3c = "iso3c",
  pmid = "pmid"
)

# Create combined dataframe based on column mapping, keeping all original columns
master_table_combined <- create_combined_df(master_table, column_mapping)

# Ensure consistent column types for merging
master_table_combined <- master_table_combined %>%
  mutate(
    study_start_year = as.character(study_start_year),
    study_end_year = as.character(study_end_year),
    x = as.double(x),
    n = as.double(n),
    prev = if_else(n != 0, x / n, NA_real_)
  )

# Add missing columns as NA to each dataframe
wwarn_res_df_cols <- c("iso3c", "admin_0", "admin_1", "site", "lat", "long", "year", 
                       "study_start_year", "study_end_year", "x", "n", "prev", 
                       "gene", "mut", "gene_mut", "annotation", "database", 
                       "pmid", "url", "source")

all_columns <- union(colnames(master_table_combined), wwarn_res_df_cols)

# Fill any columns still required for wwarn rbind which have no comparable data in geoff with na to allow rbind later on
for (col in setdiff(all_columns, colnames(master_table_combined))) {
  master_table_combined[[col]] <- NA
}

# Reorder columns to match the wwarn_res_df_cols structure
master_table_combined <- master_table_combined[, all_columns]

# Save the final merged_df as an RDS file
saveRDS(master_table_combined, here("analysis", "data-derived", "02_clean_geoffs_output_table.rds"))
print("Data saved. Ready for further analysis.")
