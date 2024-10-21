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

# Load the combined all studies geoff data table created in the first script
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
master_table_renamed <- create_renamed_df(master_table, column_mapping)

# Ensure consistent column types for merging
master_table_renamed <- master_table_renamed %>%
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

all_columns <- union(colnames(master_table_renamed), wwarn_res_df_cols)

# Fill any columns still required for wwarn rbind which have no comparable data in geoff with na to allow rbind later on
for (col in setdiff(all_columns, colnames(master_table_renamed))) {
  master_table_renamed[[col]] <- NA
}

# Reorder columns to match the wwarn_res_df_cols structure
master_table_renamed <- master_table_renamed[, all_columns]

# Save the final merged_df as an RDS file
saveRDS(master_table_renamed, here("analysis", "data-derived", "02_clean_geoffs_output_table.rds"))
print("Data saved. Ready for further analysis.")
