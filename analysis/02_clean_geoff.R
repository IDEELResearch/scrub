# Load required libraries
library(dplyr)
library(lubridate)
library(here)

# Load the combined geoff data table created in the first script
master_table <- readRDS(here("analysis", "data-derived", "combined_geoff_data_df.rds"))

# Filter the combined geoff data table for untreated data only
master_table <- master_table %>% filter(!substudy %in% c("treatedextracted", "treatedcalculated"))

# Expand gene mutation ranges for reference range syntax)
rows_to_transform <- master_table %>%
  filter(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(gene_mutation)))

indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table$gene_mutation)))

master_table$gene_mutation[indices_to_transform] <- sapply(
  master_table$gene_mutation[indices_to_transform],
  collapse_k13_range
)

# Adjust and add dates and survey IDs
master_table <- master_table %>%
  mutate(
    collection_start = adjust_invalid_date(date_start, is_start = TRUE),
    collection_end = adjust_invalid_date(date_end, is_start = FALSE),
    collection_day = case_when(
      !is.na(collection_start) & !is.na(collection_end) ~ as.Date((as.numeric(collection_start) + as.numeric(collection_end)) / 2, origin = "1970-01-01"),
      !is.na(collection_start) & is.na(collection_end) ~ collection_start,
      is.na(collection_start) & !is.na(collection_end) ~ collection_end,
      TRUE ~ NA_Date_
    ),
    collection_start = as.character(collection_start),
    collection_end = as.character(collection_end),
    collection_day = as.character(collection_day),
    survey_ID = paste0(study_uid, "_", first_author_surname, "_", site_name, "_", publication_year, "_", sample(1:1000000, size = nrow(master_table), replace = FALSE)),
    survey_ID = gsub("[^a-zA-Z0-9_]", "", survey_ID),
    survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT")
  )

# Extract gene and mutation identifiers
master_table <- master_table %>%
  mutate(
    gene = sub(":.*", "", gene_mutation),
    mut = sapply(strsplit(gene_mutation, ":"), function(x) paste(tail(x, 2), collapse = ":"))
  )

# Define a column mapping for geoff so that it matches columns in wwarn for merging later by Gina
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

# Create combined dataframe based on column mapping
create_combined_df <- function(df, mapping, default_database = "GEOFF") {
  result_df <- data.frame(matrix(ncol = length(mapping), nrow = nrow(df)))
  colnames(result_df) <- names(mapping)
  
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

# Apply transformation
master_table_combined <- create_combined_df(master_table, column_mapping)

# Ensure consistent column types for binding
master_table_combined <- master_table_combined %>%
  mutate(
    study_start_year = as.character(study_start_year),
    study_end_year = as.character(study_end_year),
    x = as.double(x),
    n = as.double(n)
  )

# Add 'prev' column to master_table_combined by dividing x by n
master_table_combined <- master_table_combined %>%
  mutate(prev = if_else(n != 0, x / n, NA_real_))

# Align columns in both dataframes for row binding
all_columns <- union(colnames(master_table_combined), colnames(wwarn_res_df))

# Add missing columns as NA to each dataframe
for (col in setdiff(all_columns, colnames(master_table_combined))) {
  master_table_combined[[col]] <- NA
}

for (col in setdiff(all_columns, colnames(wwarn_res_df))) {
  wwarn_res_df[[col]] <- NA
}

# Reorder columns and bind rows
master_table_combined <- master_table_combined[, all_columns]

# Below line is Just for testing column name and col type compatibility after 
# cleaning geoff to match OJ/Ginas warn cleaning:
# wwarn_res_df <- wwarn_res_df[, all_columns]
# merged_df <- bind_rows(master_table_combined, wwarn_res_df) # Merge geoff combined and wwarn into a single dataframe
# saveRDS(merged_df, here("analysis", "data-derived", "merged_df.rds"))

# Save the final merged_df as an RDS file
saveRDS(master_table_combined, here("analysis", "data-derived", "master_table_combined.rds"))
print("Data saved. Ready for further analysis.")
