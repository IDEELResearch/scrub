library(tidyverse)
library(here)

# source functions
devtools::load_all()

# Read each file if it exists
safe_read <- function(path) {
  if (file.exists(path)) {
    clean <- readRDS(path) %>% 
      as.data.frame() 
    if ("study_type" %in% names(clean)) {
      clean <- clean %>% 
        dplyr::filter(study_type != "unpublished") %>% 
        dplyr::filter(study_type != "data_incomplete")
        
    }
  } else {
    clean <- data.frame()
  }
  return(clean)
}

# Read in each cleaned file
clean_geoff <- safe_read(here("analysis", "data-derived", "geoff_clean.rds"))
clean_wwarn <- safe_read(here("analysis", "data-derived", "wwarn_clean.rds"))
clean_pf7k <- safe_read(here("analysis", "data-derived", "pf7k_clean.rds"))
clean_who <- safe_read(here("analysis", "data-derived", "who_clean.rds"))

# Define African iso3 countries
africa_iso3 <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD",
  "COM", "COD", "COG", "DJI", "EGY", "GNQ", "ERI", "SWZ", "ETH", "GAB",
  "GMB", "GHA", "GIN", "GNB", "CIV", "KEN", "LSO", "LBR", "LBY", "MDG",
  "MWI", "MLI", "MRT", "MUS", "MAR", "MOZ", "NAM", "NER", "NGA", "RWA",
  "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TGO", "TUN",
  "UGA", "TZA", "ZMB", "ZWE", "ESH"
)

# Combine all cleaned df into on dataframe
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% filter(iso3c %in% africa_iso3) %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_wwarn %>% filter(iso3c %in% africa_iso3) %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% filter(iso3c %in% africa_iso3) %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% filter(iso3c %in% africa_iso3) %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

# Obtain collection year for deduplication
full_bind$collection_year_start <- lubridate::year(full_bind$collection_start)
full_bind$collection_year_end <- lubridate::year(full_bind$collection_end)

# ---- Remove duplicates based on pmid ----
duplicate_pmid <- deduplicate_pmid(full_bind)
duplicate_pmid_list <- duplicate_pmid$tagged_list
duplicate_pmid_df <- duplicate_pmid$deduplicated_df

# Remove data points that have duplicate IDs
filtered_pmid_df <- duplicate_pmid_df %>%
  filter(is.na(keep_row) | keep_row == "keep")

# Identify data points that should be removed
removed_pmid_df <- duplicate_pmid_df %>%
  filter(keep_row == "remove")

# ---- Remove duplicates based on data ----
# Identify studies that may be duplicates
dedup_output = deduplicate(filtered_pmid_df)
# Final dataframe with added column indicating if a row should be keeped or removed
dedup_df = dedup_output$df 
# Dataframe showing duplicates within the same study_ID
duplicate_same_list = dedup_output$list_same
# Dataframe showing duplicates across studies
duplicate_diff_list = dedup_output$list_diff

# Combine the tagged duplictaed rows with the duplictate pmid rows
dedup_df_overall <- rbind(dedup_output$df, removed_pmid_df)

# save ready to go to stave
saveRDS(dedup_df_overall, here("analysis/data-derived/final_data.rds"))

##### SUMMARY STATS
create_summary <- function(duplicate_list) {
  # Step 1: Tag each list item with its database_signature
  database_signatures <- sapply(duplicate_list, function(df) {
    paste(sort(unique(df$database)), collapse = "|")
  })
  
  # Step 2: Get all unique combinations
  unique_combinations <- unique(database_signatures)
  
  # Step 3: For each combination, bind relevant list entries and count unique values
  combo_summary <- do.call(rbind, lapply(unique_combinations, function(combination) {
    # Get indices of list entries with this combination
    matching_indices <- which(database_signatures == combination)
    
    # *** FIX IS APPLIED HERE ***
    # 1. Select the data frames to be combined
    dfs_to_combine <- duplicate_list[matching_indices]
    
    # 2. Convert lat and lon to numeric in each data frame before binding
    # This ensures consistent types for bind_rows()
    dfs_to_combine_fixed <- lapply(dfs_to_combine, function(df) {
      df %>%
        mutate(
          lat = as.numeric(lat),
          lon = as.numeric(lon)
        )
    })
    
    # Bind rows from all matching entries
    combined_df <- bind_rows(dfs_to_combine_fixed)
    # *** END FIX ***
    
    # Keep only rows marked as "keep"
    kept_df <- combined_df %>% filter(keep_row == "keep")
    remove_df <- combined_df %>% filter(keep_row == "remove")
    
    # Count unique studies and surveys (rest of function remains the same)
    data.frame(
      database_signature = combination,
      total_unique_studies = n_distinct(combined_df$study_ID),
      total_keep_studies = n_distinct(kept_df$study_ID),
      total_remove_studies = n_distinct(remove_df$study_ID),
      total_unique_surveys = n_distinct(combined_df$survey_ID),
      total_keep_surveys = n_distinct(kept_df$survey_ID),
      total_remove_surveys = n_distinct(remove_df$survey_ID),
      total_rows = nrow(combined_df),
      total_kept_rows = nrow(kept_df),
      total_remove_rows = nrow(remove_df), 
      n_list_entries = length(matching_indices),
      stringsAsFactors = FALSE
    )
  }))
  
  # Optional: arrange
  combo_summary <- combo_summary %>%
    arrange(desc(total_unique_studies))
  
  return(combo_summary)
}

# PubMed ID
pmid_summary_df <- create_summary(duplicate_pmid_list)

# Different Studies
diff_study_summary_df <- create_summary(duplicate_diff_list)

# Same studies
same_study_summary_df <- create_summary(duplicate_same_list)

geoff_entries <- duplicate_pmid_list[
  sapply(duplicate_pmid_list, function(df) any(df$database == "GEOFF"))
]

# Obtain overall deduplication numbers
dedup_df_overall <- bind_rows(dedup_output$df %>% filter(keep_row == "remove"), removed_pmid_df)

# Calculate the totals for removed items
total_removed_results <- dedup_df_overall  %>%
  summarise(
    total_removed_datapoints = n(),
    total_removed_surveys = n_distinct(survey_ID),
    total_removed_studies = n_distinct(study_ID)
  )

print(total_removed_results)

# Generate overall summary 
total_dedup_list <- c(duplicate_pmid_list, duplicate_diff_list, duplicate_same_list)
total_dedup_df <- create_summary(total_dedup_list)
