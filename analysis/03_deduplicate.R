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

africa_iso3 <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD",
  "COM", "COD", "COG", "DJI", "EGY", "GNQ", "ERI", "SWZ", "ETH", "GAB",
  "GMB", "GHA", "GIN", "GNB", "CIV", "KEN", "LSO", "LBR", "LBY", "MDG",
  "MWI", "MLI", "MRT", "MUS", "MAR", "MOZ", "NAM", "NER", "NGA", "RWA",
  "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TGO", "TUN",
  "UGA", "TZA", "ZMB", "ZWE", "ESH"
)

### TO-DO CECILE: DELETE IF clean_pf7k, clean_who will only be Africa countries in the future
clean_pf7k <- clean_pf7k %>% filter(iso3c %in% africa_iso3)
clean_who <- clean_who %>% filter(iso3c %in% africa_iso3)

# make our full bind across
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

# TO-DO CECILE: improve this function
# Obtain collection year for deduplication
full_bind$collection_year_start <- lubridate::year(full_bind$collection_start)
full_bind$collection_year_end <- lubridate::year(full_bind$collection_end)

# identify studies that may be duplicates
dedup_output = deduplicate(full_bind)
dedup_df = dedup_output$df
summary_same_study = dedup_output$summary_same
summary_diff_study = dedup_output$summary_diff

dim(full_bind)
dim(dedup_df)

dim(summary_same_study)
dim(summary_diff_study)

# save ready to go to stave
saveRDS(dedup_df, here("analysis/data-derived/final_data.rds"))
