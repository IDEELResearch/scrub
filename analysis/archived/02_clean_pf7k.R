# Load required libraries
library(tidyverse)
library(lubridate)
library(here)

################################################################################
#
# Script: 02_clean_pf7k.R
# Purpose: Perform various data cleaning to allow deduplication and combining
# with other data sources within {stave} object
# 
################################################################################

# read in data
# Load the combined geoff data table created in the first script
pf7k <- readRDS(here("analysis", "data-derived", "pf7k_res.rds"))

# rename the columns so that they make sense
column_names <- get_column_names_for_clean()

# first add in extra information needed within column_names
pf7k_edit <- pf7k %>% 
  ungroup() %>% 
  dplyr::mutate(database = "Pf7k",
                prev = variant_num/total_num,
                gene = gsub("(^[[:alnum:]]*)_(.*)", "\\1", name),
                iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c"),
                study_ID = paste0("pf7k_", gsub("-", "_", study_ID)),
                survey_ID = paste0("pf7k_", gsub("-", "_", gsub(" ", "", survey_ID)))) %>% 
  dplyr::mutate(survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(continent = countrycode::countrycode(sourcevar = iso3c, "iso3c", "continent"),
                spatial_notes = "lat/lon sourced from pf7k meta file",
                time_notes = "Automated midpoint") %>% 
  dplyr::rename(country_name = admin_0,
                site_name = admin_1, 
                lon = long)

# sort out collecton dates
convstart <- adjust_invalid_date(unique(pf7k_edit$year), is_start = TRUE)
pf7k_edit$collection_start <- convstart[match(pf7k_edit$year, unique(pf7k_edit$year))]
convend <- adjust_invalid_date(unique(pf7k_edit$year), is_start = FALSE)
pf7k_edit$collection_end <- convend[match(pf7k_edit$year, unique(pf7k_edit$year))]

# Compute collection day
pf7k_edit <- pf7k_edit %>% rowwise() %>% 
  mutate(collection_day = median(c(collection_start, collection_end))) %>% 
  ungroup()

# Grab just the columns we need for pairing with WWARN etc
pf7k_clean <- pf7k_edit %>% 
  select(all_of(column_names))

# Save the final merged_df as an RDS file
saveRDS(pf7k_clean, here("analysis", "data-derived", "pf7k_clean.rds"))

