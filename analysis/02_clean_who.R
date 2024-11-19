# Load required libraries
library(tidyverse)
library(lubridate)
library(here)
devtools::load_all()
################################################################################
#
# Script: 02_clean_who.R
# Purpose: Perform various data cleaning to allow deduplication and combining
# with other data sources within {stave} object
# 
################################################################################

# read in data
# Load the combined geoff data table created in the first script
who <- readRDS(here("analysis", "data-derived", "who_res.rds"))

# rename the columns so that they make sense
column_names <- get_column_names_for_clean()

# first add in extra information needed within column_names
who_edit <- who %>% 
  ungroup() %>% 
  dplyr::mutate(database = "WHO",
                iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c"),
                study_ID = paste0("who_", gsub("-", "_", ID)),
                survey_ID = paste0("who_", paste0(ID, "_", gsub("-", "_", admin_1), "_", study_start_year))) %>% 
  dplyr::mutate(survey_ID = iconv(gsub(" |\\,|\\'", "", survey_ID), from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(study_ID = iconv(study_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(continent = countrycode::countrycode(sourcevar = iso3c, "iso3c", "continent"),
                spatial_notes = "lat/lon sourced from who meta file",
                time_notes = "Automated midpoint", 
                study_name = source,
                study_type = "other",
                authors = NA # to source some other way/time
                ) %>% 
  dplyr::mutate(study_type = replace(study_type, pmid != "", "peer_reviewed")) %>% 
  dplyr::rename(country_name = admin_0,
                site_name = admin_1, 
                lon = long, 
                variant_num = x, 
                total_num = n, 
                variant_string = mut)

# sort out collecton dates
convstart <- adjust_invalid_date(unique(who_edit$study_start_year), is_start = TRUE)
who_edit$collection_start <- convstart[match(who_edit$study_start_year, unique(who_edit$study_start_year))]
convend <- adjust_invalid_date(unique(who_edit$study_start_year), is_start = FALSE)
who_edit$collection_end <- convend[match(who_edit$study_start_year, unique(who_edit$study_start_year))]

# Compute collection day
who_edit <- who_edit %>% rowwise() %>% 
  mutate(collection_day = median(c(collection_start, collection_end))) %>% 
  ungroup()

# Grab just the columns we need for pairing with WWARN etc
who_clean <- who_edit %>% 
  select(all_of(column_names))

# Save the final merged_df as an RDS file
saveRDS(who_clean, here("analysis", "data-derived", "who_clean.rds"))

