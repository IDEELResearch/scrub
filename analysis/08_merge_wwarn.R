# 08_merge_wwarn_k13.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-derived/wwarn_k13_clean.rds
# - data-derived/wwarn_pd_clean.rds
#
# Outputs:
# - data-derived/WWARN_STAVE.rds
#
# Purpose:
# Reads in WWARN data on k13 and partner drug. Merges, wrangles into STAVE
# format, and saves STAVE object to file.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)
#remotes::install_github("mrc-ide/STAVE@v2.0.2")
library(STAVE)

# ------------------------------------------------------------------------
# read in and combine data

wwarn_k13 <- readRDS(here("analysis", "data-derived", "wwarn_k13_clean.rds"))
wwarn_pd <- readRDS(here("analysis", "data-derived", "wwarn_pd_clean.rds"))

wwarn_combined <- wwarn_k13 |>
  mutate(study_start = year,
         study_end = year) |>
  select(-year) |>
  bind_rows(wwarn_pd) |>
  filter(tested != 0)

# ------------------------------------------------------------------------
# make STAVE objects

# study-level data frame
df_studies <- wwarn_combined |>
  group_by(study_id) |>
  summarise(study_label = study_name[1],
            contributors = authors[1],
            reference = url_clean[1],
            reference_year = publication_year[1],
            PMID = PMID[1],
            .groups = "drop") |>
  ungroup() |>
  mutate(description = NA,
         access_level = "public") |>
  select(study_id, study_label, description, access_level, contributors, reference, reference_year, PMID)

# survey-level data frame
df_surveys <- wwarn_combined |>
  group_by(study_id, survey_id, country, site, study_start, study_end) |>
  summarise(latitude = lat[1],
            longitude = lon[1],
            .groups = "drop") |>
  mutate(location_method = "WWARN coordinates",
         location_notes = NA,
         collection_start = as.Date(sprintf("%s-01-01", study_start)),
         collection_end = as.Date(sprintf("%s-12-31", study_end)),
         collection_day = collection_start + (collection_end - collection_start) / 2,
         time_method = "Midpoint of WWARN recorded year",
         time_notes = NA) |>
  rename(country_name = country,
         site_name = site) |>
  select(study_id, survey_id, country_name, site_name, latitude, longitude, location_method, location_notes,
         collection_start, collection_end, collection_day, time_method, time_notes)

# counts-level data frame
df_counts <- wwarn_combined |>
  select(study_id, survey_id, variant_string, variant_num = present, total_num = tested, notes)

# make STAVE object and append data
s <- STAVE::STAVE_object$new()
s$append_data(studies_dataframe = df_studies,
              surveys_dataframe = df_surveys,
              counts_dataframe = df_counts)
s

# write to file
saveRDS(s, file = here("analysis", "data-derived", "WWARN_STAVE.rds"))
