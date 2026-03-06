# 05_merge_geoff.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-geoff/study_name/study_name_STAVE_imputed.xlsx (for all study_name within data-geoff)
#
# Outputs:
# - data-derived/geoff_STAVE.rds
#
# Purpose:
# Reads in the _STAVE_imputed.xlsx objects for each study from the GEOFF
# systematic review, and combines into a single STAVE object. Saves this object
# to file.
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

# ------------------------------------------------------------------

# Iterate over each extracted data file
study_folders <- list.dirs(here("analysis", "data-geoff"), recursive = FALSE) |>
  basename()

# create empty STAVE object
s <- STAVE_object$new()

for (i in seq_along(study_folders)) {
  project_dir <- study_folders[i]
  message(paste("Processing project:", project_dir))
  
  # Read in data from xlsx
  xl_path <- here("analysis", "data-geoff", sprintf("%s/%s_STAVE_imputed.xlsx", project_dir, project_dir))
  studies_table <-  read_excel(xl_path, sheet = "studies")
  surveys_table <-  read_excel(xl_path, sheet = "surveys") |>
    mutate(collection_start = as.Date(collection_start),
           collection_end = as.Date(collection_end),
           collection_day = as.Date(collection_day))
  counts_table <-  read_excel(xl_path, sheet = "counts")
  
  # append to STAVE object
  s$append_data(studies_dataframe = studies_table,
                surveys_dataframe = surveys_table,
                counts_dataframe = counts_table)
}

# combine both MSMT studies (no reason for years to be separated)
MSMT_studies <- s$get_studies() |>
  filter(study_id == "s0008a_msmt_tza_22") |>
  mutate(study_id = "s0008_msmt_tza")

MSMT_surveys <- s$get_surveys() |>
  filter(study_id %in% c("s0008a_msmt_tza_22", "s0008b_msmt_tza_23")) |>
  mutate(study_id = "s0008_msmt_tza")

MSMT_counts <- s$get_counts() |>
  filter(study_id %in% c("s0008a_msmt_tza_22", "s0008b_msmt_tza_23")) |>
  mutate(study_id = "s0008_msmt_tza")

s$drop_study(c("s0008a_msmt_tza_22", "s0008b_msmt_tza_23"))
s$append_data(studies_dataframe = MSMT_studies,
              surveys_dataframe = MSMT_surveys,
              counts_dataframe = MSMT_counts)

s

# write to file
saveRDS(s, file = here("analysis", "data-derived", "geoff_STAVE.rds"))

