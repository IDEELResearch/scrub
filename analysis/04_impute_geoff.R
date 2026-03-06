# 04_impute_geoff.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-raw/k13_ref_protein_codon_dictionary.csv
# - data-derived/geoff_to_k13_coverage.csv
# - data-geoff/study_name/study_name_STAVE.xlsx (for all study_name within data-geoff)
#
# Outputs:
# - data-geoff/study_name/study_name_STAVE_imputed.xlsx (for all study_name within data-geoff)
#
# Purpose:
# Reads in the _STAVE.xlsx objects, which are already in STAVE-compatible
# format. However, these tables do not currently contain entries for all markers
# of interest, i.e. crt:76, mdr1:86, and WHO candidate and validated k13
# markers. Therefore expand to have entries for all markers - even if these
# entries reflect that only wildtype was observed. This is done by reading in
# the known k13 coverage for each study and subsetting to the positions that are
# theoretically covered by this range.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)

# ------------------------------------------------------------------

# read k13 dictionary and filter to target mutations
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET)) |>
  select(gene = PROTEIN, pos = CODON, AA = REF)

# get k13 coverage for each study
k13_coverage <- read.csv(here("analysis", "data-derived", "geoff_to_k13_coverage.csv")) |>
  select(project_dir, k13_min, k13_max)

# Iterate over each extracted data file
study_folders <- list.dirs(here("analysis", "data-geoff"), recursive = FALSE) |>
  basename()

for (i in seq_along(study_folders)) {
  project_dir <- study_folders[i]
  message(paste("Processing project:", project_dir))
  
  # read in STAVE data
  stave_path <- here("analysis", "data-geoff", sprintf("%s/%s_STAVE.xlsx", project_dir, project_dir))
  df_studies <-  read_excel(stave_path, sheet = "studies")
  df_surveys <-  read_excel(stave_path, sheet = "surveys") |>
    mutate(collection_start = as.Date(collection_start),
           collection_end = as.Date(collection_end),
           collection_day = as.Date(collection_day))
  df_counts <-  read_excel(stave_path, sheet = "counts") |>
    select(study_id, survey_id, variant_string, variant_num, total_num) |>
    mutate(notes = NA) |>
    separate(variant_string, into = c("gene", "pos", "AA"), sep = ":") |>
    mutate(pos = as.numeric(pos),
           from_data = TRUE)
  
  # get k13 coverage range
  project_coverage <- k13_coverage |>
    filter(project_dir == .env$project_dir)
  
  # impute k13 if covered
  df_k13_combined <- NULL
  if (!is.na(project_coverage$k13_min)) {
    
    # filter dictionary to this range
    df_WT <- k13_dictionary |>
      filter(pos >= project_coverage$k13_min,
             pos <= project_coverage$k13_max)
    if (nrow(df_WT) > 0) {
      
      # focus on k13
      df_counts_k13 <- df_counts |>
        filter(gene == "k13")
      
      # expand WT data frame over surveys ready to bind
      df_WT_expanded <- expand_grid(select(df_surveys, study_id, survey_id), df_WT) |>
        mutate(variant_num = NA,
               total_num = NA,
               notes = NA,
               from_data = FALSE)
      
      # bind data with imputed WT
      df_k13_combined <- bind_rows(df_counts_k13, df_WT_expanded) |>
        arrange(survey_id, pos)
      
      # drop imputed WT rows where the WT is already recorded in the data
      df_k13_combined <- df_k13_combined |>
        group_by(survey_id, pos, AA) |>
        mutate(superceded = (from_data == FALSE) & any(from_data)) |>
        ungroup() |>
        filter(!superceded) |>
        select(-superceded)
      
      # for each position, calculate WT numerator and denominator from data where available
      df_k13_combined <- df_k13_combined |>
        group_by(survey_id, pos) |>
        mutate(total_num = ifelse((!from_data) & any(from_data),
                                  min(total_num, na.rm = TRUE),
                                  total_num),
               variant_num = ifelse(!(from_data) & any(from_data),
                                    total_num - sum(variant_num, na.rm = TRUE),
                                    variant_num)) |>
        ungroup()
      
      # for all non-recorded positions, impute as WT with denominator set as
      # minimum over the gene (for this survey)
      df_k13_combined <- df_k13_combined |>
        group_by(survey_id) |>
        mutate(total_num = ifelse(is.na(total_num), min(total_num, na.rm = TRUE), total_num),
               variant_num = ifelse(is.na(variant_num), total_num, variant_num)) |>
        ungroup()
      
      # for k13 we aren't interested in recording zeros (although we will keep these for crt and mdr1)
      df_k13_combined <- df_k13_combined |>
        filter(variant_num > 0)
      
    }
  }
  
  # impute crt if covered
  df_crt_combined <- NULL
  if ("crt" %in% df_counts$gene) {
    
    # define WT data frame
    df_WT <- data.frame(gene = "crt", pos = 76, AA = "K")
    
    # focus on crt
    df_counts_crt <- df_counts |>
      filter(gene == "crt")
    
    # expand WT data frame over surveys ready to bind
    df_WT_expanded <- expand_grid(select(df_surveys, study_id, survey_id), df_WT) |>
      mutate(variant_num = NA,
             total_num = NA,
             notes = NA,
             from_data = FALSE)
    
    # bind data with imputed WT
    df_crt_combined <- bind_rows(df_counts_crt, df_WT_expanded) |>
      arrange(survey_id, pos)
    
    # drop imputed WT rows where the WT is already recorded in the data
    df_crt_combined <- df_crt_combined |>
      group_by(survey_id, pos, AA) |>
      mutate(superceded = (from_data == FALSE) & any(from_data)) |>
      ungroup() |>
      filter(!superceded) |>
      select(-superceded)
    
    # for each position, calculate WT numerator and denominator from data where available
    df_crt_combined <- df_crt_combined |>
      group_by(survey_id, pos) |>
      mutate(total_num = ifelse((!from_data) & any(from_data),
                                min(total_num, na.rm = TRUE),
                                total_num),
             variant_num = ifelse(!(from_data) & any(from_data),
                                  total_num - sum(variant_num, na.rm = TRUE),
                                  variant_num)) |>
      ungroup()
    
    # for all non-recorded positions, impute as WT with denominator set as minimum over the gene
    df_crt_combined <- df_crt_combined |>
      mutate(total_num = ifelse(is.na(total_num), min(total_num, na.rm = TRUE), total_num),
             variant_num = ifelse(is.na(variant_num), total_num, variant_num))
    
  }
  
  # impute crt if covered
  df_mdr1_combined <- NULL
  if ("mdr1" %in% df_counts$gene) {
    
    # define WT data frame
    df_WT <- data.frame(gene = "mdr1", pos = 86, AA = "N")
    
    # focus on mdr1
    df_counts_mdr1 <- df_counts |>
      filter(gene == "mdr1")
    
    # expand WT data frame over surveys ready to bind
    df_WT_expanded <- expand_grid(select(df_surveys, study_id, survey_id), df_WT) |>
      mutate(variant_num = NA,
             total_num = NA,
             notes = NA,
             from_data = FALSE)
    
    # bind data with imputed WT
    df_mdr1_combined <- bind_rows(df_counts_mdr1, df_WT_expanded) |>
      arrange(survey_id, pos)
    
    # drop imputed WT rows where the WT is already recorded in the data
    df_mdr1_combined <- df_mdr1_combined |>
      group_by(survey_id, pos, AA) |>
      mutate(superceded = (from_data == FALSE) & any(from_data)) |>
      ungroup() |>
      filter(!superceded) |>
      select(-superceded)
    
    # for each position, calculate WT numerator and denominator from data where available
    df_mdr1_combined <- df_mdr1_combined |>
      group_by(survey_id, pos) |>
      mutate(total_num = ifelse((!from_data) & any(from_data),
                                min(total_num, na.rm = TRUE),
                                total_num),
             variant_num = ifelse(!(from_data) & any(from_data),
                                  total_num - sum(variant_num, na.rm = TRUE),
                                  variant_num)) |>
      ungroup()
    
    # for all non-recorded positions, impute as WT with denominator set as minimum over the gene
    df_mdr1_combined <- df_mdr1_combined |>
      mutate(total_num = ifelse(is.na(total_num), min(total_num, na.rm = TRUE), total_num),
             variant_num = ifelse(is.na(variant_num), total_num, variant_num))
    
  }
  
  # combine and tidy up
  df_imputed <- bind_rows(df_crt_combined,
                          df_mdr1_combined,
                          df_k13_combined) |>
    arrange(survey_id, gene, pos, AA) |>
    mutate(variant_string = sprintf("%s:%s:%s", gene, pos, AA)) |>
    select(study_id, survey_id, variant_string, variant_num, total_num, notes)
  
  # save back into new xlsx
  stave_path_new <- here("analysis", "data-geoff", sprintf("%s/%s_STAVE_imputed.xlsx", project_dir, project_dir))
  write_xlsx(list(studies = df_studies,
                  surveys = df_surveys,
                  counts = df_imputed),
             path = stave_path_new)
  
}
