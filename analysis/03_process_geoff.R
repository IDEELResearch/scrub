# 03_process_geoff.R
#
# Author: Bob Verity
# Date: 2025-12-14
#
# Inputs:
# - data-geoff/study_name/study_name.xlsx (for all study_name within data-geoff)
# - data-raw/k13_ref_protein_codon_dictionary.csv
#
# Outputs:
# - data-geoff/study_name/study_name_STAVE.xlsx (for all study_name within data-geoff)
#
# Purpose:
# Reads in the basic xlsx data, which has been extracted and cleaned from the
# raw .tsv.gz files. Organizes survey- and count-level data into STAVE format
# and saves these back into the _STAVE.xlsx object, overwriting the existing
# corresponding sheets. Note that the "studies" sheet of the _STAVE.xlsx object
# is not overwritten.
#
# ------------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)
library(writexl)

# Iterate over each extracted data file
study_folders <- list.dirs(here("analysis", "data-geoff"), recursive = FALSE) |>
  basename()

# skip over some folders, where STAVE format was extracted directly
study_folders <- setdiff(study_folders, c("s231_ali_2022_v01",
                                          "s232_wamae_2022_v01",
                                          "s233_issa_2022_v01"))

# read k13 dictionary and filter to target mutations
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET))

for (i in seq_along(study_folders)) {
  project_dir <- study_folders[i]
  message(paste("Processing project:", project_dir))
  
  # make study ID from project directory
  study_id <- gsub("_v01", "", project_dir)
  study_id <- gsub("-", "_", study_id)
  study_id <- tolower(study_id)
  study_id_root <- strsplit(study_id, "_")[[1]][1]
  
  # read in extracted data
  old_xl_path <- here("analysis", "data-geoff", project_dir, sprintf("%s.xlsx", project_dir))
  old_site <- readxl::read_excel(old_xl_path, sheet = "site")
  old_prev <- readxl::read_excel(old_xl_path, sheet = "prevalence")
  
  # drop NA dates
  old_prev <- old_prev |>
    filter(date_start != "NA")
  
  # merge tables and populate where possible
  df_all <- left_join(old_prev, old_site, by = join_by(site_uid)) |>
    mutate(study_id = study_id,
           survey_id = NA,
           country_name = country,
           site_name = site_name,
           latitude = lat_n,
           longitude = lon_e,
           location_method = NA,
           location_notes = NA,
           collection_start = date_start,
           collection_end = date_end,
           collection_day = NA,
           time_method = "imputed as midpoint of collection range",
           time_notes = NA,
           variant_string = gene_mutation,
           variant_num = mutant_num,
           total_num = total_num) |>
    select(study_id, survey_id,	country_name,	site_name,	latitude,	longitude,	location_method,	location_notes,
           collection_start,	collection_end,	collection_day,	time_method,	time_notes,
           variant_string,	variant_num,	total_num)
  
  # fix dates
  get_date <- function(x) {
    s <- strsplit(as.character(x), "-")[[1]]
    if (length(s) == 1) {
      ret <- sprintf("%s-01-01", s[1])
    } else if (length(s) == 2) {
      ret <- sprintf("%s-%s-01", s[1], s[2])
    } else {
      ret <- x
    }
    as.Date(ret)
  }
  df_all <- df_all |>
    rowwise() |>
    mutate(collection_start = get_date(collection_start),
           collection_end = get_date(collection_end)) |>
    ungroup() |>
    mutate(collection_day = collection_start + (collection_end - collection_start) / 2)
  
  # fix survey_id
  df_all <- df_all |>
    group_by(site_name) |>
    mutate(period_code = sprintf("%s_%s", collection_start, collection_end),
           period_code = match(period_code, unique(period_code)),
           n_periods = max(period_code),
           period_name = ifelse(n_periods == 1, "", sprintf("_period%s", period_code)),
           survey_id = sprintf("%s_%s%s", study_id_root, site_name, period_name),
           survey_id = gsub("\\s+", "_", survey_id),
           survey_id = gsub("-", "_", survey_id),
           survey_id = gsub("’", "", survey_id)) |>
    ungroup()
  
  # split multi-locus variants if needed
  multi_variants <- mapply(function(x) {
    s1 <- x[1]
    s2 <- strsplit(x[2], "_")[[1]]
    s3 <- strsplit(gsub("_", "", x[3]), "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]]
    sprintf("%s:%s:%s", s1, s2, s3)
  }, strsplit(df_all$variant_string, ":"), SIMPLIFY = FALSE)
  m_length <- mapply(length, multi_variants)
  if (any(m_length > 1)) {
    df_all <- df_all[rep(seq_along(m_length), times = m_length),]
    df_all$variant_string <- unlist(multi_variants)
  }
  
  # fix variants
  df_all <- df_all |>
    separate(variant_string, into = c("gene", "pos", "AA"), sep = ":") |>
    filter(gene %in% c("CRT", "MDR1", "K13")) |>
    filter(!((gene == "CRT") & (pos != 76)),
           !((gene == "MDR1") & (pos != 86)),
           !((gene == "K13") & !(pos %in% k13_dictionary$CODON))) |>
    mutate(variant_string = sprintf("%s:%s:%s", tolower(gene), pos, AA)) |>
    arrange(survey_id, variant_string)
  
  # make survey data frame
  df_surveys <- df_all |>
    group_by(study_id, survey_id, country_name, site_name, latitude, longitude, location_method,
             location_notes, collection_start, collection_end, collection_day, time_method, time_notes) |>
    summarise(.groups = "drop")
  
  # make counts data frame
  df_counts <- df_all |>
    select(study_id, survey_id, collection_start, collection_end, variant_string, variant_num, total_num)
  
  # read in STAVE xlsx
  stave_xl_path <- here("analysis", "data-geoff", project_dir, sprintf("%s_STAVE.xlsx", project_dir))
  df_studies <- read_excel(stave_xl_path, sheet = "studies")
  
  # overwrite with new tables
  write_xlsx(list(studies = df_studies,
                  surveys = df_surveys,
                  counts = df_counts),
             path = stave_xl_path)
}

