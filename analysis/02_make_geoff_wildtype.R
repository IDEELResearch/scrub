# 02_make_geoff_wildtype.R
#
# Author: Bob Verity
# Date: 2025-12-09
#
# Inputs:
# - data-geoff/ (subfolders and xlsx files within)
# - data-raw/k13_ref_protein_codon_dictionary.csv
#
# Outputs:
# - data-derived/geoff_to_k13_coverage.csv
#
# Purpose:
# Reads in geoff data and extracts the start and end codons sequenced for k13.
# These are used to generate the corresponding wild-type variant for each study,
# i.e. the expected AA code if no non-synonymous mutations were observed at any
# of these sequenced loci. This information is saved in the data-derived folder.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)

# read k13 dictionary and filter to target mutations
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET))

# Iterate over each extracted data file
study_folders <- list.dirs(here("analysis", "data-geoff"), recursive = FALSE) |>
  basename()

l <- list()
for (i in seq_along(study_folders)) {
  project_dir <- study_folders[i]
  #message(paste("Processing project:", project_dir))
  
  # Read in data from xlsx
  xl_path <- here("analysis", "data-geoff", sprintf("%s/%s.xlsx", project_dir, project_dir))
  study_table <-  read_excel(xl_path, sheet = "study")
  
  k13_min <- study_table |>
    filter(FIELDS == "K13_min") |>
    pull(DATA) |>
    as.numeric()
  k13_max <- study_table |>
    filter(FIELDS == "K13_max") |>
    pull(DATA) |>
    as.numeric()
  
  l[[i]] <- data.frame(project_dir = project_dir,
                       k13_min = k13_min,
                       k13_max = k13_max)
}
df_k13_cov <- bind_rows(l)

# get WT in variantstring
df_k13_cov$WT_variant <- NA
for (i in which(!is.na(df_k13_cov$k13_min))) {
  tmp <- k13_dictionary |>
    filter(CODON >= df_k13_cov$k13_min[i],
           CODON <= df_k13_cov$k13_max[i])
  df_k13_cov$WT_variant[i] <- sprintf("k13:%s:%s", paste(tmp$CODON, collapse = "_"), paste(tmp$REF, collapse = ""))
}

# save to file
write.csv(df_k13_cov, file = here("analysis", "data-derived", "geoff_to_k13_coverage.csv"), row.names = FALSE)
