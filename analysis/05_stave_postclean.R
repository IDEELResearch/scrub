# .R
#
# Author: Bob Verity
# Date: 2025-11-25
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# (this is an example header)
#
# ------------------------------------------------------------------
# PACKAGES AND FUNCTIONS

library(tidyverse)
library(here)

# load specific versions of STAVE and variantstring to avoid
# backward-compatibility issues
#remotes::install_github("mrc-ide/STAVE@v1.1.0")
#remotes::install_github("mrc-ide/variantstring@1.8.2")
library(STAVE)
library(variantstring)

# ------------------------------------------------------------------

# read in STAVE object
s <- readRDS(here("analysis/data-out/stave_data.rds"))

# extract tables
df_studies <- s$get_studies()
df_surveys <- s$get_surveys()
df_counts <- s$get_counts()

# replace info for specific studies manually
df_fix_study <- rbind.data.frame(list(study_id = "geoff_S0011JacquesMariGmmsUnpub",
                                      url_new = "https://www.medrxiv.org/content/10.1101/2025.08.07.25333205v1"))

df_studies <- df_studies |>
  left_join(df_fix_study) |>
  mutate(url = coalesce(url_new, url)) |>
  select(-url_new)

# replace all fields based on infor read in from file
df_study_info <- read.csv(here("analysis/data-raw", "stave_study_info.csv"))

df_studies <- df_studies |>
  select(study_id, url) |>
  left_join(df_study_info) |>
  select(study_id, study_name, study_type, authors, publication_year, url_clean) |>
  rename(url = url_clean)

# ------------------------------------------------------------------

# work out if each field study is: geoff, pf7k, who, wwarn
df_studies <- df_studies |>
  mutate(origin = mapply(function(x) x[[1]], strsplit(study_id, "_")))

# find duplicates
dup_table <- table(df_studies$url, df_studies$origin) |>
  as.data.frame() |>
  pivot_wider(names_from = Var2, values_from = Freq) |>
  mutate(num_origins = (geoff > 0) + (pf7k > 0) + (who > 0) + (wwarn > 0))

dup_table |>
  head()
