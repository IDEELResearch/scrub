# 11_merge_stave.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-derived/geoff_STAVE.rds
# - data-derived/WWARN_STAVE.rds
# - data-derived/WHO_STAVE.rds
# - data-derived/pf8_STAVE.rds
#
# Outputs:
# - data-out/stave_data (time stamped)
#
# Purpose:
# Reads in the various STAVE objects and binds them together into a single
# outward-facing object.
# 
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


s_geoff <- readRDS(file = here("analysis", "data-derived", "geoff_STAVE.rds"))
s_wwarn <- readRDS(file = here("analysis", "data-derived", "WWARN_STAVE.rds"))
s_who <- readRDS(file = here("analysis", "data-derived", "WHO_STAVE.rds"))
s_pf8 <- readRDS(file = here("analysis", "data-derived", "pf8_STAVE.rds"))

s <- STAVE_object$new()
s$append_data(studies_dataframe = s_geoff$get_studies(),
              surveys_dataframe = s_geoff$get_surveys(),
              counts_dataframe = s_geoff$get_counts())
s$append_data(studies_dataframe = s_wwarn$get_studies(),
              surveys_dataframe = s_wwarn$get_surveys(),
              counts_dataframe = s_wwarn$get_counts())
s$append_data(studies_dataframe = s_who$get_studies(),
              surveys_dataframe = s_who$get_surveys(),
              counts_dataframe = s_who$get_counts())
s$append_data(studies_dataframe = s_pf8$get_studies(),
              surveys_dataframe = s_pf8$get_surveys(),
              counts_dataframe = s_pf8$get_counts())

# save combined object to file
saveRDS(s, file = here("analysis", "data-out", "stave_data_2026.03.08.rds"))

# ------------------------------------------------------------------
# Some basic stats on the final object

s

# number of surveys, vs. number of distinct lat/lon/time points
df_surveys <- s$get_surveys()
df_surveys |>
  group_by(latitude, longitude, collection_day) |>
  summarise() |>
  dim()

# get WHO positions
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET)) |>
  mutate(WT_variant = sprintf("k13:%s:%s", CODON, REF))

# get prevalence (of WT) over all positions
l <- list()
for (i in 1:nrow(k13_dictionary)) {
  message(sprintf("%s of %s", i, nrow(k13_dictionary)))
  l[[i]] <- s$get_prevalence(k13_dictionary$WT_variant[i]) |>
    select(survey_id, collection_day, denominator) |>
    mutate(target_variant = k13_dictionary$WT_variant[i],
           collection_year = year(collection_day))
}

# also get prevalence of PD
l[[22]] <- s$get_prevalence("crt:76:T") |>
  select(survey_id, collection_day, denominator) |>
  mutate(target_variant = "crt:76:T",
         collection_year = year(collection_day))
l[[23]] <- s$get_prevalence("mdr1:86:Y") |>
  select(survey_id, collection_day, denominator) |>
  mutate(target_variant = "mdr1:86:Y",
         collection_year = year(collection_day))

df_l <- bind_rows(l)

# find max denominator over all markers
df_comb <- df_l |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator))

# get total samples sequenced at any of our positions of interest
sum(df_comb$denom_max)

# get the equivalent number for k13 561 only
bind_rows(l[16]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()

# get the equivalent number for k13 622 only
bind_rows(l[20]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()

# get the equivalent number for k13 675 only
bind_rows(l[21]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()

# get the equivalent number for all k13 only
bind_rows(l[1:21]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()

# get the equivalent number for crt only
bind_rows(l[22]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()

# get the equivalent number for mdr1 only
bind_rows(l[23]) |>
  group_by(survey_id, collection_year) |>
  summarise(denom_max = max(denominator)) |>
  pull(denom_max) |>
  sum()
