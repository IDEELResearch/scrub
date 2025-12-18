# 11_merge_stave.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-derived/geoff_STAVE.rds
# - data-derived/WWARN_STAVE.rds
# - data-derived/WHO_STAVE.rds
# - data-derived/pf7_STAVE.rds
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
s_pf7 <- readRDS(file = here("analysis", "data-derived", "pf7_STAVE.rds"))

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
s$append_data(studies_dataframe = s_pf7$get_studies(),
              surveys_dataframe = s_pf7$get_surveys(),
              counts_dataframe = s_pf7$get_counts())

# save combined object to file
saveRDS(s, file = here("analysis", "data-out", "stave_data_2025.12.17.rds"))

# ------------------------------------------------------------------------
# sanity check plots

# get prevalence of a target mutation
p <- s$get_prevalence("k13:675:V")

# split studies into sources
p$source <- mapply(function(x) {
  match(x[1], c("WWARN", "WHO", "pf7"))
}, strsplit(p$study_id, "_"))
p$source[is.na(p$source)] <- 4
p$source <- c("WWARN", "WHO", "pf7", "GEOFF")[p$source]

# plot samples over time
p |>
  mutate(collection_year = year(collection_day)) |>
  group_by(collection_year, source) |>
  summarise(n = sum(denominator)) |>
  ggplot() + theme_bw() +
  geom_col(aes(x = collection_year, y = n, fill = source)) +
  labs(x = "Collection Year", y = "Total Samples Sequenced")

