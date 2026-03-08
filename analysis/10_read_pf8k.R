# 10_read_pf8k.R
#
# Author: Bob Verity
# Date: 2026-03-06
#
# Inputs:
# - data-raw/pf8_raw.csv
# - data-raw/pf8_samples.csv
# - data-raw/pf8_k13_reformat.csv
# - data-raw/pf8_study_details.csv
#
# Outputs:
# - data-derived/pf8_STAVE.rds
#
# Purpose:
# Reads in the Pf8 data. Wrangles into STAVE format for crt, mdr1, and k13. Saves to file.
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
library(data.table)

# ---------------------------------------------------- o
# 1. Malariagen wrangle
# ---------------------------------------------------- o

mdf <- data.table::fread(here("analysis", "data-raw", "pf8_raw.csv")) |>
  as.data.frame()

# grab the meta information for getting lat long and sample year
meta <- read.csv(here("analysis", "data-raw", "pf8_samples.csv"))

# and join
mdf <- left_join(mdf, meta, by = "Sample")

# sort this into a consistent format that is the same as the WHO/WWWARN information

# 1. Rename annoying markers
mdf <- mdf |> janitor::clean_names()

# select all that are needed
mdf <- mdf |> select(
  sample,
  study_id = study,
  country,
  admin1 = admin_level_1,
  lat = admin_level_1_latitude,
  lon = admin_level_1_longitude,
  year,
  crt_76_k,
  mdr1_86_n,
  kelch13_349_726_ns_changes
) |> 
  rename(k13_all = kelch13_349_726_ns_changes) |>
  mutate(study_id = sprintf("MalariaGEN_%s", gsub("-", "_", study_id)))

# count studies
length(unique(mdf$study_id))

# ----------------------------
# drop non-African samples

pf8_Africa <- c("Benin", "Burkina Faso", "Cameroon", "Cote d'Ivoire", "Democratic Republic of the Congo", 
                "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Kenya", "Madagascar", 
                "Malawi", "Mali", "Mauritania", "Mozambique", "Nigeria", "Senegal", 
                "Sudan", "Tanzania", "Uganda")

mdf <- mdf |>
  mutate(country = ifelse(country == "C√¥te d'Ivoire", "Cote d'Ivoire", country)) |>
  filter(country %in% pf8_Africa)

# count studies
length(unique(mdf$study_id))

# ----------------------------
# crt 76

unique(mdf$crt_76_k)

# get into variantstring format
mdf$crt_76_k <- mapply(function(x) {
  if (x[1] == "-") {
    return(NA)
  }
  if (length(x) == 1) {
    return(sprintf("crt:76:%s", x))
  }
  return(sprintf("crt:76:%s", paste(sort(x), collapse = "/")))
}, strsplit(mdf$crt_76_k, ","))

# ----------------------------
# mdr1 86

unique(mdf$mdr1_86_n)

# get into variantstring format
mdf$mdr1_86_n <- mapply(function(x) {
  if (x[1] == "-") {
    return(NA)
  }
  if (length(x) == 1) {
    return(sprintf("mdr1:86:%s", x))
  }
  return(sprintf("mdr1:86:%s", paste(sort(x), collapse = "/")))
}, strsplit(mdf$mdr1_86_n, ","))

# ----------------------------
# k13

# when we filter to African samples and WHO-validated loci, there are only a
# very small number of unique Pf8 k13 variants. These can be converted to
# variantstring format using a simple lookup table
pf8_reformat <- read.csv(here("analysis", "data-raw", "pf8_k13_reformat.csv"), tryLogical = FALSE)
pf8_reformat$new_format <- apply(pf8_reformat, 1, function(x) {
  paste(x[-1], collapse = "")
})

k13_stem <- "k13:441_446_449_458_469_476_481_493_515_527_537_538_539_543_553_561_568_574_580_622_675"
pf8_reformat$k13_variantstring <- sprintf("%s:%s", k13_stem, pf8_reformat$new_format)

mdf <- left_join(mdf, pf8_reformat |>
            select(k13_all = old_format,
                   k13_variantstring))

# all variants that did not match the lookup can be coded as WT. However, we
# then need to go back in and ensure that missing data are correctly coded as NA
k13_WT <- sprintf("%s:PFGNCMAYRPNGRIPRVPCRA", k13_stem)
mdf <- mdf |>
  mutate(k13_variantstring = ifelse(is.na(k13_variantstring), k13_WT, k13_variantstring),
         k13_variantstring = ifelse(k13_all %in% c("-", "!", "*"), NA, k13_variantstring))

# ----------------------------
# get into long format over genes and aggregate samples

df_long <- mdf |>
  select(-k13_all) |>
  pivot_longer(cols = c(crt_76_k, mdr1_86_n, k13_variantstring),
               names_to = "gene", values_to = "variant_string") |>
  filter(!is.na(variant_string)) |>
  group_by(study_id, country, admin1, lat, lon, year, gene, variant_string) |>
  summarise(variant_num = n()) |>
  group_by(study_id, country, admin1, lat, lon, year, gene) |>
  mutate(total_num = sum(variant_num)) |>
  ungroup() |>
  select(-gene)

# create survey_id from study_id, admin1 and year
df_long <- df_long |>
  mutate(survey_id = sprintf("%s_%s_%s", study_id, janitor::make_clean_names(admin1, allow_dupes = TRUE), year))

# check that this method of producing survey_id corresponds to unique lat/lon
# (would fail if there were multiple sites in the same admin1 in the same year)
df_long |>
  group_by(survey_id) |>
  mutate(latlon = sprintf("%s_%s", lat, lon)) |>
  summarise(n = length(unique(latlon))) |>
  pull(n)

# ----------------------------
# get into STAVE format

# read in details of individual Pf8 studies
pf8_details <- read.csv(here("analysis", "data-raw", "pf8_study_details.csv"))

# counts data frame
df_counts <- df_long |>
  select(study_id, survey_id, variant_string, variant_num, total_num) |>
  mutate(notes = "")

# surveys data frame
df_surveys <- df_long |>
  select(study_id,
         survey_id,
         country_name = country,
         site_name = admin1,
         latitude = lat,
         longitude = lon,
         year = year) |>
  mutate(location_method = NA,
         location_notes = NA,
         collection_start = NA,
         collection_end = NA,
         collection_day = as.Date(sprintf("%s-07-01", year)),
         time_method = "Collection day imputed as midpoint of year",
         time_notes = NA) |>
  group_by(study_id, survey_id, country_name, site_name, latitude, longitude, location_method, location_notes,
         collection_start, collection_end, collection_day, time_method, time_notes) |>
  summarise() |>
  ungroup()

# studies data frame
df_studies <- df_long |>
  group_by(study_id) |>
  summarise() |>
  left_join(pf8_details, by = "study_id") |>
  mutate(description = NA,
         access_level = "public",
         reference = "https://doi.org/10.12688/wellcomeopenres.24031.1",
         reference_year = 2025,
         PMID = "41267740") |>
  select(study_id, study_label, description, access_level, contributors, reference,
         reference_year, PMID)

# make STAVE object and append data
s <- STAVE::STAVE_object$new()
s$append_data(studies_dataframe = df_studies,
              surveys_dataframe = df_surveys,
              counts_dataframe = df_counts)
s

# write to file
saveRDS(s, file = here("analysis", "data-derived", "pf8_STAVE.rds"))

