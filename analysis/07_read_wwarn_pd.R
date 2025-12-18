# 07_read_wwarn_pd.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-raw/WWARN_partnerdrug_database_04-12-2023.xls
# - data-raw/WWARN_pd_PMID.csv
# - data-raw/PMID_pd_replace.xlsx
# - data-raw/paper_info.csv
#
# Outputs:
# - data-derived/wwarn_pd_clean.rds
#
# Purpose:
# Reads in WWARN partner drug data (crt and mdr1). Filters, deduplicates, excludes PMIDs already covered
# by GEOFF data. Wrangles into format needed by STAVE. Saves data object to
# file, ready to combine with k13 info.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
library(janitor)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)

# ------------------------------------------------------------------

# read in data and tidy up
wwarn_pd <- readxl::read_xls(here("analysis", "data-raw", "WWARN_partnerdrug_database_04-12-2023.xls"), 
                              sheet = 1) |>
  janitor::clean_names() |>
  mutate(tested = as.numeric(tested),
         present = as.numeric(present),
         study_start = as.numeric(study_start),
         study_end = as.numeric(study_end),
         lat = as.numeric(lat),
         lon = as.numeric(lon),
         url = publication_url) |>
  select(study_id, country, site, lon, lat, study_start, study_end, marker_type,
         tested, present, mixed_present, url)

# filter to Africa countries
wwarn_pd <- wwarn_pd |>
  filter(country %in% c("Angola", "Benin", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
                        "Chad", "Comoros", "Congo", "Côte d'Ivoire", "Democratic Republic of the Congo",
                        "Sao Tome and Principe", "Equatorial Guinea", "Ethiopia", "Gabon", "Ghana", "Guinea",
                        "Guinea-Bissau", "Kenya", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania",
                        "Mozambique", "Niger", "Nigeria", "Rwanda", "Senegal", "Sierra Leone", "South Africa",
                        "South Sudan", "Sudan", "Swaziland", "Tanzania", "Gambia", "Togo", "Uganda", "Zambia",
                        "Zimbabwe", "Eritrea", "South sudan", "Libya", "Djibouti", "Algeria", "Botswana", 
                        "Somalia"))

# merge with data on PMID
df_PMID <- read.csv(here("analysis", "data-raw", "WWARN_pd_PMID.csv"))
wwarn_pd <- wwarn_pd |>
  left_join(df_PMID, by = join_by(url))

# remove studies with no PMID
wwarn_pd <- wwarn_pd |>
  filter(!is.na(PMID))

# remove specific PMIDs
wwarn_pd <- wwarn_pd |>
  filter(PMID != 26236581) |> # data has some issues (numerator doesn't sum to denominator) and the original paper does not provide sufficient detail to correct
  filter(PMID != 22453078) |> # not enough information in original paper to work out which mixed infections are present in Ghana
  filter(PMID != 22208458) |> # not enough information in original paper to work out the numbers of different mutations
  filter(PMID != 20199676) |> # issues with data, and cannot verify because original paper does not give breakdown by site
  filter(PMID != 19718439) |> # issues with data, and cannot verify because original paper only has values in charts
  filter(PMID != 17313506) # issues with data, and cannot verify because original paper only has values in charts

# replace site with country if NA
wwarn_pd <- wwarn_pd |>
  mutate(site = ifelse(is.na(site), country, site))

# make study and survey IDs
wwarn_pd <- wwarn_pd |>
  left_join(data.frame(site = unique(wwarn_pd$site)) |>
              mutate(site_clean = janitor::make_clean_names(site)),
            by = join_by(site))  |>
  mutate(study_id = sprintf("WWARN_%s", PMID),
         survey_id = sprintf("%s_%s_%s", study_id, site_clean, study_start))

# filter markers and exclude duplicates
wwarn_pd <- wwarn_pd |>
  filter(marker_type %in% c("pfmdr1 N86", "pfmdr1 86Y", "pfmdr1 86N/Y",
                            "pfcrt K76", "pfcrt 76T", "pfcrt 76K/T")) |>
  distinct(survey_id, marker_type, .keep_all = TRUE)

# replace marker with variantstring format
wwarn_pd <- wwarn_pd |>
  left_join(data.frame(marker_type = c("pfmdr1 N86", "pfmdr1 86Y", "pfmdr1 86N/Y",
                                       "pfcrt K76", "pfcrt 76T", "pfcrt 76K/T"),
                       variant_string = c("mdr1:86:N", "mdr1:86:Y", "mdr1:86:N/Y",
                                          "crt:76:K", "crt:76:T", "crt:76:K/T")),
            by = join_by(marker_type)) |>
  select(study_id, survey_id, country, site, lon, lat, study_start, study_end,
         variant_string, present, tested, PMID)

# replace some values manually from file where there are mistakes in the original
sheet_names <- excel_sheets(here("analysis", "data-raw", "PMID_pd_replace.xlsx"))
for (i in seq_along(sheet_names)) {
  wwarn_pd <- wwarn_pd |>
    filter(PMID != as.numeric(sheet_names[i])) |>
    bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_pd_replace.xlsx"), sheet = sheet_names[i]))
}

# merge with info on paper
df_paper <- read.csv(here("analysis", "data-raw", "paper_info.csv")) |>
  distinct(PMID, .keep_all = TRUE)
wwarn_pd <- wwarn_pd |>
  left_join(df_paper,
            by = join_by(PMID))

# ------------------------------------------------------------------------

# tidy up and save to file
wwarn_pd <- wwarn_pd |>
  select(-c(url, doi))

saveRDS(wwarn_pd, file = here("analysis", "data-derived", "wwarn_pd_clean.rds"))
