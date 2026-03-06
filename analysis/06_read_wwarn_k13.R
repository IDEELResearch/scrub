# 06_read_wwarn_k13.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-raw/WWARN_K13_database_04-12-2023.xls
# - data-derived/geoff_STAVE.rds
# - data-raw/PMID_k13_replace.xlsx
# - data-raw/paper_info.csv
# - data-derived/PMID_to_k13_coverage.csv
# - data-raw/k13_ref_protein_codon_dictionary.csv
#
# Outputs:
# - data-derived/wwarn_k13_clean.rds
#
# Purpose:
# Reads in WWARN k13 data. Filters, deduplicates, excludes PMIDs already covered
# by GEOFF data. Wrangles into format needed by STAVE. Saves data object to
# file, ready to combine with partner drug info.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)

# ------------------------------------------------------------------------
# read in data and apply initial filters

wwarn_k13 <- readxl::read_xls(here("analysis", "data-raw", "WWARN_K13_database_04-12-2023.xls"), sheet = 1) |>
  mutate(tested = as.numeric(tested),
         present = as.numeric(present),
         year = as.numeric(year),
         lat = as.numeric(lat),
         lon = as.numeric(lon),
         PMID = as.numeric(pubMedId),
         notes = as.character(NA))

# count studies
length(unique(wwarn_k13$authors))

# filter to Africa
wwarn_k13 <- wwarn_k13 |>
  filter(continent == "Africa") |>
  arrange(country, site, year, mutation) |>
  select(country, site, year, mutation, present, tested, lat, lon, PMID, authors)

# count studies
length(unique(wwarn_k13$authors))

# drop studies with no PMID
wwarn_k13 <- wwarn_k13 |>
  filter(!is.na(PMID)) |>
  select(-authors)

# count studies
length(unique(wwarn_k13$PMID))

# drop studies already in GEOFF
s <- readRDS(here("analysis", "data-derived", "geoff_STAVE.rds"))
PMID_geoff <- s$get_studies() |>
  filter(!is.na(PMID)) |>
  pull(PMID)

wwarn_k13 <- wwarn_k13|>
  filter(!(PMID %in% PMID_geoff))

# count studies
length(unique(wwarn_k13$PMID))

# ------------------------------------------------------------------------
# remove some studies and fix errors in others

# remove duplicate entries
wwarn_k13 <- wwarn_k13 |>
  distinct(PMID, country, site, year, mutation, .keep_all = TRUE)

# drop some studies because these are frequencies based on reads, not prevalence
wwarn_k13 <- wwarn_k13 |>
  filter(!(PMID %in% c(25180240, 31591113)))

# replace some studies with data re-entered manually due to mistakes
wwarn_k13 <- wwarn_k13 |>
  filter(PMID != 25918205) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "25918205")) |>
  filter(PMID != 26667053) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "26667053")) |>
  filter(PMID != 28039354) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "28039354")) |>
  filter(PMID != 28797235) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "28797235")) |>
  filter(PMID != 29436339) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "29436339")) |>
  filter(PMID != 29582728) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "29582728")) |>
  filter(PMID != 30559133) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "30559133")) |>
  filter(PMID != 31132213) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "31132213")) |>
  filter(PMID != 32795367) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "32795367")) |>
  filter(PMID != 33789667) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "33789667"))

# replace some bad characters in site names
wwarn_k13 <- wwarn_k13 |>
  mutate(site = gsub("√©", "e", site),
         site = gsub("é", "e", site),
         site = gsub(",", "", site))

# add study and survey IDs
wwarn_k13 <- wwarn_k13 |>
  mutate(study_id = sprintf("WWARN_%s", PMID),
         survey_id = sprintf("%s_%s_%s", study_id, site, year),
         survey_id = gsub("-", "_", survey_id),
         survey_id = gsub(" ", "_", survey_id),
         survey_id = gsub("'", "_", survey_id))

# ------------------------------------------------------------------------

# link to data on paper based on PMID
paper_info <- read.csv(here("analysis", "data-raw", "paper_info.csv")) |>
  filter(!duplicated(PMID))

wwarn_k13 <- wwarn_k13 |>
  left_join(paper_info, by = join_by(PMID))

# get mutations in variantstring
wwarn_k13$variant_string <- mapply(function(x) {
  if (x == "wildtype") {
    return(NA)
  }
  x <- strsplit(x, "&")[[1]]
  x_parts <- str_match(x, "^([A-Za-z])(\\d+)([A-Za-z/]+)$")
  sprintf("k13:%s:%s", paste(x_parts[,3], collapse = "_"), paste(x_parts[,4], collapse = "_"))
}, wwarn_k13$mutation)

# read in WT for each PMID (calculated based on coverage)
k13_coverage <- read.csv(here("analysis", "data-derived", "PMID_to_k13_coverage.csv"))

wwarn_k13 <- wwarn_k13 |>
  left_join(k13_coverage |>
              select(PMID, WT_variant),
            by = join_by(PMID))

# count studies
length(unique(wwarn_k13$PMID))

# drop studies for which WT is NA, implying that none of the target loci are covered
wwarn_k13 <- wwarn_k13 |>
  filter(!is.na(WT_variant))

# count studies
length(unique(wwarn_k13$PMID))

# overlay variant on top of WT
wwarn_k13 <- wwarn_k13 |>
  mutate(variant_string_combined = overlay_variant(variant_string, WT_variant))

# ------------------------------------------------------------------------

# get each variant back into single-position format
l <- wwarn_k13 |>
  rowwise() |>
  group_split()
for (i in seq_along(l)) {
  vl <- variant_to_long(l[[i]]$variant_string_combined)[[1]]
  vs <- vl |>
    group_by(gene, pos) |>
    summarise(aa = paste(aa, collapse = "/"),
              .groups = "drop") |>
    ungroup() |>
    mutate(variant_string = sprintf("%s:%s:%s", gene, pos, aa)) |>
    pull(variant_string)
  l[[i]] <- l[[i]][rep(1, length(vs)),]
  l[[i]]$variant_string <- vs
}
wwarn_k13 <- bind_rows(l) |>
  select(study_id, survey_id, country, site, year, variant_string, present, tested, lat, lon, PMID,
         notes, url_clean, study_name, study_type, authors, publication_year)

# group and sum numerators for same mutation
wwarn_k13 <- wwarn_k13 |>
  group_by(study_id, survey_id, variant_string) |>
  summarise(present = sum(present, na.rm = TRUE),
            tested = tested[1],
            across(-c(present, tested), first),
            .groups = "drop")

# filter to WHO codons
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET))

wwarn_k13 <- wwarn_k13 |>
  mutate(variant_string_split = variant_string) |>
  separate(variant_string_split, into = c("gene", "pos", "AA"), sep = ":") |>
  filter(pos %in% k13_dictionary$CODON) |>
  select(-c(gene, pos, AA))

# ------------------------------------------------------------------------

# replace Asua et al. data
wwarn_k13 <- wwarn_k13 |>
  filter(PMID != 33146722)

df_asua <- readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "33146722")

# merge with data on paper and tidy up to match wwarn_k13 format 
asua_info <- paper_info |>
  filter(PMID == 33146722)

df_asua <- df_asua |>
  mutate(study_id = "WWARN_33146722",
         survey_id = sprintf("%s_%s_%s", study_id, site, year),
         url_clean = asua_info$url_clean,
         study_name = asua_info$study_name,
         study_type = "peer_reviewed",
         authors = asua_info$authors,
         publication_year = asua_info$publication_year)

wwarn_k13 <- bind_rows(wwarn_k13, df_asua)

# count studies
length(unique(wwarn_k13$PMID))

# ------------------------------------------------------------------------
# import further cleaned datasets

import_PMIDs <- c(28594879, 35144535, 33107096, 26483118, 31296223, 33864801, 33350925,
                  27573632, 34216470, 32747827, 32822392, 34551228, 35477399, 31427297,
                  31034031, 29458380, 32103124)

wwarn_k13 <- wwarn_k13 |>
  filter(!(PMID %in% import_PMIDs))
for (i in seq_along(import_PMIDs)) {
  wwarn_k13 <- wwarn_k13 |>
    bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"),
                                 sheet = as.character(import_PMIDs[i])))
}

# count studies
length(unique(wwarn_k13$PMID))

# ------------------------------------------------------------------------

# save to file
saveRDS(wwarn_k13, file = here("analysis", "data-derived", "wwarn_k13_clean.rds"))

