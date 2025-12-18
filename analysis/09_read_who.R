# 09_read_who.R
#
# Author: Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-raw/WHO_res_database_02-01-2024.xlsx
# - data-raw/paper_info.csv
# - data-derived/geoff_STAVE.rds
# - data-derived/WWARN_STAVE.rds
# - data-raw/PMID_k13_replace.xlsx
# - data-derived/PMID_to_k13_coverage.csv
# - data-raw/k13_ref_protein_codon_dictionary.csv
#
# Outputs:
# - data-derived/WHO_STAVE.rds
#
# Purpose:
# Reads in WHO data, which ultimately only contains info on k13 once filters
# have been applied. Deduplicates, fixes data issues, and wrangles into STAVE
# format. Then saves to file.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)

# ------------------------------------------------------------------------

# Read in WHO data and merge
who_df <- readxl::read_xlsx(here("analysis", "data-raw", "WHO_res_database_02-01-2024.xlsx"), sheet = 2)
who_res <- readxl::read_xlsx(here("analysis", "data-raw", "WHO_res_database_02-01-2024.xlsx"), sheet = 3)
who <- left_join(who_res, who_df, by = "ID") |>
  mutate(LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE))

# Subset to African samples
who_countries <- c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Côte d’Ivoire", 
                   "Cabo Verde", "Cameroon", "Central African Republic", "Chad", 
                   "Comoros", "Congo", "Democratic Republic of the Congo", "Equatorial Guinea", 
                   "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", 
                   "Guinea-Bissau", "Kenya", "Liberia", "Madagascar", "Malawi", 
                   "Mali", "Mauritania", "Mayotte", "Mozambique", "Niger", "Nigeria", 
                   "Rwanda", "Sao Tome and Principe", "Senegal", "Sierra Leone", 
                   "Somalia", "South Africa", "South Sudan", "Sudan", "Togo", "Uganda", 
                   "United Republic of Tanzania", "Zambia", "Zimbabwe")

who <- who |>
  filter(COUNTRY_NAME %in% who_countries)

# Filter to k13, as mdr1 is copy number info only, and no crt left after previous filters
who <- who |>
  filter(MM_TYPE == "Pfkelch13")

# filter to those with citation
who <- who |>
  filter(!is.na(CITATION_URL))

# merge with paper info
paper_info <- read.csv(here("analysis", "data-raw", "paper_info.csv"))

who <- who |>
  rename(url = CITATION_URL) |>
  left_join(paper_info, by = join_by(url))

# drop studies already in GEOFF
s <- readRDS(here("analysis", "data-derived", "geoff_STAVE.rds"))
PMID_geoff <- s$get_studies()$PMID
PMID_geoff <- PMID_geoff[!is.na(PMID_geoff)]

who <- who |>
  filter(!(PMID %in% PMID_geoff))

# drop studies already in WWARN
s <- readRDS(here("analysis", "data-derived", "WWARN_STAVE.rds"))
PMID_wwarn <- s$get_studies()$PMID

who <- who |>
  filter(!(PMID %in% PMID_wwarn))

# get variant_num from proportion
who <- who |>
  mutate(variant_num = round(as.numeric(PROPORTION) / 100 * SAMPLE_SIZE))

# ------------------------------------------------------------------------
# manual fixes to data

who <- who |>
  mutate(notes = NA) |>
  filter(PMID != 28069653) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "28069653")) |>
  filter(PMID != 33549122) |>
  bind_rows(readxl::read_excel(here("analysis", "data-raw", "PMID_k13_replace.xlsx"), sheet = "33549122"))

# replace some bad characters in site names
who <- who |>
  mutate(SITE_NAME = gsub("√©", "e", SITE_NAME),
         SITE_NAME = gsub("é", "e", SITE_NAME),
         SITE_NAME = gsub(",", "", SITE_NAME))

# add study and survey IDs
who <- who |>
  mutate(study_id = sprintf("WHO_%s", PMID),
         survey_id = sprintf("%s_ID%s", study_id, ID),
         survey_id = gsub("-", "_", survey_id),
         survey_id = gsub(" ", "_", survey_id),
         survey_id = gsub("'", "_", survey_id))

# ------------------------------------------------------------------------

# get mutations in variantstring
who$variant_string <- mapply(function(x) {
  if (x == "WT") {
    return(NA)
  }
  x <- strsplit(x, "&")[[1]]
  x_parts <- str_match(x, "^([A-Za-z])(\\d+)([A-Za-z/]+)$")
  sprintf("k13:%s:%s", paste(x_parts[,3], collapse = "_"), paste(x_parts[,4], collapse = "_"))
}, who$GENOTYPE)

# read in WT for each PMID (calculated based on coverage)
k13_coverage <- read.csv(here("analysis", "data-derived", "PMID_to_k13_coverage.csv"))

who <- who |>
  left_join(k13_coverage |>
              select(PMID, WT_variant),
            by = join_by(PMID))

# drop studies for which WT is NA, implying that none of the target loci are covered
who <- who |>
  filter(!is.na(WT_variant))

# overlay variant on top of WT
who <- who |>
  mutate(variant_string_combined = overlay_variant(variant_string, WT_variant))

# ------------------------------------------------------------------------

# get each variant back into single-position format
l <- who |>
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
who <- bind_rows(l) |>
  select(study_id, survey_id, COUNTRY_NAME, SITE_NAME, YEAR_START, variant_string, variant_num,
         SAMPLE_SIZE, LATITUDE, LONGITUDE, PMID, notes, url_clean, study_name, study_type,
         authors, publication_year)

# group and sum numerators for same mutation
who <- who |>
  group_by(study_id, survey_id, variant_string) |>
  summarise(variant_num = sum(variant_num, na.rm = TRUE),
            total_num = SAMPLE_SIZE[1],
            across(-c(variant_num, SAMPLE_SIZE), first),
            .groups = "drop")

# filter to WHO codons
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET))

who <- who |>
  mutate(variant_string_split = variant_string) |>
  separate(variant_string_split, into = c("gene", "pos", "AA"), sep = ":") |>
  filter(pos %in% k13_dictionary$CODON) |>
  select(-c(gene, pos, AA))

# ------------------------------------------------------------------------
# make STAVE objects

# study-level data frame
df_studies <- who |>
  group_by(study_id) |>
  summarise(study_label = study_name[1],
            contributors = authors[1],
            reference = url_clean[1],
            reference_year = publication_year[1],
            PMID = PMID[1]) |>
  ungroup() |>
  mutate(description = NA,
         access_level = "public") |>
  select(study_id, study_label, description, access_level, contributors, reference, reference_year, PMID)

# survey-level data frame
df_surveys <- who |>
  group_by(study_id, survey_id, COUNTRY_NAME, SITE_NAME, YEAR_START) |>
  summarise(latitude = LATITUDE[1],
            longitude = LONGITUDE[1],
            .groups = "drop") |>
  mutate(location_method = "WHO coordinates",
         location_notes = NA,
         collection_start = as.Date(sprintf("%s-01-01", YEAR_START)),
         collection_end = as.Date(sprintf("%s-12-31", YEAR_START)),
         collection_day = as.Date(sprintf("%s-07-01", YEAR_START)),
         time_method = "Midpoint of WHO recorded year",
         time_notes = NA) |>
  rename(country_name = COUNTRY_NAME,
         site_name = SITE_NAME) |>
  select(study_id, survey_id, country_name, site_name, latitude, longitude, location_method, location_notes,
         collection_start, collection_end, collection_day, time_method, time_notes)

# counts-level data frame
df_counts <- who |>
  mutate(notes = notes) |>
  select(study_id, survey_id, variant_string, variant_num, total_num, notes)

# make STAVE object and append data
s <- STAVE::STAVE_object$new()
s$append_data(studies_dataframe = df_studies,
              surveys_dataframe = df_surveys,
              counts_dataframe = df_counts)
s

# write to file
saveRDS(s, file = here("analysis", "data-derived", "WHO_STAVE.rds"))
