# Load required libraries
library(tidyverse)
library(lubridate)
library(here)
library(purrr)
devtools::load_all()

################################################################################
#
# Script: 02_clean_geoff.R
# Purpose: Perform various data cleaning to allow deduplication and combining
# with other data sources within {stave} object
# 
################################################################################

# read in data
# Load the combined geoff data table created in the first script
master_table <- readRDS(here("analysis", "data-derived", 
                             "geoff_res.rds"))

# Load mutation key for k13 reference ranges
mutation_key <- read.csv(here("analysis", "data-raw", 
                              "k13_ref_protein_codon_dictionary.csv"))

## start by looking at the data - explore this in the console
master_unique <- lapply(master_table, unique)

# Initial main clean
master_table_clean <- master_table |>
  dplyr::mutate(substudy = gsub(" ","", substudy)) |>
  dplyr::mutate(substudy = gsub("_","", substudy)) |>
  dplyr::mutate(substudy = if_else(substudy == "untreadextracted",
                                   "untreatedextracted", substudy)) |>
  dplyr::mutate(substudy = if_else(substudy == "extracted",
                                   "untreatedextracted", substudy)) |>
  dplyr::mutate(substudy = if_else((substudy == "day0" & study_uid == "s0029_warsame_2019"),
                                   "day0extracted", substudy)) |>
  dplyr::mutate(
    lat_n = if_else(study_uid == "s0020_some_2024", "11.79479482", lat_n),
    lon_e = if_else(study_uid == "s0020_some_2024", "-2.919282118", lon_e)
  ) |>
  # OJ: I checked this and there is no issue here - suggest deleting
  # mutate(lat_n = as.numeric(if_else(grepl("^‚àí", lat_n), 
  #                                   gsub("^‚àí", "-", lat_n), 
  #                                   lat_n))) |> # does not appear to work - Ronald will need to fix manually in S0147Jeang2024 entry
  dplyr::mutate(site_name = gsub(" ", "_", tolower(site_name))) |>
  dplyr::mutate(collection_location = gsub(" ", "_", tolower(collection_location))) |>
  dplyr::mutate(collection_location = gsub(",", "", tolower(collection_location))) |>
  dplyr::mutate(country = stringr::str_to_title(country)) |>
  dplyr::mutate(country = if_else(country %in% c("Democratic Republic Of The Congo", "Drc"),
                                  "Democratic Republic of the Congo",
                                  country)) |>
  dplyr::mutate(country = if_else(country == "Rukara", "Rwanda", country)) |>
  dplyr::mutate(pretreatment_samples = if_else(pretreatment_samples == "",
                                               NA, pretreatment_samples)) |>
  dplyr::mutate(study_design_age_min_years = if_else(study_design_age_min_years == "na",
                                                     NA, study_design_age_min_years)) |>
  dplyr::mutate(study_design_age_min_years = if_else(study_design_age_min_years == "1year",
                                                     "1", study_design_age_min_years)) |>
  dplyr::mutate(study_design_age_max_years = if_else(study_design_age_max_years == "na",
                                                     NA, study_design_age_max_years)) |>
  dplyr::mutate(site_study_type = if_else(study_uid == "s0034_uwimana", 
                                          "TES", site_study_type)) |>
  dplyr::mutate(site_study_type = toupper(gsub("_","", site_study_type))) |>
  dplyr::mutate(publication_status = gsub(" ", "_", publication_status)) |>
  dplyr::mutate(publication_status = if_else(publication_status == "na", NA, publication_status)) |>
  dplyr::mutate(publication_status = if_else(str_detect(study_uid, "unpub"), "unpublished", publication_status)) |>
  dplyr::mutate(publication_status = if_else(publication_status == "draft_paper", "unpublished", publication_status)) |>
  dplyr::mutate(publication_status = if_else(study_uid %in% c("s0031_bergmann_2021", #GCD: manually checked that these were pub'd
                                                              "s0057_koko", 
                                                              "s0059_osborne"),
                                             "peer_reviewed", publication_status)) |>
  dplyr::mutate(publication_status = if_else(publication_status == "published" & pmid == "37420265", "peer_reviewed", publication_status)) |>
  dplyr::mutate(publication_year = if_else(publication_status == "unpublished", NA, publication_year)) |>
  dplyr::mutate(countries_covered = gsub(",", ", ", countries_covered)) |>
  dplyr::mutate(countries_covered = stringr::str_to_title(countries_covered)) |>
  dplyr::mutate(countries_covered = if_else(countries_covered %in% c("Democratic Republic Of The Congo", "Drc"),
                                            "Democratic Republic of the Congo",
                                            country)) |>
  dplyr::mutate(data_processing_pipeline = toupper(data_processing_pipeline)) |>
  dplyr::mutate(first_author_surname = stringr::str_to_title(first_author_surname)) |>
  dplyr::mutate(iso3c = replace(iso3c, iso3c == "DRC", "COD")) |> 
  dplyr::mutate(iso3c = replace(iso3c, iso3c == "ERT", "ERI")) |>
  dplyr::mutate(country = replace(country, country == "Gqn", "Gnq")) |>
  dplyr::mutate(iso3c = replace(iso3c, iso3c == "GQN", "GNQ")) %>% 
  dplyr::mutate(site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "32459360", "COMMUNITY")) %>% 
  dplyr::mutate(site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "37670357", "COMMUNITY"))

# additional cleaning to correct for malformed iso3cs
master_table_clean <- master_table_clean %>% 
  mutate(country = if_else(country %in% c("Ago", "Ben", "Cmr", "Caf", "Cog", "Cod", "Eth", "Gha", "Stp", "Sen", "Sdn", "Uga", "Gnq", "Gab", "Mli", "Nga", "Sle", "Som", "Zmb", "Zaf"),
                           suppressWarnings(countrycode::countrycode(toupper(country), "iso3c", "country.name.en")),
                           country))

# additional cleaning to correct for malformed pretreatment_samples
master_table_clean <- master_table_clean %>% 
  mutate(pretreatment_samples = replace(pretreatment_samples, pretreatment_samples %in% c("pre_treatment", "Y"), "yes"))

# additional cleaning to correct for malformed pretreatment_samples
master_table_clean <- master_table_clean %>% 
  mutate(publication_year = replace(publication_year, publication_year %in% c("none"), NA))

# perform tests
african_countries <- data.frame(country = countrycode::codelist$country.name.en,
                                continent = countrycode::codelist$continent,
                                iso3c = countrycode::codelist$iso3c) |>
                      dplyr::filter(continent == "Africa") |> dplyr::distinct()
allowed_countries <- c(african_countries$country, "Democratic Republic of the Congo")
allowed_iso3c <- c(african_countries$iso3c)
allowed_substudy <- c("untreatedextracted", "untreatedcalculated",
                      "treatedextracted", "treatedcalculated",
                      "day0extracted", "day0calculated", 
                      "tesday0treatmentfailure","day0reinfectionextracted",
                      "day0recrudescenceextracted") # 
allowed_pretreatment <- c("no", "yes", NA)
allowed_publication <- c("peer_reviewed","preprint", "unpublished", "data_incomplete")
allowed_pub_year <- c(as.character(seq(from = 2000, to = (lubridate::year(Sys.Date())))), NA) # automatically update as years change
allowed_site_types <- c("HEALTHFACILITY", "COMMUNITY", "TES", "CCS", "DHS", "REFUGEERECEPTIONCENTER", "UNKNOWN", "SURVEILLANCE")

# run tests
check_values_in_column(master_table_clean, "country", allowed_countries)
check_values_in_column(master_table_clean, "iso3c", allowed_iso3c)
check_values_in_column(master_table_clean, "substudy", allowed_substudy)
check_values_in_column(master_table_clean, "pretreatment_samples", allowed_pretreatment)
check_values_in_column(master_table_clean, "publication_status", allowed_publication)
check_values_in_column(master_table_clean, "publication_year", allowed_pub_year)
check_values_in_column(master_table_clean, "site_study_type", allowed_site_types)

# Clean mutation names

# Expand all gene mutation ranges for reference range syntax
indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))

master_table_clean$gene_mutation[indices_to_transform] <- sapply(
  master_table_clean$gene_mutation[indices_to_transform],
  collapse_k13_range, 
  mutation_key = mutation_key
)

# TODO: Check with Isabela about the s0054 extraction for 1034
master_table_clean$gene_mutation[master_table_clean$gene_mutation == "MDR1:1034:FC"] <- "MDR1:1034:C"

# Adjust and add dates and survey IDs
convstart <- adjust_invalid_date(unique(master_table_clean$date_start), is_start = TRUE)
master_table_clean$collection_start <- convstart[match(master_table_clean$date_start, unique(master_table_clean$date_start))]
convend <- adjust_invalid_date(unique(master_table_clean$date_end), is_start = FALSE)
master_table_clean$collection_end <- convend[match(master_table_clean$date_end, unique(master_table_clean$date_end))]

# Compute collection day
master_table_clean <- master_table_clean %>% rowwise() %>% 
  mutate(collection_day = median(c(collection_start, collection_end))) %>% 
  ungroup()

# Fix long and lat which got read in as formula vs values
master_table_clean$lat_n <- clean_excel_formulas(master_table_clean$lat_n)
master_table_clean$lon_e <- clean_excel_formulas(master_table_clean$lon_e)

# Correct gene mutation format
convmutation <- format_variants_for_stave(unique(master_table_clean$gene_mutation))
master_table_clean$gene_mutation <- convmutation[match(master_table_clean$gene_mutation, unique(master_table_clean$gene_mutation))]

# change variable classes
master_table_clean <- master_table_clean |>
  dplyr::mutate(collection_day = as.Date(collection_day),
                collection_start = as.Date(collection_start),
                collection_end = as.Date(collection_end),
                lat_n = as.numeric(lat_n),
                lon_e = as.numeric(lon_e),
                mutant_num = as.numeric(mutant_num),
                total_num = as.numeric(total_num),
                publication_year = as.numeric(publication_year)) |>
  dplyr::mutate(gene = str_extract(gene_mutation, "^[^:]+"),
                mut = str_extract(gene_mutation, "(?<=:).*")) |>
  dplyr::filter(total_num > 0) # zeroes are because these codons weren't genotyped

if(sum(master_table_clean$total_num < master_table_clean$mutant_num) > 0) {
  errorCondition("Entries with more mutant samples than total samples")
}
if(sum(master_table_clean$total_num == 0) > 0) {
  errorCondition("Entries with zero total samples")
}

# rename the columns so that they make sense
column_names <- get_column_names_for_clean()

# first add in extra information needed within column_names - !IMPORTANT NOTE! we are losing 12 studies here... may not be appropriate to dplyr::filter(grepl("calculated", substudy) == 0)
master_table_formatted <- master_table_clean %>% 
  group_by(study_uid) %>% 
  dplyr::mutate(study_uid = janitor::make_clean_names(study_uid[1], "big_camel")) %>% 
  ungroup() %>% 
  dplyr::group_by(study_uid, site_name, collection_start) %>% 
  dplyr::mutate(survey_ID = paste0("geoff_", study_uid[1], "_", 
                                   janitor::make_clean_names(site_name[1], "big_camel"), "_", 
                                   lubridate::year(collection_start[1])),
                continent = countrycode::countrycode(sourcevar = iso3c[1], "iso3c", "continent")) %>% 
  ungroup() %>% 
  dplyr::mutate(database = "GEOFF",
                study_nm = study_uid,
                prev = mutant_num/total_num,
                spatial_notes = "Literature review sourced lat/lon",
                time_notes = "Automated midpoint"
                )  %>% 
  # second check we have all the correctly named columns 
  dplyr::rename(study_ID = "study_uid",
                study_name = "study_nm",
                study_type = "publication_status",
                authors = "first_author_surname", 
                publication_year = "publication_year",
                url = "study_url",
                survey_ID = "survey_ID",
                country_name = "country", 
                site_name = "site_name",
                lat = "lat_n",
                lon = "lon_e",
                spatial_notes = "spatial_notes",
                collection_start = "collection_start", 
                collection_end = "collection_end", 
                collection_day = "collection_day", 
                time_notes = "time_notes",
                variant_string = "gene_mutation",
                variant_num = "mutant_num",
                total_num = "total_num",
                iso3c = "iso3c", 
                continent = "continent",
                pmid = "pmid",
                prev = "prev",
                gene = "gene",
                database = "database"
  ) %>% 
  # filter out any entries that are calculated - these reflect extra data entries for total genotype counts
  # which are already captured within haplotype counts
  dplyr::filter(grepl("calculated", substudy) == 0)

################################################################################
# Impute k13 reference survey records not explicitly entered in data entry
################################################################################

# Extract K13 mutations from the mutation dictionary
k13_mutations <- read_csv(here("analysis", "data-raw", "mutation_dictionary.csv")) %>%
  filter(gene == "k13") %>%
  mutate(mut = sub("^[A-Z]", "", mut)) %>%  # Remove the ref(first) codon letter (e.g., P553L → 553L)
  select(mut)

# Define K13 reference range
k13_range <- 1:726

# for each unique study_uid, impute k13 validated and candidate reference call counts using mean k13 sample coverage per study and k13 sequenced range before adding to add to master_table_formatted

# Get unique survey IDs
unique_surveys <- unique(master_table_formatted$survey_ID)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(unique_surveys), style = 3)  # Style 3 = moving bar

# Run the processing sequentially using lapply and show progress
imputed_records <- lapply(unique_surveys, process_survey)

# Close progress bar
close(pb)

imputed_data <- bind_rows(imputed_records)

# Ensure no duplicate variant_string entries for the same survey_ID
duplicates <- inner_join(imputed_data, master_table_formatted, 
                         by = c("survey_ID", "variant_string"))

if (nrow(duplicates) > 0) {
  message("Warning: ", nrow(duplicates), " duplicate variant_string entries found and removed.")
  imputed_data <- anti_join(imputed_data, master_table_formatted, 
                            by = c("survey_ID", "variant_string"))
}

# Append imputed records to master_table_formatted
master_table_formatted <- bind_rows(master_table_formatted, imputed_data)

# TODO: add check to confirm all k13_alleles (24 validated and candidate) are present for each survey (year/location)

# remove site="none" records from study_ID==S0006WrairKenreadKen which were generated automatically from MIP pipeline reading metadata for samples without site info which have no lat/lon. To be fixed manually in entry eventually.
master_table_formatted <- master_table_formatted[!(master_table_formatted$study_ID == "S0006WrairKenreadKen" & 
                                                     is.na(master_table_formatted$lat) & 
                                                     is.na(master_table_formatted$lon)), ]

# flag any study which has NA for lat or long after format conversion
lat_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lat)])
lon_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lon)])
# flag any study which has NA for substudy which is selectively used for inclusion of data in next scripts (ex we must use this field to exclude categories like tesday0treatmentfailure or treatedextracted for population prevalence calcs)
substudy_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$substudy)])
cat("Studies with missing or improper latitudes:\n", paste(lat_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper longitudes:\n", paste(lon_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper substudy values (critical for inclusion/exclusion in downstream analysis):\n", 
    paste(substudy_missing_or_improper, collapse = ", "), "\n\n")

# Grab just the columns we need for pairing with WWARN etc
master_table_simplified <- master_table_formatted %>% 
  select(all_of(column_names))

# flag all study uids which report less than 40% prevalence of validated k13 reference allele indicating entry error (using 469C arbitratily to detect erroneous user behavior) 
low_prev_studies <- master_table_formatted %>%
  filter(variant_string == "k13:469:C", prev < 0.40) %>%
  distinct(study_ID, data_entry_author)  # Get unique pairs

print(low_prev_studies)

# Save the final merged_df as an RDS file
saveRDS(master_table_simplified, here("analysis", "data-derived", "geoff_clean.rds"))
saveRDS(master_table_formatted, here("analysis", "data-derived", "geoff_clean_complete.rds"))
