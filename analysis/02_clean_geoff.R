# Load required libraries
library(tidyverse)
library(lubridate)
library(here)
library(purrr)
library(stringr)
devtools::load_all()

################################################################################
#
# Script: 02_clean_geoff.R
# Purpose: Perform various data cleaning to allow deduplication and combining
# with other data sources within {stave} object
# 
################################################################################

# 1. Step 1 - Read In Data and do broad cleaning # -----------------------------

# read in data
# Load the combined geoff data table created in the first script
master_table <- readRDS(here("analysis", "data-derived", 
                             "geoff_res.rds"))

# Load mutation key for k13 reference ranges
mutation_key <- read.csv(here("analysis", "data-raw", 
                              "k13_ref_protein_codon_dictionary.csv"))

## start by looking at the data - explore this in the console
master_unique <- lapply(master_table, unique)

# Initial main data cleaning
master_table_clean <- master_table |>
  dplyr::mutate(
    # Clean substudy field
    substudy = gsub(" ","", substudy),
    substudy = gsub("_","", substudy),
    substudy = if_else(substudy == "untreadextracted", "untreatedextracted", substudy),
    substudy = if_else(substudy == "extracted", "untreatedextracted", substudy),
    substudy = if_else((substudy == "day0" & study_uid == "s0029_warsame_2019"), "day0extracted", substudy),
    
    # Fix coordinates for specific study
    lat_n = if_else(study_uid == "s0020_some_2024", "11.79479482", lat_n),
    lon_e = if_else(study_uid == "s0020_some_2024", "-2.919282118", lon_e),
    
    # Clean site and location names
    site_name = gsub(" ", "_", tolower(site_name)),
    collection_location = gsub(" ", "_", tolower(collection_location)),
    collection_location = gsub(",", "", tolower(collection_location)),
    
    # Standardize country names
    country = stringr::str_to_title(country),
    country = if_else(country %in% c("Democratic Republic Of The Congo", "Drc"), 
                     "Democratic Republic of the Congo", country),
    country = if_else(country == "Rukara", "Rwanda", country),
    
    # Clean pretreatment samples
    pretreatment_samples = if_else(pretreatment_samples == "", NA, pretreatment_samples),
    
    # Standardize age fields
    study_design_age_min_years = if_else(study_design_age_min_years == "na", NA, study_design_age_min_years),
    study_design_age_min_years = if_else(study_design_age_min_years == "1year", "1", study_design_age_min_years),
    study_design_age_max_years = if_else(study_design_age_max_years == "na", NA, study_design_age_max_years),
    
    # Clean site study type
    site_study_type = if_else(study_uid == "s0034_uwimana", "TES", site_study_type),
    site_study_type = toupper(gsub("_","", site_study_type)),
    
    # Standardize publication status
    publication_status = gsub(" ", "_", publication_status),
    publication_status = if_else(publication_status == "na", NA, publication_status),
    publication_status = if_else(str_detect(study_uid, "unpub"), "unpublished", publication_status),
    publication_status = if_else(publication_status == "draft_paper", "unpublished", publication_status),
    publication_status = if_else(study_uid %in% c("s0031_bergmann_2021", "s0057_koko", "s0059_osborne",
                                                 "s0006_wrair_kenread_ken", "s0011_jacques-mari_gmms_unpub", 
                                                 "s0014_young_rwa_unpub"),
                                "peer_reviewed", publication_status),
    publication_status = if_else(publication_status == "published" & pmid == "37420265", "peer_reviewed", publication_status),
    publication_year = if_else(publication_status == "unpublished", NA, publication_year),
    
    # Clean countries covered
    countries_covered = gsub(",", ", ", countries_covered),
    countries_covered = stringr::str_to_title(countries_covered),
    countries_covered = if_else(countries_covered %in% c("Democratic Republic Of The Congo", "Drc"),
                               "Democratic Republic of the Congo", country),
    
    # Standardize data processing pipeline
    data_processing_pipeline = toupper(data_processing_pipeline),
    
    # Standardize author names
    first_author_surname = stringr::str_to_title(first_author_surname),
    
    # Fix ISO3C codes
    iso3c = replace(iso3c, iso3c == "DRC", "COD"),
    iso3c = replace(iso3c, iso3c == "ERT", "ERI"),
    iso3c = replace(iso3c, iso3c == "GQN", "GNQ"),
    country = replace(country, country == "Gqn", "Gnq"),
    
    # Fix missing site study types
    site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "32459360", "COMMUNITY"),
    site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "37670357", "COMMUNITY")
  )

# Additional country and data cleaning
master_table_clean <- master_table_clean |>
  dplyr::mutate(
    # Fix malformed iso3c codes used as country names
    country = if_else(country %in% c("Ago", "Ben", "Cmr", "Caf", "Cog", "Cod", "Eth", "Gha", 
                                    "Stp", "Sen", "Sdn", "Uga", "Gnq", "Gab", "Mli", 
                                    "Nga", "Sle", "Som", "Zmb", "Zaf"),
                     suppressWarnings(countrycode::countrycode(toupper(country), "iso3c", "country.name.en")),
                     country),
    
    # Fix pretreatment samples
    pretreatment_samples = replace(pretreatment_samples, 
                                  pretreatment_samples %in% c("pre_treatment", "Y"), "yes"),
    
    # Fix publication year
    publication_year = replace(publication_year, publication_year %in% c("none"), NA)
  )

# Data validation tests
african_countries <- data.frame(
  country = countrycode::codelist$country.name.en,
  continent = countrycode::codelist$continent,
  iso3c = countrycode::codelist$iso3c
) |>
  dplyr::filter(continent == "Africa") |> 
  dplyr::distinct()

# Define allowed values for validation
allowed_countries <- c(african_countries$country, "Democratic Republic of the Congo")
allowed_iso3c <- c(african_countries$iso3c)
allowed_substudy <- c("untreatedextracted", "untreatedcalculated",
                     "treatedextracted", "treatedcalculated",
                     "day0extracted", "day0calculated", 
                     "tesday0treatmentfailure", "day0reinfectionextracted",
                     "day0recrudescenceextracted")
allowed_pretreatment <- c("no", "yes", NA)
allowed_publication <- c("peer_reviewed", "preprint", "unpublished", "data_incomplete")
allowed_pub_year <- c(as.character(seq(from = 2000, to = lubridate::year(Sys.Date()))), NA)
allowed_site_types <- c("HEALTHFACILITY", "COMMUNITY", "TES", "CCS", "DHS", 
                       "REFUGEERECEPTIONCENTER", "UNKNOWN", "SURVEILLANCE")

# Run validation tests
check_values_in_column(master_table_clean, "country", allowed_countries)
check_values_in_column(master_table_clean, "iso3c", allowed_iso3c)
check_values_in_column(master_table_clean, "substudy", allowed_substudy)
check_values_in_column(master_table_clean, "pretreatment_samples", allowed_pretreatment)
check_values_in_column(master_table_clean, "publication_status", allowed_publication)
check_values_in_column(master_table_clean, "publication_year", allowed_pub_year)
check_values_in_column(master_table_clean, "site_study_type", allowed_site_types)

# 2. Step 2 - Clean mutation names and data entry issues # --------------------

# Clean mutation names and expand K13 ranges
# Extract K13 mutations from the mutation dictionary
k13_mutations <- read_csv(here("analysis", "data-raw", "mutation_dictionary.csv")) |>
  dplyr::filter(gene == "k13") |>
  dplyr::mutate(mut = sub("^[A-Z]", "", mut)) |>  # Remove ref codon letter (e.g., P553L → 553L)
  dplyr::select(mut)

mutation_positions <- as.numeric(gsub("[^0-9]", "", k13_mutations$mut))

# Identify K13 range mutations that need expansion
indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))

# Expand K13 range mutations to individual positions
if (length(indices_to_transform) > 0) {
  rows_to_expand <- master_table_clean[indices_to_transform, ]
  
  # Expand each using mapply
  expanded_rows <- mapply(
    expand_k13_range_to_rows,
    gene_mutation = rows_to_expand$gene_mutation,
    row_data = split(rows_to_expand, seq_len(nrow(rows_to_expand))),
    MoreArgs = list(mutation_key = mutation_key, mutation_positions = mutation_positions),
    SIMPLIFY = FALSE
  )
  
  # Combine all into one data frame
  expanded_rows_df <- dplyr::bind_rows(expanded_rows)
  
  # Create final master table by removing collapsed ones and adding expanded
  master_table_clean <- dplyr::bind_rows(
    master_table_clean[-indices_to_transform, ],
    expanded_rows_df
  )
}

# Handle individual K13 reference positions with asterisk
single_ref_indices_to_transform <- which(grepl("^k13:[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)) &
                                         master_table_clean$mutant_num == 0)

if (length(single_ref_indices_to_transform) > 0) {
  master_table_clean$gene_mutation[single_ref_indices_to_transform] <- as.character(sapply(
    master_table_clean$gene_mutation[single_ref_indices_to_transform],
    convert_k13_asterisk, 
    mutation_key = mutation_key
  ))
  
  master_table_clean$mutant_num[single_ref_indices_to_transform] <- master_table_clean$total_num[single_ref_indices_to_transform]
}

# Remove remaining K13 asterisk entries (STOP codons that can't be handled)
indices_to_transform <- which(grepl("^k13:[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))
if (length(indices_to_transform) > 0) {
  master_table_clean <- master_table_clean[-indices_to_transform, ]
}

# Fix specific mutation entries
# Fix MDR1 FC insertion (convert to wildtype S)
mdr1_fc_indices <- which(master_table_clean$gene_mutation == "MDR1:1034:FC")
if (length(mdr1_fc_indices) > 0) {
  master_table_clean$mutant_num[mdr1_fc_indices] <- master_table_clean$total_num[mdr1_fc_indices]
  master_table_clean$gene_mutation[mdr1_fc_indices] <- "MDR1:1034:S"
}

# Remove K13 insertion mutations (not SNPs)
master_table_clean <- master_table_clean |> 
  dplyr::filter(!(gene_mutation %in% c("K13:142:NN", "K13:142:NNN")))

# Correct gene mutation format for STAVE compatibility
convmutation <- format_variants_for_stave(unique(master_table_clean$gene_mutation))
master_table_clean$gene_mutation <- convmutation[match(master_table_clean$gene_mutation, unique(master_table_clean$gene_mutation))]

# 3. Step 3 - Clean dates, coordinates, and data types # ----------------------

# Adjust and standardize dates
convstart <- adjust_invalid_date(unique(master_table_clean$date_start), is_start = TRUE)
master_table_clean$collection_start <- convstart[match(master_table_clean$date_start, unique(master_table_clean$date_start))]

convend <- adjust_invalid_date(unique(master_table_clean$date_end), is_start = FALSE)
master_table_clean$collection_end <- convend[match(master_table_clean$date_end, unique(master_table_clean$date_end))]

# Compute collection day as median of start and end dates
master_table_clean <- master_table_clean |> 
  dplyr::rowwise() |> 
  dplyr::mutate(collection_day = median(c(collection_start, collection_end))) |> 
  dplyr::ungroup()

# Fix coordinates that may have been read as Excel formulas
master_table_clean$lat_n <- clean_excel_formulas(master_table_clean$lat_n)
master_table_clean$lon_e <- clean_excel_formulas(master_table_clean$lon_e)

# Convert variable classes and extract gene information
master_table_clean <- master_table_clean |>
  dplyr::mutate(
    collection_day = as.Date(collection_day),
    collection_start = as.Date(collection_start),
    collection_end = as.Date(collection_end),
    lat_n = as.numeric(lat_n),
    lon_e = as.numeric(lon_e),
    mutant_num = as.numeric(mutant_num),
    total_num = as.numeric(total_num),
    publication_year = as.numeric(publication_year),
    gene = str_extract(gene_mutation, "^[^:]+"),
    mut = str_extract(gene_mutation, "(?<=:).*")
  ) |>
  dplyr::filter(total_num > 0)  # Remove entries where codons weren't genotyped

# 4. Step 4 - Fix specific data entry issues # --------------------------------

# Fix specific data entry error in Koulamoutou study
master_table_clean <- master_table_clean |> 
  dplyr::mutate(mutant_num = replace(mutant_num, 
                                    mutant_num == 31 & total_num == 24 & site_uid == "koulamoutou", 21))

# Remove records without coordinates from specific study
master_table_clean <- master_table_clean |>
  dplyr::filter(!(study_uid == "s0006_wrair_kenread_ken" & 
                  is.na(lat_n) & 
                  is.na(lon_e)))

# Remove day0calculated entries (duplicate genotype counts already captured in haplotypes)
master_table_clean <- master_table_clean |> 
  dplyr::filter(grepl("day0calculated", substudy) == 0)


# 5. Step 5 - Data validation and formatting # --------------------------------

# Data integrity checks
if (sum(master_table_clean$total_num < master_table_clean$mutant_num) > 0) {
  errorCondition("Entries with more mutant samples than total samples")
}
if (sum(master_table_clean$total_num == 0) > 0) {
  errorCondition("Entries with zero total samples")
}

# Get column names for final output
column_names <- get_column_names_for_clean()

# Format data for STAVE compatibility
master_table_formatted <- master_table_clean |> 
  dplyr::group_by(study_uid) |> 
  dplyr::mutate(study_uid = paste0("geoff_", janitor::make_clean_names(study_uid[1], "big_camel"))) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(study_uid, site_name, collection_start) |> 
  dplyr::mutate(
    survey_ID = paste0(study_uid[1], "_", 
                      janitor::make_clean_names(site_name[1], "big_camel"), "_", 
                      lubridate::year(collection_start[1])),
    continent = countrycode::countrycode(sourcevar = iso3c[1], "iso3c", "continent")
  ) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(
    database = "GEOFF",
    study_nm = study_uid,
    prev = mutant_num/total_num,
    spatial_notes = "Literature review sourced lat/lon",
    time_notes = "Automated midpoint"
  ) |> 
  # Rename columns for consistency
  dplyr::rename(
    study_ID = "study_uid",
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
  )

# Replace missing URLs with placeholder to satisfy downstream URL checkers
studies_missing_url <- master_table_formatted |>
  dplyr::filter(is.na(url)) |>
  dplyr::distinct(study_ID) |>
  dplyr::pull(study_ID)

if (length(studies_missing_url) > 0) {
  cat("Replacing missing URLs with placeholder for", length(studies_missing_url), "studies:\n")
  for (study in studies_missing_url) {
    cat("- ", study, "\n")
  }
  
  # Replace missing URLs with placeholder
  master_table_formatted <- master_table_formatted |>
    dplyr::mutate(url = if_else(is.na(url), "https://placeholder-url-missing.com", url))
  
  cat("Placeholder URL used: https://placeholder-url-missing.com\n\n")
} else {
  cat("All studies have URLs - no placeholder replacement needed.\n\n")
}

# 6. Step 6 - Impute K13 reference records # ----------------------------------

# Define K13 reference range
k13_range <- 1:726

# Impute K13 validated and candidate reference call counts for each survey
unique_surveys <- unique(master_table_formatted$survey_ID)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(unique_surveys), style = 3)

# Process surveys and impute missing K13 records
imputed_records <- lapply(unique_surveys, process_survey)

# Close progress bar
close(pb)

# Combine imputed records
imputed_data <- bind_rows(imputed_records)

# Check for and remove duplicate variant_string entries
duplicates <- inner_join(imputed_data, master_table_formatted, 
                        by = c("survey_ID", "variant_string"))

if (nrow(duplicates) > 0) {
  message("Warning: ", nrow(duplicates), " duplicate variant_string entries found and removed.")
  imputed_data <- anti_join(imputed_data, master_table_formatted, 
                           by = c("survey_ID", "variant_string"))
}

# Append imputed records to main table
master_table_formatted <- bind_rows(master_table_formatted, imputed_data)

# Validate K13 allele completeness (should return 0 rows if all alleles present)
k13_validation <- master_table_formatted |> 
  dplyr::filter(gene == "k13") |> 
  dplyr::group_by(survey_ID) |> 
  dplyr::summarise(
    n = str_extract_all(variant_string, "(?<=k13:)([\\d_]+)") |>
        unlist() |>
        str_split("_") |>
        unlist() |>
        as.numeric() |>
        unique() |> 
        length(),
    min = as.integer(unique(k13_min)), 
    max = as.integer(unique(k13_max)), 
    pos = sum(as.integer(unique(substr(k13_mutations$mut, 1, 3))) < max & 
              as.integer(unique(substr(k13_mutations$mut, 1, 3))) > min),
    .groups = "drop"
  ) |> 
  dplyr::filter(n < pos)

# 7. Step 7 - Data quality checks # -------------------------------------------

# Identify studies with missing coordinates
lat_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lat)])
lon_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lon)])

# Identify studies with missing substudy values (critical for downstream filtering)
substudy_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$substudy)])

# Report missing data
cat("Studies with missing or improper latitudes:\n", paste(lat_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper longitudes:\n", paste(lon_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper substudy values (critical for inclusion/exclusion in downstream analysis):\n", 
    paste(substudy_missing_or_improper, collapse = ", "), "\n\n")

# Flag studies with low K13 reference allele prevalence (< 40% indicates potential entry error)
low_prev_studies <- master_table_formatted |>
  dplyr::filter(variant_string == "k13:469:C", prev < 0.40) |>
  dplyr::distinct(study_ID, data_entry_author)
print(low_prev_studies)

# 8. Step 8 - Any further specific data entry issues, likely identified when trying to read to STAVE ---------

# Error - Now fixed: number of amino acid loci (5) must equal the number of codon positions (3)
# This has now been fixed on 02_clean_geoff.R
# e.g. crt:72_73_74_75_76:CVMNK/CVIET/SVMNK which should be crt:72_73_74_75_76:C/S_V_M/I_N/E_K/T

# last fix for some incorrectly formatted loci
issue <- c(
  "mdr1:86_184_1246:YFD/YYD",
  "mdr1:86_184_1246:YYD/NFD",
  "mdr1:86_184_1246:NFD/YYD",
  "mdr1:86_184_1246:YYY/NYD",
  "mdr1:86_184_1246:NFY/NYD",
  "mdr1:86_184_1246:NYD/NFD",
  "mdr1:86_184_1246:NFD/NYD",
  "mdr1:86_184_1246:YFD/YYD",
  "mdr1:86_184_1246:YYD/NFD",
  "mdr1:86_184_1246:NFD/YYD",
  "mdr1:86_184_1246:YYY/NYD",
  "mdr1:86_184_1246:NFY/NYD",
  "mdr1:86_184_1246:NYD/NFD",
  "mdr1:86_184_1246:NFD/NYD",
  "crt:72_73_74_75_76:CVMNK/CVIET",
  "crt:72_73_74_75_76:CVMNK/SVMNK",
  "crt:72_73_74_75_76:CVMNK/CVIET/SVMNK"
)

# Create a helper function for readability
fix_mutation_string <- function(issue_str) {
  parts <- str_split(issue_str, ":")[[1]]
  loci <- parts[1:2] %>% paste(collapse = ":")
  
  alleles <- parts[3] %>%
    str_split("/") %>%
    unlist() %>%
    str_split("") %>%
    simplify2array()
  
  allele_summary <- apply(alleles, 1, function(pos) paste(unique(pos), collapse = "/")) %>%
    paste(collapse = "_")
  
  paste(loci, allele_summary, sep = ":")
}

# Use map_chr for clarity
fixed <- map_chr(issue, fix_mutation_string)

# Simple indexing to update master_table_formatted
indices <- match(master_table_formatted$variant_string, issue)
master_table_formatted$variant_string[!is.na(indices)] <- fixed[indices[!is.na(indices)]]

# Then two other issues
# "k13:578_579_580:BND" "crt:76:TK" 
#https://malariajournal.biomedcentral.com/articles/10.1186/s12936-021-03713-2#Sec7
master_table_formatted$variant_num[master_table_formatted$variant_string == "k13:578_579_580:BND"] <- 
master_table_formatted$total_num[master_table_formatted$variant_string == "k13:578_579_580:BND"]
master_table_formatted$variant_string[master_table_formatted$variant_string == "k13:578_579_580:BND"] <- "k13:578_579_580:A_M_C"
#https://malariajournal.biomedcentral.com/articles/10.1186/s12936-017-1777-0/tables/2
master_table_formatted$variant_string[master_table_formatted$variant_string == "crt:76:TK"] <- "crt:76:T/K"

# 9. Step 9 - Placeholder for final column selection # -----------------------
# Note: master_table_simplified will be created after all cleaning is complete

# 10. Step 10 - Report missing studies # --------------------------------------

# Load expected study UIDs and compare with cleaned data
all_expected_sample_uids <- read.table(
  here("analysis", "data-raw", "entered_geoff_study_ids.txt"),
  stringsAsFactors = FALSE
)[[1]]

# Helper function to capitalize underscore-separated parts
capitalize_underscore_parts <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  parts_cap <- paste0(toupper(substring(parts, 1, 1)), substring(parts, 2))
  paste0(parts_cap, collapse = "")
}

# Clean and normalize expected study UIDs
all_expected_sample_uids_cleaned <- all_expected_sample_uids |>
  lapply(capitalize_underscore_parts) |>
  unlist(use.names = FALSE) |>
  (\(x) gsub(" ", "", x))() |>
  (\(x) paste0("geoff_", x))() |>
  as.character()

# Find study IDs in cleaned data that are not in expected list
unexpected_study_ids <- master_table_formatted |>
  dplyr::distinct(study_ID) |>
  dplyr::filter(!study_ID %in% all_expected_sample_uids_cleaned)

# Find expected study IDs that are missing from cleaned data
missing_study_ids <- data.frame(study_ID = all_expected_sample_uids_cleaned) |>
  dplyr::filter(!study_ID %in% unique(master_table_formatted$study_ID))

# Extract study numbers for flexible matching (e.g., S0148 from geoff_S0148Martínez-pérez2018)
extract_study_number <- function(study_id) {
  stringr::str_extract(study_id, "(?<=geoff_)(S[0-9]+)(?=[A-Za-z])")
}

# Get study numbers from both datasets
expected_study_numbers <- sapply(all_expected_sample_uids_cleaned, extract_study_number)
actual_study_numbers <- sapply(unique(master_table_formatted$study_ID), extract_study_number)

# Find missing studies by study number (more flexible matching)
missing_by_number <- data.frame(
  study_ID = all_expected_sample_uids_cleaned,
  study_number = expected_study_numbers
) |>
  dplyr::filter(!is.na(study_number) & !study_number %in% actual_study_numbers)

# Report some final study inclusion checks. Using the more flexible check for now as some study ids will be different in the cleaned data than in the entry sheet.
# Report unexpected study IDs
#cat("The following study_IDs are present in master_table_formatted but missing from entered_geoff_study_ids.txt (after normalization):\n")
#print(unexpected_study_ids)

# Report missing study IDs (exact match)
#cat("\nThe following study_IDs are in entered_geoff_study_ids.txt but missing from master_table_formatted (exact match):\n")
#print(missing_study_ids)

# Report missing study IDs (flexible study number match)
cat("\nThe following study numbers are in entered_geoff_study_ids.txt but missing from master_table_formatted (flexible match by study number):\n")
print(missing_by_number)

# 11. Step 11 - Final validation checks # -------------------------------------

# Check for missing variant numbers
na_variant_num <- master_table_formatted |>
  dplyr::filter(is.na(variant_num)) |>
  dplyr::distinct(survey_ID) |>
  dplyr::mutate(problem_column = "variant_num")

# Check for missing total numbers
na_total_num <- master_table_formatted |>
  dplyr::filter(is.na(total_num)) |>
  dplyr::distinct(survey_ID) |>
  dplyr::mutate(problem_column = "total_num")

# Check for missing URLs
na_url <- master_table_formatted |>
  dplyr::filter(is.na(url)) |>
  dplyr::distinct(survey_ID) |>
  dplyr::mutate(problem_column = "url")

# Combine and report all missing value findings
problem_surveys <- bind_rows(na_variant_num, na_total_num, na_url) |>
  dplyr::arrange(survey_ID, problem_column)

if (nrow(problem_surveys) > 0) {
  message("Surveys with missing values:")
  print(problem_surveys)
} else {
  message("No NAs found in variant_num, total_num, or url.")
}



# Step 11.5 - Update URLs for specific studies # ------------------------------

# Define URL updates for studies with missing or incorrect URLs
url_updates <- c(
  "geoff_S0002AyelawEth2023"   = "https://pubmed.ncbi.nlm.nih.gov/40666313/",
  "geoff_S0006WrairKenreadKen" = "https://www.medrxiv.org/content/10.1101/2025.07.15.25331603v1",
  "geoff_S0014YoungRwaUnpub"   = "https://academic.oup.com/jid/article/231/1/269/7811784",
  "geoff_S0007Connelly2024Zim" = "https://www.medrxiv.org/content/10.1101/2025.03.22.25323829v1"
)

# Update URLs in the formatted table
master_table_formatted <- master_table_formatted |>
  dplyr::mutate(url = if_else(
    study_ID %in% names(url_updates),
    url_updates[study_ID],
    url
  ))

# Step 11.6 - Additional data cleaning patches # ------------------------------

# 1) Fix variant string ordering for study S0046 (geoff_S0046Bwire)
# Check before fix
s0046_before <- master_table_formatted |>
  dplyr::filter(study_ID == "geoff_S0046Bwire" & variant_string == "k13:578_565:S/C") |>
  nrow()

master_table_formatted <- master_table_formatted |>
  dplyr::mutate(variant_string = if_else(
    study_ID == "geoff_S0046Bwire" & variant_string == "k13:578_565:S/C",
    "k13:565_578:C_S",
    variant_string
  ))

# Check after fix
s0046_after <- master_table_formatted |>
  dplyr::filter(study_ID == "geoff_S0046Bwire" & variant_string == "k13:565_578:C_S") |>
  nrow()

# Report the fix
cat("S0046 variant string fix - Records changed from 'k13:578_565:S/C' to 'k13:565_578:C_S':", s0046_before, "->", s0046_after, "\n")

# 2) Check codon and base pair counts match in variant strings
codon_bp_mismatch <- master_table_formatted |>
  dplyr::mutate(
    # Extract codon positions (between first : and second :)
    codon_part = stringr::str_extract(variant_string, "(?<=:)[^:]+(?=:)"),
    # Extract allele part (after second :)
    allele_part = stringr::str_extract(variant_string, "(?<=:)[^:]+$"),
    # Count codons (split by _ and count elements)
    codon_count = lengths(stringr::str_split(codon_part, "_")),
    # Count alleles (split by _ and count elements)
    allele_count = lengths(stringr::str_split(allele_part, "_"))
  ) |>
  dplyr::filter(codon_count != allele_count) |>
  dplyr::select(study_ID, variant_string, codon_count, allele_count)

# Report mismatches
if (nrow(codon_bp_mismatch) > 0) {
  cat("\nVariant strings with codon/allele count mismatches:\n")
  print(codon_bp_mismatch)
} else {
  cat("\nAll variant strings have matching codon and allele counts.\n")
}

# 3) Fix variant strings by adding underscores between allele letters
if (nrow(codon_bp_mismatch) > 0) {
  # Function to add underscores between individual letters in allele part
  fix_allele_format <- function(variant_string, codon_count, allele_count) {
    # Only fix if we have more codons than alleles and alleles appear to be concatenated letters
    if (codon_count > allele_count) {
      parts <- stringr::str_split(variant_string, ":")[[1]]
      if (length(parts) == 3) {
        gene <- parts[1]
        codons <- parts[2]
        alleles <- parts[3]
        
        # If alleles are letters without separators, add underscores
        if (stringr::str_detect(alleles, "^[A-Z]+$") && nchar(alleles) == codon_count) {
          # Split into individual letters and join with underscores
          allele_letters <- stringr::str_split(alleles, "")[[1]]
          fixed_alleles <- paste(allele_letters, collapse = "_")
          return(paste(gene, codons, fixed_alleles, sep = ":"))
        }
      }
    }
    return(variant_string)
  }
  
  # Apply the fix to mismatched variant strings
  master_table_formatted <- master_table_formatted |>
    dplyr::mutate(
      # Extract codon and allele counts again for the fix
      codon_part_temp = stringr::str_extract(variant_string, "(?<=:)[^:]+(?=:)"),
      allele_part_temp = stringr::str_extract(variant_string, "(?<=:)[^:]+$"),
      codon_count_temp = lengths(stringr::str_split(codon_part_temp, "_")),
      allele_count_temp = lengths(stringr::str_split(allele_part_temp, "_")),
      
      # Fix variant strings where needed
      variant_string = mapply(fix_allele_format, variant_string, codon_count_temp, allele_count_temp)
    ) |>
    dplyr::select(-codon_part_temp, -allele_part_temp, -codon_count_temp, -allele_count_temp)
  
  # Re-check after fix
  codon_bp_mismatch_after <- master_table_formatted |>
    dplyr::mutate(
      codon_part = stringr::str_extract(variant_string, "(?<=:)[^:]+(?=:)"),
      allele_part = stringr::str_extract(variant_string, "(?<=:)[^:]+$"),
      codon_count = lengths(stringr::str_split(codon_part, "_")),
      allele_count = lengths(stringr::str_split(allele_part, "_"))
    ) |>
    dplyr::filter(codon_count != allele_count) |>
    dplyr::select(study_ID, variant_string, codon_count, allele_count)
  
  cat("\nAfter fixing allele formatting - remaining mismatches:\n")
  if (nrow(codon_bp_mismatch_after) > 0) {
    print(codon_bp_mismatch_after)
  } else {
    cat("All variant strings now have matching codon and allele counts.\n")
  }
}

# 4) Fix specific variant strings for study geoff_S0122Li2015
master_table_formatted <- master_table_formatted |>
  dplyr::mutate(variant_string = case_when(
    study_ID == "geoff_S0122Li2015" & variant_string == "crt:72_73_74_75_76:CVM/IN/EK/T" ~ "crt:72_73_74_75_76:C_V_M/I_N/E_K/T",
    study_ID == "geoff_S0122Li2015" & variant_string == "mdr1:86_184:NY/F" ~ "mdr1:86_184:N_Y/F",
    study_ID == "geoff_S0122Li2015" & variant_string == "mdr1:86_184:YY/F" ~ "mdr1:86_184:Y_Y/F",
    TRUE ~ variant_string
  ))

# Final check after all fixes
final_codon_bp_check <- master_table_formatted |>
  dplyr::mutate(
    codon_part = stringr::str_extract(variant_string, "(?<=:)[^:]+(?=:)"),
    allele_part = stringr::str_extract(variant_string, "(?<=:)[^:]+$"),
    codon_count = lengths(stringr::str_split(codon_part, "_")),
    allele_count = lengths(stringr::str_split(allele_part, "_"))
  ) |>
  dplyr::filter(codon_count != allele_count) |>
  dplyr::select(study_ID, variant_string, codon_count, allele_count)

cat("\nFinal check - variant strings with codon/allele count mismatches:\n")
if (nrow(final_codon_bp_check) > 0) {
  print(final_codon_bp_check)
} else {
  cat("All variant strings now have matching codon and allele counts.\n")
}

# 5) Check for NA values in the URL column
na_url_check <- master_table_formatted |>
  dplyr::filter(is.na(url)) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

cat("\nStudies with missing URLs:\n")
if (nrow(na_url_check) > 0) {
  print(na_url_check)
  cat("Number of studies with missing URLs:", nrow(na_url_check), "\n")
} else {
  cat("All studies have URLs.\n")
}

# 6) Check for missing date fields
missing_dates_check <- master_table_formatted |>
  dplyr::summarise(
    missing_collection_day = sum(is.na(collection_day)),
    missing_collection_start = sum(is.na(collection_start)),
    missing_collection_end = sum(is.na(collection_end)),
    .groups = "drop"
  )

cat("\nMissing date field counts:\n")
print(missing_dates_check)

# Identify studies with missing dates
studies_missing_dates <- master_table_formatted |>
  dplyr::filter(is.na(collection_day) | is.na(collection_start) | is.na(collection_end)) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

if (nrow(studies_missing_dates) > 0) {
  cat("\nStudies with missing date fields:\n")
  print(studies_missing_dates)
  cat("Number of studies with missing dates:", nrow(studies_missing_dates), "\n")
  
  # Detailed examination of studies with missing collection days
  studies_missing_collection_day <- master_table_formatted |>
    dplyr::filter(is.na(collection_day)) |>
    dplyr::distinct(study_ID) |>
    dplyr::pull(study_ID)
  
  if (length(studies_missing_collection_day) > 0) {
    cat("\nDETAILED EXAMINATION - Studies with Missing Collection Days:\n")
    cat("Studies with missing collection_day:", length(studies_missing_collection_day), "\n")
    
    # Remove rows with missing collection_day from specific vetted studies
    # These are remnants from MIP processing pipeline that can be safely omitted
    vetted_studies_to_clean <- c("geoff_S0007Connelly2024Zim", "geoff_S0006WrairKenreadKen")
    
    rows_before <- nrow(master_table_formatted)
    master_table_formatted <- master_table_formatted |>
      dplyr::filter(!(study_ID %in% vetted_studies_to_clean & is.na(collection_day)))
    rows_after <- nrow(master_table_formatted)
    
    removed_count <- rows_before - rows_after
    cat("Removed", removed_count, "rows with missing collection_day from vetted studies:\n")
    cat("- geoff_S0007Connelly2024Zim\n")
    cat("- geoff_S0006WrairKenreadKen\n\n")
    
    # Summary of date field completeness for problematic studies
    date_completeness <- master_table_formatted |>
      dplyr::filter(study_ID %in% studies_missing_collection_day) |>
      dplyr::group_by(study_ID) |>
      dplyr::summarise(
        total_records = n(),
        missing_collection_day = sum(is.na(collection_day)),
        missing_collection_start = sum(is.na(collection_start)),
        missing_collection_end = sum(is.na(collection_end)),
        missing_date_start = sum(is.na(date_start)),
        missing_date_end = sum(is.na(date_end)),
        .groups = "drop"
      )
    
    cat("Date field completeness summary:\n")
    print(date_completeness)
  }
} else {
  cat("All studies have complete date information.\n")
}

# 7) Check for missing geographic coordinates
missing_coords_check <- master_table_formatted |>
  dplyr::summarise(
    missing_lat = sum(is.na(lat)),
    missing_lon = sum(is.na(lon)),
    missing_both = sum(is.na(lat) & is.na(lon)),
    .groups = "drop"
  )

cat("\nMissing coordinate field counts:\n")
print(missing_coords_check)

# Identify studies with missing coordinates
studies_missing_coords <- master_table_formatted |>
  dplyr::filter(is.na(lat) | is.na(lon)) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

if (nrow(studies_missing_coords) > 0) {
  cat("\nStudies with missing geographic coordinates:\n")
  print(studies_missing_coords)
  cat("Number of studies with missing coordinates:", nrow(studies_missing_coords), "\n")
} else {
  cat("All studies have complete geographic coordinates.\n")
}

# 8) Check for NA values in variant_num and total_num
missing_counts_check <- master_table_formatted |>
  dplyr::summarise(
    missing_variant_num = sum(is.na(variant_num)),
    missing_total_num = sum(is.na(total_num)),
    missing_both_counts = sum(is.na(variant_num) & is.na(total_num)),
    .groups = "drop"
  )

cat("\nMissing count field summary:\n")
print(missing_counts_check)

# Identify studies with missing variant_num
studies_missing_variant_num <- master_table_formatted |>
  dplyr::filter(is.na(variant_num)) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

if (nrow(studies_missing_variant_num) > 0) {
  cat("\nStudies with missing variant_num:\n")
  print(studies_missing_variant_num)
  cat("Number of studies with missing variant_num:", nrow(studies_missing_variant_num), "\n")
} else {
  cat("All studies have variant_num values.\n")
}

# Identify studies with missing total_num
studies_missing_total_num <- master_table_formatted |>
  dplyr::filter(is.na(total_num)) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

if (nrow(studies_missing_total_num) > 0) {
  cat("\nStudies with missing total_num:\n")
  print(studies_missing_total_num)
  cat("Number of studies with missing total_num:", nrow(studies_missing_total_num), "\n")
} else {
  cat("All studies have total_num values.\n")
}

# Flag problematic studies that may need attention
problematic_studies <- c(
  "geoff_S0006WrairKenreadKen",
  "geoff_S0013MsmtTza21", 
  "geoff_S0016Connelly2024Elife",
  "geoff_S0026Verity",
  "geoff_S0062Vonwowern2024",
  "geoff_S0067Makenga2022",
  "geoff_S0108Natama2020",
  "geoff_S0114Lepiscopia2021",
  "geoff_S0128VoumboMatoumona2018",
  "geoff_S0148MartinezPerez2018"
)

# Check which of these problematic studies have missing count data
problematic_with_missing_counts <- master_table_formatted |>
  dplyr::filter(study_ID %in% problematic_studies & (is.na(variant_num) | is.na(total_num))) |>
  dplyr::distinct(study_ID) |>
  dplyr::arrange(study_ID)

cat("\nProblematic studies (from expected list) with missing count data:\n")
if (nrow(problematic_with_missing_counts) > 0) {
  print(problematic_with_missing_counts)
} else {
  cat("None of the flagged problematic studies have missing count data.\n")
}

# 12. Step 12 - Final column selection and save cleaned data # ----------------

# Create simplified table with only required columns after all cleaning is complete
master_table_simplified <- master_table_formatted |> 
  dplyr::select(all_of(column_names))

# Comprehensive check for missing values in master_table_simplified
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("FINAL DATA QUALITY CHECK - Missing Values in master_table_simplified\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Count missing values by column
missing_summary <- master_table_simplified |>
  dplyr::summarise(dplyr::across(everything(), ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(everything(), names_to = "column", values_to = "missing_count") |>
  dplyr::filter(missing_count > 0) |>
  dplyr::arrange(desc(missing_count))

if (nrow(missing_summary) > 0) {
  cat("Columns with missing values:\n")
  print(missing_summary)
  cat("\nTotal records in dataset:", nrow(master_table_simplified), "\n")
  
  # Show percentage missing for each column
  missing_summary_pct <- missing_summary |>
    dplyr::mutate(percent_missing = round(missing_count / nrow(master_table_simplified) * 100, 2))
  
  cat("\nMissing values with percentages:\n")
  print(missing_summary_pct)
  
  # Flag columns with high missingness (>10%)
  high_missing <- missing_summary_pct |>
    dplyr::filter(percent_missing > 10)
  
  if (nrow(high_missing) > 0) {
    cat("\nWARNING: Columns with >10% missing data:\n")
    print(high_missing)
  }
} else {
  cat("✓ No missing values found in any column of master_table_simplified\n")
}

cat(paste(rep("=", 60), collapse = ""), "\n\n")

# 13. Imputing the partner drugs
# 13a. Starting with mdr1

# first split the dataframes so that I don't break the k13 which is working fine
gene_keep <- master_table_simplified %>% 
  filter(!(gene %in% c("mdr1", "crt")))

gene_simplified <- master_table_simplified %>% 
  filter(gene %in% c("mdr1", "crt")) %>%
  split(.$gene)

# need to filter out rows that don't contant mdr1 86 at all
gene_simplified$mdr1 <- gene_simplified$mdr1 %>%
  filter(grepl("^mdr1:[^:]*86[^:]*:", variant_string)) %>%
  dplyr::filter(variant_string != "mdr1:86:F")

# assign to uuids 
assign_ids <- function(x) {
  x %>%
    group_by(across(c(-variant_num, -total_num, -prev, -variant_string))) %>%
    mutate(uuid = cur_group_id()) %>%
    ungroup %>% 
    group_by(across(c(-variant_num, -prev, -variant_string))) %>%
    mutate(nid = cur_group_id()) %>%
    ungroup() %>% 
    mutate(iid = seq_len(n()))
}

# fix the sheer number of mutants we are dealing with -- combine all 86s together
# TODO: make this a proper function with documentatin and tests
standardise_mdr1_86 <- function(x) {
  # split into 3 parts: "mdr1", positions, alleles
  parts <- str_split_fixed(x, ":", 3)
  # alleles for *all positions*
  alleles_all <- parts[,3]
  # extract only the allele block for the first position (86)
  first_block <- str_extract(alleles_all, "^[A-Z/]+")
  # keep only A or A/B
  first_allele <- str_extract(first_block, "^[A-Z](?:/[A-Z])?")
  # ----- NEW RULE: drop F if mixed -----
  # N/F → N
  # Y/F → Y
  first_allele <- gsub("^N/F$", "N", first_allele)
  first_allele <- gsub("^Y/F$", "Y", first_allele)
  
  # normalise N/Y ordering
  first_allele <- ifelse(first_allele == "Y/N", "N/Y", first_allele)
  
  paste0("mdr1:86:", first_allele)
}

# TODO: test the function more thoroughly
# test this function
standardise_mdr1_86("mdr1:86:N")
# "mdr1:86:N"
standardise_mdr1_86("mdr1:86:Y/N")
# "mdr1:86:N/Y"
standardise_mdr1_86("mdr1:86_184_1246:N_Y_D")
# "mdr1:86:N"
standardise_mdr1_86("mdr1:86_184:Y/F_N")
# "mdr1:86:Y" 

# current approach:
# convert full variant strings into simplified version focused on codon 86
# group by study info and 86 variant, num = sum(variant_num), denom = unique(total_num) if there is only one unique value
# clean based on the length of unique(86_string)

# checked that the length of the denominator is always 1 when we group by nid
# mdr1 %>%
# mutate(string_86 = standardise_mdr1_86(variant_string)) %>% 
# filter(string_86 != "mdr1:86:F") %>%
# group_by(survey_ID, collection_day, nid, string_86) %>%
# mutate(denom = length(unique(total_num))) %>%
# filter(denom > 1) %>% View()
# mutate(x = sum(variant_num), n = unique(total_num))

# make the original easier to examine and add ids to the 
mdr1_original <- assign_ids(gene_simplified$mdr1) 

# some issues with different ns in different rows -- these need manual fixing
# find these studies: same survey ID and different total_num
# mdr1_fix -- these have been sent to Neeva to manually exmaine
mdr1_fix <- mdr1_original %>%
  dplyr::group_by(survey_ID, collection_day) %>%
  mutate(n = length(unique(total_num))) %>%
  filter(n > 1) %>%
  arrange(survey_ID, total_num) %>%
  split(.$study_name)

# fix these and add them back into the cleaning

# keep the bottom row -- impute the others
mdr1_fix$geoff_S0036Schreidah <- mdr1_fix$geoff_S0036Schreidah[3,]

# keep bottom 3 rows
mdr1_fix$geoff_S0062Vonwowern2024 <- mdr1_fix$geoff_S0062Vonwowern2024[4:6,]

# Plucinski needs new entries
mdr1_fix$geoff_S0101Plucinski2017 <- mdr1_fix$geoff_S0101Plucinski2017 %>%
  distinct(site_name, collection_day, total_num, .keep_all = TRUE) %>%
  mutate(variant_string = "mdr1:86:N") %>%
  ungroup() %>%
  mutate(variant_num = c(1, 13, 9, 22)) %>% 
  mutate(prev = variant_num / total_num)

# fix incorrect denominator
mdr1_fix$geoff_S0115Tuedomagb2021$total_num <- 739
mdr1_fix$geoff_S0115Tuedomagb2021$prev <- mdr1_fix$geoff_S0115Tuedomagb2021$variant_num / 
  mdr1_fix$geoff_S0115Tuedomagb2021$total_num

# drop
mdr1_fix$geoff_S0177Hussien2020 <- NULL

# keep specific entries
mdr1_fix$geoff_S0204Nibaptn2023 <- mdr1_fix$geoff_S0204Nibaptn2023 %>%
  filter(variant_string == "mdr1:86:Y") # these are the correct rows

mdr1_fix <- do.call(rbind, mdr1_fix) # removed 34 rows

# make a dataset with the manual fixes combined
mdr1 <- mdr1_original %>% 
  mutate(variant_string = standardise_mdr1_86(variant_string)) %>% 
  filter(variant_string != "mdr1:86:F") %>%
  dplyr::group_by(survey_ID, collection_day) %>%
  mutate(n = length(unique(total_num))) %>%
  filter(n == 1) %>% # exclude surveys with different Ns
  ungroup() %>%
  rbind(mdr1_fix) %>%
  group_by(survey_ID, collection_day, total_num, variant_string) %>%
  mutate(variant_num = sum(variant_num),
         total_num = unique(total_num)) %>%
  ungroup() %>%
  distinct(survey_ID, collection_day, total_num, variant_string,.keep_all = TRUE) %>% 
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(num_variants = n())

# split by the number of variants in the survey
mdr1_1 <- mdr1 %>% filter(num_variants == 1)
mdr1_2 <- mdr1 %>% filter(num_variants == 2)
mdr1_3 <- mdr1 %>% filter(num_variants == 3)

# start with cleaning the surveys where we only have one variant
mdr1_1 %>% group_by(variant_string) %>% summarise(n = n()) # all are mdr1:86:Y only

mdr1_1_n <- mdr1_1 %>% filter(variant_string == "mdr1:86:N")
mdr1_1_y <- mdr1_1 %>% filter(variant_string == "mdr1:86:Y")

# for each row in mdr1_1_n add an additional row that is the complement
complement <- NULL
for(i in 1:nrow(mdr1_1_n)) {
  df <- mdr1_1_n[i,]
  df$variant_num <- df$total_num - df$variant_num
  df$variant_string <- "mdr1:86:Y"
  df$prev <- df$variant_num / df$total_num
  complement <- rbind(complement, df)
}

# checked that sum(variants) == denom
mdr1_1_n <- rbind(mdr1_1_n, complement) %>%
  arrange(survey_ID, collection_day, total_num) %>%
  mutate(prev = variant_num/total_num)

# for each row in mdr1_1_y add an additional row that is the complement
complement <- NULL
for(i in 1:nrow(mdr1_1_y)) {
  df <- mdr1_1_y[i,]
  df$variant_num <- df$total_num - df$variant_num
  df$variant_string <- "mdr1:86:N"
  df$prev <- df$variant_num / df$total_num
  complement <- rbind(complement, df)
}

mdr1_1_y <- rbind(mdr1_1_y, complement) %>%
  arrange(survey_ID, collection_day, total_num) %>%
  mutate(prev = variant_num/total_num)

mdr1_1 <- rbind(mdr1_1_n, mdr1_1_y) %>%
  mutate(prev = variant_num / total_num) 

# checked that sum(variants) == denom
mdr1_1 %>%
  group_by(collection_day, survey_ID, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>%
  nrow() 
# [1] 0

# now clean those that report all three because this is another simple case
mdr1_3 <- mdr1 %>% 
  filter(num_variants == 3) %>%
  group_by(collection_day, survey_ID, total_num) %>%
  mutate(xs = sum(variant_num)) 

mdr1_3_correct <- mdr1_3 %>%
  filter(xs == total_num)

# TODO: one study that still needs fixing

mdr1_3 <- mdr1_3_correct

# now clean those with 2 mutations
# 3 possible iterations of 2 mutations
# add a column that indicates which of the three classes it falls into so that we can treat them separately
mdr1_2_groups <- mdr1_2 %>%
  group_by(nid) %>%
  mutate(
    pair = unique(variant_string) |> str_c(collapse = "+"),
    pair_simple = case_when(
      pair %in% c("mdr1:86:N+mdr1:86:Y", "mdr1:86:Y+mdr1:86:N")   ~ "N+Y",
      pair %in% c("mdr1:86:N+mdr1:86:N/Y", "mdr1:86:N/Y+mdr1:86:N") ~ "N+N/Y",
      pair %in% c("mdr1:86:Y+mdr1:86:N/Y", "mdr1:86:N/Y+mdr1:86:Y") ~ "Y+N/Y",
      TRUE                                                          ~ pair
    )
  ) %>%
  ungroup() 


mdr1_2_groups <-  mdr1_2_groups %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  mutate(correct = if_else(xs == total_num, "yes", "no"))

# these rows are already correct and not missing imputation -- leave as is 
mdr1_2_groups_correct <- mdr1_2_groups %>% filter(correct == "yes")
  
# these are the rows that need fixing -- fix based on the classification
mdr1_2_groups <- mdr1_2_groups %>% filter(correct == "no") 
  
# split the dataframe with two unique variants into what combination of variants - M ~ mixed
mdr1_2_NY <- mdr1_2_groups %>% filter(pair_simple == "N+Y")
mdr1_2_NM <- mdr1_2_groups %>% filter(pair_simple == "N+N/Y")
mdr1_2_YM <- mdr1_2_groups %>% filter(pair_simple == "Y+N/Y")

# NM = 2 rows
mdr1_2_NM[3,] <- mdr1_2_NM[2,]
mdr1_2_NM[3,]$variant_num <- mdr1_2_NM[3,]$total_num - sum(mdr1_2_NM$variant_num[1:2])
mdr1_2_NM[3,]$variant_string <- "mdr1:86:Y"
# check that sum numerator = denominator now 
mdr1_2_NM %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow()

# now looking at N & Y
mdr1_2_NY <- mdr1_2_NY %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(sid = cur_group_id()) %>%
  arrange(sid)

mdr1_2_NY_impute <- NULL
for(i in 1:length(unique(mdr1_2_NY$sid))) {
  df <- mdr1_2_NY[c((2*i)-1, 2*i),]
  df[3,] <- df[2,]
  df[3,]$variant_string <- "mdr1:86:N/Y"
  df[3,]$variant_num <- df[3,]$total_num - df[3,]$xs
  df[3,]$prev <- df[3,]$variant_num / df[3,]$total_num
  mdr1_2_NY_impute <- rbind(df, mdr1_2_NY_impute)
}
mdr1_2_NY_impute %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow() # check these are all fixed

mdr1_2_NY <- mdr1_2_NY_impute

# now clean YM
mdr1_2_YM <- mdr1_2_YM %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(sid = cur_group_id()) %>%
  arrange(sid)

mdr1_2_YM_impute <- NULL
for(i in 1:length(unique(mdr1_2_YM$sid))) {
  df <- mdr1_2_YM[c((2*i)-1, 2*i),]
  df[3,] <- df[2,]
  df[3,]$variant_string <- "mdr1:86:N/Y"
  df[3,]$variant_num <- df[3,]$total_num - df[3,]$xs
  df[3,]$prev <- df[3,]$variant_num / df[3,]$total_num
  mdr1_2_YM_impute <- rbind(df, mdr1_2_YM_impute)
}

mdr1_2_YM_impute %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow() # check these are all fixed

mdr1_2_YM <- mdr1_2_YM_impute

mdr1 <- rbind(mdr1_1,
                mdr1_3, 
                mdr1_2_NM,
                mdr1_2_NY,
                mdr1_2_YM) %>%
  select(names(master_table_simplified)) %>%
  mutate(prev = variant_num / total_num)
  
gene_keep <- rbind(gene_keep, mdr1)

# Save the cleaned and formatted data
saveRDS(gene_keep, here("analysis", "data-derived", "geoff_clean.rds"))
saveRDS(master_table_formatted, here("analysis", "data-derived", "geoff_clean_complete.rds"))
