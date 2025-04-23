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
  # OJ checked these and worked out their site study type
  dplyr::mutate(site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "32459360", "COMMUNITY")) %>% 
  dplyr::mutate(site_study_type = replace(site_study_type, is.na(site_study_type) & pmid == "37670357", "COMMUNITY"))

# additional cleaning to correct for malformed iso3cs
master_table_clean <- master_table_clean %>%
  mutate(country = if_else(country %in%  c("Ago", "Ben", "Cmr", "Caf", "Cog", 
                                           "Cod", "Eth", "Gha", "Stp", "Sen", 
                                           "Sdn", "Uga", "Gnq", "Gab", "Mli", 
                                           "Nga", "Sle", "Som", "Zmb", "Zaf"),
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

# 2. Step 2 - Clean mutation names and data entry issues in these # ------------

# Clean mutation names

# Expand all gene mutation ranges for reference range syntax
# TODO: These are inconcistently entered. Some use mutant_num to refer to REF ranges
# Some go the other way. Plus many data entry issues here on review. 
# Leaving for now and will pass back to data entry
indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))
master_table_clean$gene_mutation[indices_to_transform] <- sapply(
  master_table_clean$gene_mutation[indices_to_transform],
  collapse_k13_range, 
  mutation_key = mutation_key
)

# For an individual asterisk with no mutants, we simply flip to ref and set all to total
indices_to_transform <- which(grepl("^k13:[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)) &
                                master_table_clean$mutant_num == 0)
master_table_clean$gene_mutation[indices_to_transform] <- as.character(sapply(
  master_table_clean$gene_mutation[indices_to_transform],
  convert_k13_asterisk, 
  mutation_key = mutation_key
))
master_table_clean$mutant_num[indices_to_transform] <- master_table_clean$total_num[indices_to_transform]

# OJ: This last one had a STOP codon based on the paper
# https://journals.plos.org/plosone/article/figure?id=10.1371/journal.pone.0235401.t006
# No way to handle so remove this - there are other data points in this location so 
# the denominator will be correctly captured still
indices_to_transform <- which(grepl("^k13:[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))
master_table_clean <- master_table_clean[-indices_to_transform,]

# OJ: They did enter FC which is what the paper says. This is then likely an insertion rather than a SNP
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-024-04945-8/tables/3
# However, all entires had 0 observations of this, so we cannot simply remove as that removes
# the denominator for this locus. So set to wildtype and mutant_num to total_num
master_table_clean$mutant_num[master_table_clean$gene_mutation == "MDR1:1034:FC"] <- 
  master_table_clean$total_num[master_table_clean$gene_mutation == "MDR1:1034:FC"]
master_table_clean$gene_mutation[master_table_clean$gene_mutation == "MDR1:1034:FC"] <- "MDR1:1034:S"

# OJ: They entered insertions, which are not SNPs so we will remove these 
# This will remove the locus entirely but for now this is okay as it is outside propeller (may want to revisit)
# https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-472/tables/1
master_table_clean <- master_table_clean %>% filter(!(gene_mutation %in% c("K13:142:NN", "K13:142:NNN")))

# Correct gene mutation format
convmutation <- format_variants_for_stave(unique(master_table_clean$gene_mutation))
master_table_clean$gene_mutation <- convmutation[match(master_table_clean$gene_mutation, unique(master_table_clean$gene_mutation))]

# 3. Step 3 - Clean dates, lat, lon # ------------

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
# This has now been manually corrected higehr up but leave for now just in case reappears
master_table_clean$lat_n <- clean_excel_formulas(master_table_clean$lat_n)
master_table_clean$lon_e <- clean_excel_formulas(master_table_clean$lon_e)

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

# 4. Step 4 - Pre formatting other data entry issues # ------------

# One manual fix to data entry issue in 
# s0128_voumbo-matoumona_2018_v01_prevalence_data_LONG_validated.tsv.gz
master_table_clean <- master_table_clean %>% 
  mutate(mutant_num = replace(mutant_num, mutant_num == 31 & total_num == 24 & site_uid == "koulamoutou", 21))

# remove site="none" records from study_ID==S0006WrairKenreadKen for which lat and lon was not available
# the meta file for this has a none location so it passed through
# automatically from MIP pipeline reading metadata for samples without site info which have no lat/lon. 
# To be fixed manually in entry eventually.
master_table_clean <- master_table_clean[!(master_table_clean$study_uid == "s0006_wrair_kenread_ken" & 
                                                     is.na(master_table_clean$lat_n) & 
                                                     is.na(master_table_clean$lon_e)), ]

# filter out any entries that are day0 calculated - these reflect extra data entries for total genotype counts
# which are already captured within haplotype counts
master_table_clean <- master_table_clean %>% 
  dplyr::filter(grepl("day0calculated", substudy) == 0)


# 5. Step 5 - Formatting # ------------

if(sum(master_table_clean$total_num < master_table_clean$mutant_num) > 0) {
  errorCondition("Entries with more mutant samples than total samples")
}
if(sum(master_table_clean$total_num == 0) > 0) {
  errorCondition("Entries with zero total samples")
}

# rename the columns so that they make sense
column_names <- get_column_names_for_clean()

# first add in extra information needed within column_names
master_table_formatted <- master_table_clean %>% 
  group_by(study_uid) %>% 
  dplyr::mutate(study_uid = paste0("geoff_", janitor::make_clean_names(study_uid[1], "big_camel"))) %>% 
  ungroup() %>% 
  dplyr::group_by(study_uid, site_name, collection_start) %>% 
  dplyr::mutate(survey_ID = paste0(study_uid[1], "_", 
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
  )


# 6. Step 6 - Impute k13 reference survey records not explicitly entered in data entry # ------------

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

# Check to confirm all k13_alleles (24 validated and candidate) are present for each survey (year/location) at least (should return 0 rows)
master_table_formatted %>% filter(gene == "k13") %>% group_by(survey_ID) %>% 
  summarise(n = str_extract_all(variant_string, "(?<=k13:)([\\d_]+)") %>%
                unlist() %>%
                str_split("_") %>%
                unlist() %>%
                as.numeric() %>%
                unique() %>% length,
            min = as.integer(unique(k13_min)), 
            max = as.integer(unique(k13_max)), 
            pos = sum(as.integer(unique(substr(k13_mutations$mut,1,3))) < max & 
                        as.integer(unique(substr(k13_mutations$mut,1,3))) > min)) %>% 
  filter(n < pos)

# 7. Step 7 - Last checks # ------------

# flag any study which has NA for lat or long after format conversion
lat_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lat)])
lon_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$lon)])

# flag any study which has NA for substudy which is selectively used for inclusion of data in next scripts (ex we must use this field to exclude categories like tesday0treatmentfailure or treatedextracted for population prevalence calcs)
substudy_missing_or_improper <- unique(master_table_formatted$study_ID[is.na(master_table_formatted$substudy)])
cat("Studies with missing or improper latitudes:\n", paste(lat_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper longitudes:\n", paste(lon_missing_or_improper, collapse = ", "), "\n\n")
cat("Studies with missing or improper substudy values (critical for inclusion/exclusion in downstream analysis):\n", 
    paste(substudy_missing_or_improper, collapse = ", "), "\n\n")

# flag all study uids which report less than 40% prevalence of validated k13 reference allele indicating entry error (using 469C arbitratily to detect erroneous user behavior) 
low_prev_studies <- master_table_formatted %>%
  filter(variant_string == "k13:469:C", prev < 0.40) %>%
  distinct(study_ID, data_entry_author)  # Get unique pairs
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

# 9. Step 9 - Save Formatted Data ---------

# Grab just the columns we need for pairing with WWARN etc
master_table_simplified <- master_table_formatted %>% 
  select(all_of(column_names))

# Save the final merged_df as an RDS file
saveRDS(master_table_simplified, here("analysis", "data-derived", "geoff_clean.rds"))
saveRDS(master_table_formatted, here("analysis", "data-derived", "geoff_clean_complete.rds"))
