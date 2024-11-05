
################################################################################
#
# Script: 02_clean_geoffs.R
# Purpose: Perform various data cleaning to allow deduplication and combining
# with other data sources within {stave} object
# 
################################################################################

# read in data
# Load the combined geoff data table created in the first script
master_table <- readRDS(here("analysis", "data-derived", 
                             "01_read_geoffs_output_table.rds"))

# Load mutation key for k13 reference ranges
mutation_key <- read.csv(here("analysis", "data-raw", 
                              "k13_ref_protein_codon_dictionary.csv"))

## start by looking at the data - explore this in the console
master_unique <- lapply(master_table, unique)

# TODO: fix collection location entry: this is to do with how it's read in 01 
# "MBanza Congo Municipal Hospital\tAngola\tAGO\t-6.265461256\t14.2530007\tMBanza Congo Municipal Hospital"
# refers to these study keys:"s0050_davlantes_v01" "s0052_ljolje_v01"  # have asked team to investigate
# TODO: fix all "na" to be NA

master_table_clean <- master_table |>
  dplyr::mutate(substudy = gsub(" ","", substudy)) |>
  dplyr::mutate(substudy = gsub("_","", substudy)) |>
  dplyr::mutate(substudy = if_else(substudy == "untreadextracted",
                                   "untreatedextracted", substudy)) |>
  dplyr::mutate(substudy = if_else((substudy == "day0" & study_uid == "s0029_warsame_2019"),
                                   "day0extracted", substudy)) |>
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
  dplyr::mutate(publication_year = if_else(publication_status == "unpublished", NA, publication_year)) |>
  dplyr::mutate(countries_covered = gsub(",", ", ", countries_covered)) |>
  dplyr::mutate(countries_covered = stringr::str_to_title(countries_covered)) |>
  dplyr::mutate(countries_covered = if_else(countries_covered %in% c("Democratic Republic Of The Congo", "Drc"),
                                            "Democratic Republic of the Congo",
                                            country)) |>
  dplyr::mutate(data_processing_pipeline = toupper(data_processing_pipeline)) |>
  dplyr::mutate(first_author_surname = stringr::str_to_title(first_author_surname))

# perform tests
african_countries <- data.frame(country = countrycode::codelist$country.name.en,
                                continent = countrycode::codelist$continent,
                                iso3c = countrycode::codelist$iso3c) |>
  dplyr::filter(continent == "Africa") |> dplyr::distinct()
allowed_countries <- c(african_countries$country, "Democratic Republic of the Congo")
allowed_iso3c <- c(african_countries$iso3c, "DRC", "ERT")
allowed_substudy <- c("untreatedextracted", "untreatedcalculated",
                      "treatedextracted", "treatedcalculated",
                      "day0extracted", "day0calculated")
allowed_pretreatment <- c("no", "yes", NA)
allowed_publication <- c("peer_reviewed", "preprint", "unpublished")
allowed_pub_year <- c(as.character(seq(from = 2000, to = (lubridate::year(Sys.Date())))), NA) # automatically update as years change
allowed_site_types <- c("HEALTHFACILITY", "COMMUNITY", "TES", "CCS", "DHS" )

# run tests
check_values_in_column(master_table_clean, "country", allowed_countries)
check_values_in_column(master_table_clean, "iso3c", allowed_iso3c)
check_values_in_column(master_table_clean, "substudy", allowed_substudy)
check_values_in_column(master_table_clean, "pretreatment_samples", allowed_pretreatment)
check_values_in_column(master_table_clean, "publication_status", allowed_publication)
check_values_in_column(master_table_clean, "publication_year", allowed_pub_year)
check_values_in_column(master_table_clean, "site_study_type", allowed_site_types)

#TODO: check the extracted vs calculated and filter accordingly



# Clean mutation names
# Expand gene mutation ranges for reference range syntax
indices_to_transform <- which(grepl("^k13:[0-9]+-[0-9]+:\\*$", tolower(master_table_clean$gene_mutation)))
master_table_clean$gene_mutation[indices_to_transform] <- sapply(
  master_table_clean$gene_mutation[indices_to_transform],
  collapse_k13_range
)

# Adjust and add dates and survey IDs
master_table_clean$collection_start <- sapply(master_table_clean$date_start, adjust_invalid_date, is_start = TRUE)
master_table_clean$collection_end <- sapply(master_table_clean$date_end, adjust_invalid_date, is_start = FALSE)

# Compute collection day
master_table_clean$collection_day <- sapply(1:nrow(master_table_clean), function(i) {
  if (!is.na(master_table_clean$collection_start[i]) & !is.na(master_table_clean$collection_end[i])) {
    return(as.Date(median(c(as.numeric(master_table_clean$collection_start[i]), as.numeric(master_table_clean$collection_end[i]))), origin = "1970-01-01"))
  } else if (!is.na(master_table_clean$collection_start[i])) {
    return(as.Date(master_table_clean$collection_start[i]))
  } else if (!is.na(master_table_clean$collection_end[i])) {
    return(as.Date(master_table_clean$collection_end[i]))
  } else {
    return(as.Date(NA))
  }
})

# TODO: fix long and lat which got read in as formula vs values
fix_lat <- which(is.na(as.numeric(master_table_clean$lat_n)))
fix_lon <- which(is.na(as.numeric(master_table_clean$lon_e)))

master_table_clean$lat_n[fix_lat] <- mean(as.numeric(str_extract_all(master_table_clean$lat_n[fix_lat], 
                                                                     "\\d+\\.\\d+")[[1]]))
master_table_clean$lon_e[fix_lat] <- mean(as.numeric(str_extract_all(master_table_clean$lon_e[fix_lat], 
                                                                     "\\d+\\.\\d+")[[1]]))
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
wwarn_names <- c("iso3c", "admin_0", "admin_1", "site", "lat", "long", "year", 
                 "study_start_year", "study_end_year", "x", "n", "prev", 
                 "gene", "mut", "gene_mut", "annotation", "database", 
                 "pmid", "url", "source")

master_table_clean <- master_table_clean |>
  dplyr::mutate(database = "GEOFF") |>
  dplyr::rename(admin_0 = "country",
                # admin_1 = NA,
                site = "site_name",
                lat = "lat_n",
                long = "lon_e",
                year = "collection_day",
                study_start_year = "collection_start", 
                study_end_year = "collection_end", 
                x = "mutant_num",
                n = "total_num",
                # prev = NA,
                gene = "gene",
                mut = "mut",
                gene_mut = "gene_mutation",
                # annotation = NA,
                database = "database",
                url = "study_url",
                source = "publication_status",
                iso3c = "iso3c",
                pmid = "pmid") |>
  dplyr::mutate(admin_1 = NA,
                prev = NA, 
                annotation = NA) |> #reorder to have the same column order as WWARN
  dplyr::relocate(all_of(wwarn_names))


# TODO: figure out which additional columns would be helpful for inclusion here to enable better deduplication
master_table_simplified <- master_table_clean[, c(wwarn_names, "study_uid", "data_entry_author", "data_processing_pipeline")]

# Save the final merged_df as an RDS file
saveRDS(master_table_simplified, here("analysis", "data-derived", "02_clean_geoffs_simple.RDS"))
saveRDS(master_table_clean, here("analysis", "data-derived", "02_clean_geoffs_complete.RDS"))
print("Data saved. Ready for further analysis.")
