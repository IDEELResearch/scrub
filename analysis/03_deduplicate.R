library(tidyverse)
library(here)

# source functions
source("R/adjust_invalid_date.R")
source("R/deduplicate.R")

# Read each file if it exists
safe_read <- function(path) {
  if (file.exists(path)) {
    clean <- readRDS(path) %>% 
      as.data.frame() 
    if ("study_type" %in% names(clean)) {
      clean <- clean %>% 
        dplyr::filter(study_type != "unpublished") %>% 
        dplyr::filter(study_type != "data_incomplete")
        
    }
  } else {
    clean <- data.frame()
  }
  return(clean)
}

# Function to add tags based on the criteria
add_tags <- function(df) {
  df %>%
    mutate(keep_delete = case_when(database == "Pf7k" & any(database == "WHO") ~ "Keep",
                                   database == "GEOFF" & (any(database == "WHO") | any(database == "Pfk7")) ~ "Keep",
                                   TRUE ~ "Delete"))}

dedpuplicate_study <- function(df){
  ### **Scenario 1: Identify Different `study_ID`s Reporting the Same Data**
  same_data_diff_study <- df %>%
    group_by(name_2, variant_string, collection_start, collection_end) %>%
    filter(n_distinct(study_ID) > 1) %>%  # Ensure multiple studies report the same data
    ungroup()
  
  # Split the dataframe into a list by the grouping columns
  duplicate_diff_study_list <- same_data_diff_study %>%
    group_by(name_2, variant_string, collection_start, collection_end) %>%
    filter(n() > 1) %>%  # Double check that groups have >1 entry
    group_split()
  
  # Apply the function to each dataframe in the list
  tagged_duplicate_diff_study_list <- lapply(duplicate_diff_study_list, add_tags)
  
  # Example of how to view the first dataframe in the list
  head(tagged_duplicate_diff_study_list[[1]])
  
  ### TO-DO CECILE: SINCE WWARN TAKE MID POINT ACROSS MULTIPLE SITES IT MIGHT NOT BE AS ACCURATE WITH ADMIN2
  ### NEED TO TALK TO OJ IF THIS METHOD SHOULD SURFFICE OR NOT
  
  ### **Scenario 3: Identify Duplicate Data Entries in the Same Study**
  same_data_same_study <- df %>%
    group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
    filter(n() > 1) %>%  # Ensure duplicate rows exist
    mutate(keep_row = dplyr::row_number() == 1) %>%
    ungroup()
  
  # Split the dataframe into a list by the grouping columns
  duplicate_same_study_list <- same_data_same_study %>%
    group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
    filter(n() > 1) %>%  # Double check that groups have >1 entry
    group_split()
  
  # TO-DO CECILE: REMOVE NON WHO STUDIES AND KEEP ONE ISTANCE OF WHO STUDY
 
  # ### **Scenario 4: Identify Duplicate Data Entries in the Same Study and unique site names (more stringent)**
  # same_data_same_study_site_name <- df %>%
  #   group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
  #   filter(n() > 1) %>%  # Ensure duplicate rows exist
  #   mutate(keep_row = dplyr::row_number() == 1) %>%
  #   ungroup()
  # 
  # # Split the dataframe into a list by the grouping columns
  # duplicate_same_study_site_name_list <- same_data_same_study_site_name %>%
  #   group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
  #   filter(n() > 1) %>%  # Double check that groups have >1 entry
  #   group_split()
  
  return(list(diff_study = duplicate_diff_study_list, 
              same_study = duplicate_same_study_list))
}

# Read in each cleaned file
clean_geoff <- safe_read(here("analysis", "data-derived", "geoff_clean.rds"))
clean_wwarn <- safe_read(here("analysis", "data-derived", "wwarn_clean.rds"))
clean_pf7k <- safe_read(here("analysis", "data-derived", "pf7k_clean.rds"))
clean_who <- safe_read(here("analysis", "data-derived", "who_clean.rds"))

africa_iso3 <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD",
  "COM", "COD", "COG", "DJI", "EGY", "GNQ", "ERI", "SWZ", "ETH", "GAB",
  "GMB", "GHA", "GIN", "GNB", "CIV", "KEN", "LSO", "LBR", "LBY", "MDG",
  "MWI", "MLI", "MRT", "MUS", "MAR", "MOZ", "NAM", "NER", "NGA", "RWA",
  "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TGO", "TUN",
  "UGA", "TZA", "ZMB", "ZWE", "ESH"
)

######## CECILE TO-DO: DELETE IF clean_pf7k, clean_who will only be Africa countries in the future ########
clean_pf7k <- clean_pf7k %>% filter(iso3c %in% africa_iso3)
clean_who <- clean_who %>% filter(iso3c %in% africa_iso3)

######## CECILE TO-DO: CHECK WITH OJ HOW WE SHOULD HANDLE LOADING OF ADMIN 1/ADMIN2 data ########
# Use `malariaAtlas` R package to get shape files for Admin 1 and Admin 2
# sf_adm1_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin1", version = "202206")
# saveRDS(sf_adm1_africa, "analysis/data-derived/sf_admin1_africa.rds")
# sf_adm2_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin2", version = "202206")
# saveRDS(sf_adm2_africa, "analysis/data-derived/sf_admin2_africa.rds")
# sf_adm0_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin0")
# saveRDS(sf_adm0_africa, "analysis/data-derived/sf_admin0_africa.rds")

# Read in malariaAtlas shape files for admin 1 and admin 2
sf_adm0_africa <- readRDS("analysis/data-derived/sf_admin0_africa.rds")
sf_adm1_africa <- readRDS("analysis/data-derived/sf_admin1_africa.rds")
sf_adm2_africa <- readRDS("analysis/data-derived/sf_admin2_africa.rds")

# make our full bind across
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

df <- full_bind
admin2_df <- sf_adm2_africa

df_sf <- df %>%
  filter(!is.na(lon) & !is.na(lat)) %>% # Delete rows where lat/long is NA - NOTE: study S0147Jeang2024 is missing lat
  mutate(lon_keep = lon, lat_keep = lat) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # Ensure lat lon coordinates use CRS 4326
  sf::st_transform(df_sf, crs = sf::st_crs(admin2_df)) # Ensure both dataframes use the same CRS (WGS 84)

# Perform spatial join: Match points (df_sf) to Admin 2 polygons (admin2_df)
sf::sf_use_s2(FALSE)
df_sf_joined <- df_sf %>%
  sf::st_join(admin2_df %>% 
                rename(admin2_iso3c = iso) %>% 
                select(name_2, admin2_iso3c, geometry), 
              left = TRUE)

#------------------ ISO3 AND ADMIN2 CHECKS ------------------
######## CECILE TO-DO: REMOVE AT SOME POINT ########
# These studies have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
mismatch_rows <- df_sf_joined %>%
  filter(iso3c != admin2_iso3c) %>%  # Compare entry_iso3c with admin0_iso3c
  select(lon_keep, lat_keep, site_name, site_name, name_2, iso3c, admin2_iso3c, study_ID, survey_ID)

######## CECILE TO-DO: REMOVE AT SOME POINT ########
# These studies have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
na_admin <- df_sf_joined %>%
  filter(is.na(admin2_iso3c)) %>%
  select(lon_keep, lat_keep, site_name, site_name, name_2, iso3c, admin2_iso3c, study_ID, survey_ID)

# Remove studies that need to be fixed and have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
df_sf_joined <- df_sf_joined %>%
  filter(!(survey_ID %in% na_admin$survey_ID | survey_ID %in% mismatch_rows$survey_ID))

#------------------ DEDUPLICATION ------------------
# Convert back to a regular dataframe
df_final <- df_sf_joined %>% 
  sf::st_drop_geometry() %>% 
  mutate(lon = lon_keep, lat = lat_keep) %>%  # Restore lat/lon
  select(-lon_keep, -lat_keep)  # Remove temporary columns

# 1. approach could be to add extra columns to help with deduplication
# e.g. if collection start and end we think is too specific then perhaps year
df_final$collection_year_start <- lubridate::year(df_final$collection_start)
df_final$collection_year_end <- lubridate::year(df_final$collection_end)

# identify studies that may be duplicates
dedup_output = dedpuplicate_study(df_final)
duplicate_diff_study_list = dedup_output$diff_study
duplicate_same_study_list = dedup_output$same_study

length(duplicate_diff_study_list)
length(duplicate_same_study_list)

#------------------ CHECK DIFFERENCES BETWEEN HAVING THE site_name BE THE SAME ------------------
# Check studies that are included in duplicate_diff_study_list but not in duplicate_diff_study_site_name_list
# Extract unique grouping keys for both lists
get_group_keys <- function(dup_list) {
  map_dfr(dup_list, ~ select(.x, name_2, variant_string, collection_start, collection_end) %>% distinct())
}

diff_study_keys <- get_group_keys(duplicate_diff_study_list)
diff_study_site_name_keys <- get_group_keys(duplicate_diff_study_site_name_list)

# Find the groups that are in diff_study_keys but not in diff_study_site_name_keys
extra_duplicate_keys <- anti_join(diff_study_keys, diff_study_site_name_keys, 
                                  by = c("name_2", "variant_string", "collection_start", "collection_end"))

# Filter `duplicate_diff_study_list` to keep only the extra duplicates
extra_duplicate_list <- keep(duplicate_diff_study_list, function(df) {
  any(inner_join(df, extra_duplicate_keys, by = c("name_2", "variant_string", "collection_start", "collection_end")) %>% nrow() > 0)
})

# Print the extra duplicates
length(extra_duplicate_list)













# ------

# 2. there is a little helper function to show what gets removed by a specific set of columns
rm_ex <- rows_removed_by_distinct(full_bind, url, country_name, lat, lon, gene, 
                                  variant_string, total_num, collection_year_start, 
                                  collection_year_end)

duplicated_sum <- generate_duplication_summary(rm_ex)

# this could then help possibly to help identify data entry issues e.g.
# (which perhaps we should just be grouping by and summing)
as.data.frame(rm_ex[[1]])

# or it could identify that perhaps filtering based on year is too stringent 
# e.g. the below data entry is probably correct from within the same study 
as.data.frame(rm_ex[[2]])

# and in fact all of these are only 2 rows long so likely these are fine
lengths(rm_ex)

# in truth Cecile - this deduplication will likely only really come in when WWARN
# is in which does have large crossover i suspect with WHO, so perhaps just using
# the approach above to help identify possible errors of the first kind that can be 
# grouped and summed would be very helpful

# 3. once we have all the columns we want then the function that is called is
# below so add into this function defined in scrub for this
# deduplication
full_bind <- deduplicate(full_bind)

# save ready to go to stave
saveRDS(full_bind, here("analysis/data-derived/final_data.rds"))

