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

### TO-DO CECILE: DELETE IF clean_pf7k, clean_who will only be Africa countries in the future
clean_pf7k <- clean_pf7k %>% filter(iso3c %in% africa_iso3)
clean_who <- clean_who %>% filter(iso3c %in% africa_iso3)

### TO-DO CECILE: CHECK WITH OJ HOW WE SHOULD HANDLE LOADING OF ADMIN 1/ADMIN2 data ########
# Use `malariaAtlas` R package to get shape files for Admin 1 and Admin 2
# sf_adm1_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin1", version = "202206")
# saveRDS(sf_adm1_africa, "analysis/data-derived/sf_admin1_africa.rds")
# sf_adm2_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin2", version = "202206")
# saveRDS(sf_adm2_africa, "analysis/data-derived/sf_admin2_africa.rds")
# sf_adm0_africa <- malariaAtlas::getShp(ISO = africa_iso3, admin_level = "admin0")
# saveRDS(sf_adm0_africa, "analysis/data-derived/sf_admin0_africa.rds")

# Read in malariaAtlas shape files for admin 1 and admin 2
# sf_adm0_africa <- readRDS("analysis/data-derived/sf_admin0_africa.rds")
# sf_adm1_africa <- readRDS("analysis/data-derived/sf_admin1_africa.rds")
# sf_adm2_africa <- readRDS("analysis/data-derived/sf_admin2_africa.rds")

# Read in admin2 data loaded from geodata
geodata_adm2_africa <- readRDS("analysis/data-derived/geodata_admin2_africa.rds")
geodata_admin2_sf <- geodata_adm2_africa %>%
  rename(admin2_iso3c = GID_0,
         name_2 = NAME_2) %>%
  select(name_2, admin2_iso3c, geometry)

# make our full bind across
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

#-------------------- Deconstruct spatial join
df <- full_bind
admin2_sf <- geodata_adm2_africa
df_sf <- df %>%
  filter(!is.na(lon) & !is.na(lat)) %>% # Delete rows where lat/long is NA - NOTE: study S0147Jeang2024 is missing lat
  mutate(lon_keep = lon, lat_keep = lat) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # Ensure lat lon coordinates use CRS 4326
  sf::st_transform(crs = sf::st_crs(admin2_sf)) # Ensure both dataframes use the same CRS (WGS 84)

# Perform spatial join: Match points (df_sf) to Admin 2 polygons (admin2_sf)
sf::sf_use_s2(FALSE)
# 1. Perform initial spatial join
df_sf_joined <- sf::st_join(df_sf, admin2_sf, left = TRUE)

# 2. Identify unmatched points
unmatched <- df_sf_joined %>% filter(is.na(name_2))

# 3. Get nearest polygon for unmatched points
if (nrow(unmatched) > 0) {
  nearest_idx <- sf::st_nearest_feature(unmatched, admin2_sf)
  nearest_matches <- admin2_sf[nearest_idx, ]
  
  # Compare iso codes and keep name_2 only if they match
  name_2_corrected <- ifelse(
    unmatched$iso == nearest_matches$admin2_iso3c,
    nearest_matches$name_2,
    NA_character_
  )
  
  # Update admin2_iso3c accordingly (optional, NA if iso doesn't match)
  admin2_iso3c_corrected <- ifelse(
    unmatched$iso == nearest_matches$admin2_iso3c,
    nearest_matches$admin2_iso3c,
    NA_character_
  )
  
  # Bind corrected columns and retain original geometry
  unmatched_filled <- unmatched %>%
    select(-name_2, -admin2_iso3c) %>%
    mutate(
      name_2 = name_2_corrected,
      admin2_iso3c = admin2_iso3c_corrected
    )
  
  # 4. Combine matched and filled unmatched results
  matched <- df_sf_joined %>% filter(!is.na(name_2))
  df_sf_joined_final <- bind_rows(matched, unmatched_filled)
  
} else {
  df_sf_joined_final <- df_sf_joined
}

# Convert back to a regular dataframe
df_final <- df_sf_joined %>% 
  sf::st_drop_geometry() %>% 
  mutate(lon = lon_keep, lat = lat_keep) %>%  # Restore lat/lon
  select(-lon_keep, -lat_keep)  # Remove temporary columns
#-------------------- END Deconstruct spatial join

# Spatially join the dataframe of all databases (GEOFF, WHO, Pf7k, WWARN) with the admin 2 region
full_bind_sf <- sf_join_admn2(full_bind, geodata_admin2_sf)

# TO-DO CECILE: improve this function
# Obtain collection year for deduplication
full_bind_sf$collection_year <- NA
for (i in 1:dim(full_bind_sf)[1]){
  y <- full_bind_sf$collection_start[i]
  full_bind_sf$collection_year[i] <- str_split(y, "-")[[1]][1]
}

#------------------------ WORKING ON DEDUPLICATE FUNCTION
df <- full_bind_sf
### Scenario 1: Identify potential duplicates reported across different study_IDs
# Group by administrative region, mutation, and collection timeframe
# Keep groups where more than one unique study_ID reports the same data
same_data_diff_study <- df %>%
  dplyr::group_by(name_2, variant_string, prev, total_num, collection_year) %>%
  dplyr::filter(n_distinct(study_ID) > 1) %>%
  dplyr::ungroup()

# Split the dataframe into a list where each list element represents a potential duplicate group
duplicate_diff_study_list <- same_data_diff_study %>%
  dplyr::group_by(name_2, variant_string, prev,  total_num, collection_year) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::group_split()

# Apply tagging function to flag which studies to keep or remove
tagged_duplicate_diff_study_list <- purrr::map(duplicate_diff_study_list, add_tags_diff_studyID)

### Scenario 2: Identify potential duplicate entries within the same study_ID
# Group by administrative region, site, mutation, collection timeframe, and study_ID
# Keep groups with more than one entry (i.e., internal duplication within the study)
same_data_same_study_prev <- df %>%
  dplyr::group_by(name_2, site_name, variant_string, prev, total_num, collection_year, study_ID) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

same_data_same_study_date <- df %>%
  dplyr::group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Split into a list where each list element is a group of internal duplicates
duplicate_same_study_list <- same_data_same_study_prev %>%
  dplyr::group_by(name_2, site_name, variant_string, prev, total_num, collection_year, study_ID) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::group_split()

# Apply logic to determine which row to keep for internal duplicates
tagged_duplicate_same_study_list <- lapply(duplicate_same_study_list, handle_same_studyID_duplicates)

### TO-DO CECILE: REMOVE EMPTY LIST ITEMS ONCE GEOFF IS FIXED
# Remove null elements from the list (these represent cases where no valid tag was applied)
tagged_duplicate_same_study_list <- tagged_duplicate_same_study_list[!sapply(tagged_duplicate_same_study_list, is.null)]

# Combine duplicates identified across studies and within studies
duplicate_list <- c(tagged_duplicate_diff_study_list, tagged_duplicate_same_study_list)

# Convert list of duplicates into one dataframe
duplicates_df <- dplyr::bind_rows(duplicate_list)

# Get all rows from df that are not in the duplicates_df
non_duplicate_rows <- df %>%
  dplyr::anti_join(duplicates_df, by = colnames(df)) %>%  # assumes all columns in `df` are relevant for identifying uniqueness
  dplyr::mutate(keep_row = "keep")

# Combine duplicates with non-duplicates into a single dataframe
deduplicate_df_with_tags <- bind_rows(duplicates_df, non_duplicate_rows)

#------------------------ END: WORKING ON DEDUPLICATE FUNCTION

# identify studies that may be duplicates
dedup_output = deduplicate(df_final)

# save ready to go to stave
saveRDS(dedup_output, here("analysis/data-derived/final_data.rds"))

#------------------ CHECKS ------------------
# Below are some checks I do to the input of the deduplication to identify any errors in data entry

#------------------ CHECK ISO3 AND ADMIN2 ------------------
### TO-DO CECILE: REMOVE AT SOME POINT ONCE WHO AND GOEFF ARE UPDATED
# These studies have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
df_sf_joined <- full_bind_sf
mismatch_rows <- df_sf_joined %>%
  filter(iso3c != admin2_iso3c) %>%  # Compare entry_iso3c with admin0_iso3c
  select(lon, lat, site_name, site_name, name_2, iso3c, admin2_iso3c, study_ID, survey_ID)

### TO-DO CECILE: REMOVE AT SOME POINT ONCE WHO AND GOEFF ARE UPDATED
### TO-DO CECILE: CHECK THAT THESE studies have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
na_admin <- df_sf_joined %>%
  filter(is.na(admin2_iso3c)) %>%
  rename(name2_geodata = name_2) %>%
  select(lon, lat, site_name, name2_geodata, iso3c, admin2_iso3c, study_ID, survey_ID)

# Convert na_admin into an sf object for spatial operations
na_admin_sf <- na_admin %>%
  filter(!is.na(lon) & !is.na(lat)) %>%  # Ensure valid coordinates
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)  # WGS 84 CRS

# Perform spatial join to find the closest available admin2 region
matched_admin <- sf::st_join(na_admin_sf, geodata_admin2_sf %>% 
                               rename(name2_malariaAtlas = name_2) %>%
                               select(name2_malariaAtlas), 
                             left = FALSE)

# Drop geometry and add matched values back into na_admin
na_admin_matched <- na_admin %>%
  left_join(matched_admin %>% 
              sf::st_drop_geometry() %>% 
              select(study_ID, survey_ID, name2_malariaAtlas), 
            by = c("study_ID", "survey_ID")) %>%
  mutate(name_2 = ifelse(!is.na(name2_malariaAtlas), name2_malariaAtlas, name2_geodata)) %>%
  select("study_ID", "survey_ID", "name2_geodata", "name2_malariaAtlas", "name_2", "site_name", "iso3c", "admin2_iso3c", "lon", "lat")

unique_iso3c <- unique(na_admin_matched$iso3c[!is.na(na_admin_matched$name_2)])
print(unique_iso3c)

# Remove studies that need to be fixed and have been added to https://docs.google.com/spreadsheets/d/1YkCO_tV6XtCUHHvTB0snqomhg0FlBUkknqn_nXhoxP8/edit?usp=sharing (2025/03/18)
df_sf_joined <- df_sf_joined %>%
  filter(!(survey_ID %in% na_admin$survey_ID | survey_ID %in% mismatch_rows$survey_ID))

### TO-DO CECILE: CHECK DIFFERENCE BETWEEN GEORGE AND MALARIAGEN ADMIN2 --> did this in analysis/create_admin2_data.R

