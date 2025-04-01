library(tidyverse)
library(here)

# source functions
devtools::load_all()

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

# make our full bind across
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

# TO-DO CECILE: improve this function
# Obtain collection year for deduplication
full_bind$collection_year_start <- lubridate::year(full_bind$collection_start)
full_bind$collection_year_end <- lubridate::year(full_bind$collection_end)

# identify studies that may be duplicates
dedup_output = deduplicate(full_bind)
dedup_df = dedup_output$df
summary_same_study = dedup_output$summary_same
summary_diff_study = dedup_output$summary_diff

dim(full_bind)
dim(dedup_df)

# save ready to go to stave
saveRDS(dedup_df, here("analysis/data-derived/final_data.rds"))

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
