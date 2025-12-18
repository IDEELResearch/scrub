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

# Read in shape file data
gpkg_africa_shp <- sf::st_read("analysis/data-raw/gpkg_africa_shape_file.gpkg") %>%
  rename(geometry = geom)

# Obtain Africa iso3 country names
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

# Make our full bind across
column_names <- get_column_names_for_clean()
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

# Rename lat lon to keep after sf::join()
full_bind <- full_bind %>% mutate(lon_keep = lon, lat_keep = lat)
# Convert the coordinate class of the geometry in the combined cleaned df and the shape file to the same class
full_bind_sf <- sf::st_as_sf(full_bind, coords = c("lon", "lat"), crs = sf::st_crs(gpkg_africa_shp)) %>%
  rename(lat = lat_keep, lon = lon_keep)

# Merge combined cleaned df and the shape filebased on the geometry
merged_full_bind_sf <- sf::st_join(full_bind_sf, gpkg_africa_shp, left = TRUE)

# iso3c does not match iso3 from geometry
mismatch_iso3c <- merged_full_bind_sf %>% filter(iso3c != GID_0)

# Check which survey_IDs are NA
na_admin2 <- merged_full_bind_sf %>% filter(is.na(NAME_2))

# Manually checked NA admin 2 and mismatch of iso3
full_bind_cleaned_lon_lat <- full_bind_sf %>%
  mutate(lon = as.numeric(lon), lat = as.numeric(lat)) %>%
  mutate(
    lon = case_when(
      survey_ID == "who_1171_NA_2011" ~ -17.321190,
      survey_ID == "who_1491_NA_2017" ~ -23.505476720688392,
      survey_ID == "who_182_NA_2013" ~ 45.22758346654867,
      survey_ID == "who_184_Ngazidja_2013" ~ 43.375,
      survey_ID == "who_246_NA_2013" ~ -16.58865735038155,
      survey_ID == "who_879_NA_2014" ~ 9.29690,
      survey_ID == "who_888_Benishangul_GumuzandAmhara_2014" ~ 38.7652126480875,
      survey_ID == "who_1223_ThiesRegion_2019" ~ -16.98590086453347,
      TRUE ~ lon
    ),
    lat = case_when(
      survey_ID == "who_1171_NA_2011" ~ 14.783380,
      survey_ID == "who_1491_NA_2017" ~ 14.920444770903352,
      survey_ID == "who_182_NA_2013" ~ -12.77878320419641,
      survey_ID == "who_184_Ngazidja_2013" ~ -11.5,
      survey_ID == "who_246_NA_2013" ~ 13.459448211037648,
      survey_ID == "who_879_NA_2014" ~ -1.685000,
      survey_ID == "who_888_Benishangul_GumuzandAmhara_2014" ~ 9.043439082761527,
      survey_ID == "who_1223_ThiesRegion_2019" ~ 14.418164041698063,
      TRUE ~ lat
    )
  ) %>%
  filter(survey_ID != "who_1223_ThiesRegion_2019")
 
# Update geometry in full_bind_cleaned_lon_lat
sf::st_geometry(full_bind_cleaned_lon_lat) <- sf::st_as_sf(
  full_bind_cleaned_lon_lat %>%
    sf::st_drop_geometry(), # remove current geometry
  coords = c("lon", "lat"),
  crs = sf::st_crs(full_bind_cleaned_lon_lat_sf)) %>%
  sf::st_geometry()
# Rename lat lon to keep after sf::join()
full_bind_cleaned_lon_lat <- full_bind_cleaned_lon_lat %>% mutate(lon_keep = lon, lat_keep = lat)
# Convert the coordinate class of the geometry in the combined cleaned df and the shape file to the same class
full_bind_cleaned_lon_lat_sf <- full_bind_cleaned_lon_lat %>%
  select(-lon, -lat) %>%
  rename(lat = lat_keep, lon = lon_keep) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = sf::st_crs(gpkg_africa_shp))

# Merge combined cleaned df and the shape filebased on the geometry
merged_full_bind_cleaned_lon_lat <- sf::st_join(full_bind_cleaned_lon_lat_sf, gpkg_africa_shp, left = TRUE)

# iso3c does not match iso3 from geometry
mismatch_iso3c_cleaned <- merged_full_bind_cleaned_lon_lat %>% filter(iso3c != GID_0)

# Check which survey_IDs are NA
na_admin2_cleaned <- merged_full_bind_cleaned_lon_lat %>% filter(is.na(NAME_2))

dim(mismatch_iso3c)
dim(mismatch_iso3c_cleaned)
dim(na_admin2)
dim(na_admin2_cleaned)
