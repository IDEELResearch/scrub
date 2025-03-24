library(tidyverse)
library(here)

# Define function to download and validate integrity of boundary geometries
download_admin_boundaries <- function(level, african_countries) {
  if (!level %in% c(0, 1, 2)) {
    stop("Invalid level: Use 0 for admin-0, 1 for admin-1, or 2 for admin-2.")
  }
  
  safe_gadm <- safely(function(country_code) {
    message("Downloading admin-", level, " data for: ", country_code)
    admin <- geodata::gadm(country = country_code, level = level, path = tempdir())
    admin_sf <- sf::st_as_sf(admin)  # Convert to sf object
    
    # Ensure valid geometries immediately after download
    admin_sf <- st_make_valid(admin_sf)  # Repair invalid geometries
    admin_sf <- st_transform(admin_sf, crs = 4326)  # Transform to EPSG:4326
    admin_sf
  })
  
  # Parallelized download with error handling
  future::plan(multisession)
  admin_safe <- future_map(african_countries, safe_gadm)
  future::plan(sequential)
  
  # Extract successful downloads
  admin_success <- purrr::map(admin_safe, "result")
  
  # Combine and validate geometries
  admin_combined <- do.call(rbind, admin_success)
  admin_combined <- st_make_valid(admin_combined)  # Final check and fix for all geometries
  return(admin_combined)
}

# Define African iso3 countries
africa_iso3 <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD",
  "COM", "COD", "COG", "DJI", "EGY", "GNQ", "ERI", "SWZ", "ETH", "GAB",
  "GMB", "GHA", "GIN", "GNB", "CIV", "KEN", "LSO", "LBR", "LBY", "MDG",
  "MWI", "MLI", "MRT", "MUS", "MAR", "MOZ", "NAM", "NER", "NGA", "RWA",
  "STP", "SEN", "SYC", "SLE", "SOM", "ZAF", "SSD", "SDN", "TGO", "TUN",
  "UGA", "TZA", "ZMB", "ZWE", "ESH"
)

# Read in malariaAtlas shape files for admin 1 and admin 2
sf_adm0_africa <- readRDS("analysis/data-derived/sf_admin0_africa.rds")
sf_adm1_africa <- readRDS("analysis/data-derived/sf_admin1_africa.rds")
sf_adm2_africa <- readRDS("analysis/data-derived/sf_admin2_africa.rds")

### TO-DO: ask OJ about this overlap and how to handle
# admin2_data <- download_admin_boundaries(level = 2, africa_iso3)
# saveRDS(admin2_data, "analysis/data-derived/geodata_admin2_africa.rds")
geodata_adm2_africa <- readRDS("analysis/data-derived/geodata_admin2_africa.rds")

#------------- CHECK ADMIN2 OVERLAP BETWEEN MALARIAATLAS AND GEODATA
# Get unique Admin2 regions from MalariaAtlas and Geodata
admin2_malariaAtlas <- sf_adm2_africa %>% distinct(name_2) %>% pull(name_2)
admin2_geodata <- geodata_adm2_africa %>% rename("name_2" = "NAME_2") %>% distinct(name_2) %>% pull(name_2)

# Find admin2 regions that exist in both sources (intersection)
admin2_overlap <- intersect(admin2_malariaAtlas, admin2_geodata)

# Find admin2 regions present in MalariaAtlas but missing in Geodata
admin2_only_malariaAtlas <- setdiff(admin2_malariaAtlas, admin2_geodata)

# Find admin2 regions present in Geodata but missing in MalariaAtlas
admin2_only_geodata <- setdiff(admin2_geodata, admin2_malariaAtlas)

# Print results
cat("Admin2 Regions in BOTH:\n", length(admin2_overlap), "\n")
cat("Admin2 ONLY in MalariaAtlas:\n", length(admin2_only_malariaAtlas), "\n")
cat("Admin2 ONLY in Geodata:\n", length(admin2_only_geodata), "\n")

#------------- MERGE MALARIAATLAS AND GEODATA BASED ON GEOMETRY
# Ensure both have the same CRS (WGS 84)
sf_adm2_africa <- st_transform(sf_adm2_africa, crs = st_crs(geodata_adm2_africa))

# Rename columns to avoid conflicts before merging
sf_adm2_africa_select <- sf_adm2_africa %>%
  select(iso, name_0, name_1, name_2, geometry) %>%
  rename(name0_malariaAtlas = name_0, name1_malariaAtlas = name_1, name2_malariaAtlas = name_2, iso3c_malariaAtlas = iso)
geodata_admin2_data_select <- geodata_adm2_africa %>%
  select(GID_0, NAME_1, NAME_2, geometry) %>%
  rename(name1_geodata = NAME_1, name2_geodata = NAME_2, iso3c_geodata = GID_0)

# Merge based on the geometry ID
merged_data <- sf::st_join(sf_adm2_africa_select, geodata_admin2_data_select, by = "geometry", left = TRUE)

View(merged_data %>% filter(is.na(name2_geodata)))
View(merged_data %>% filter(is.na(name2_malariaAtlas)))

#------------- FIND ADMIN2 REGIONS WHERE GEOMETRY OVERLAPS BUT NAME DOESNT MATCH
# Filter cases where the geometry matches but names do not
geometry_overlap_name_mismatch <- merged_data %>%
  filter(!is.na(name2_malariaAtlas) & !is.na(name2_geodata) & name2_malariaAtlas != name2_geodata) %>%
  select(iso3c_malariaAtlas, name2_malariaAtlas, iso3c_geodata, name2_geodata, geometry)

iso_mismatch <- geometry_overlap_name_mismatch %>%
  filter(iso3c_malariaAtlas != iso3c_geodata) %>%
  select(iso3c_malariaAtlas, name2_malariaAtlas, iso3c_geodata, name2_geodata, geometry)

#------------- FIND ADMIN2 REGIONS WHERE GEOMETRY WAS NOT FOUND
# Find geometries that exist in one dataset but not the other
missing_in_geodata <- sf_adm2_africa_select %>%
  filter(!geometry %in% geodata_admin2_data_select$geometry) %>%
  select(geometry, name2_malariaAtlas, iso3c_malariaAtlas) %>%
  mutate(dataset = "MalariaAtlas")

missing_in_malariaAtlas <- geodata_admin2_data_select %>%
  filter(!geometry %in% sf_adm2_africa_select$geometry) %>%
  select(geometry, name2_geodata, iso3c_geodata) %>%
  mutate(dataset = "Geodata")

geometry_not_found <- bind_rows(
  missing_in_geodata %>% rename(name2 = name2_malariaAtlas),
  missing_in_malariaAtlas %>% rename(name2 = name2_geodata)
)

cat("Admin2 Regions in Malaria Atlas that dont match with geodata:\n", dim(geometry_not_found %>% filter(dataset == "MalariaAtlas"))[1], "\n")
cat("Admin2 Regions in geodata that dont match with Malaria Atlas:\n", dim(geometry_not_found %>% filter(dataset == "MalariaAtlas"))[1], "\n")
