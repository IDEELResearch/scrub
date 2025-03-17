library(tidyverse)
library(here)

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

# consistent set of names we want
column_names <- get_column_names_for_clean()

# make our full bind across
full_bind <- rbind(
  clean_geoff %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)), 
  # clean_wwarn %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_pf7k %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character)),
  clean_who %>% select(all_of(column_names)) %>% mutate(across(everything(), as.character))
)

# TODO: Cecile:

# 1. approach could be to add extra columns to help with deduplication
# e.g. if collection start and end we think is too specific then perhaps year
full_bind$collection_year_start <- lubridate::year(full_bind$collection_start)
full_bind$collection_year_end <- lubridate::year(full_bind$collection_end)

# 2. there is a little helper function to show what gets removed by a specific set of columns
rm_ex <- rows_removed_by_distinct(full_bind, url, country_name, lat, lon, gene, 
                                  variant_string, total_num, collection_year_start, 
                                  collection_year_end)

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

