library(tidyverse)
library(here)

# Read each file if it exists
safe_read <- function(path) {
  if (file.exists(path)) {
    clean <- readRDS(path) %>% 
      as.data.frame() 
    if ("study_type" %in% names(clean)) {
      clean <- clean %>% 
        dplyr::filter(study_type != "unpublished")
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

# deduplication
full_bind <- deduplicate(full_bind)

# save ready to go to stave
saveRDS(full_bind, here("analysis/data-derived/final_data.rds"))
