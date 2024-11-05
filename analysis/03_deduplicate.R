library(tidyverse)
library(here)

# Read each file if it exists
safe_read <- function(path) {
  if (file.exists(path)) {
    clean <- readRDS(path) %>% 
      as.data.frame() %>% 
      dplyr::filter(source != "unpublished")
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
comnms <- names(clean_wwarn)

# make our full bind across
full_bind <- rbind(
  # clean_geoff %>% select(all_of(comnms)) %>% mutate(across(everything(), as.character)), 
  clean_wwarn %>% select(all_of(comnms)) %>% mutate(across(everything(), as.character))
  # clean_pf7k %>% select(all_of(comnms)) %>% mutate(across(everything(), as.character)),
  # clean_who %>% select(all_of(comnms)) %>% mutate(across(everything(), as.character))
)

# deduplication
full_bind <- deduplicate(full_bind)

# save ready to go to stave
saveRRDS(full_bind, here("analysis/data-derived/final_data.rds"))
