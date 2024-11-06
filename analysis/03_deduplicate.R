library(tidyverse)
library(here)

<<<<<<< HEAD
# read data sets
clean_geoff <- readRDS(here("analysis", "data-derived", "02_clean_geoffs_output_table.rds"))
clean_wwarn <- readRDS(here("analysis", "data-derived", "wwarn_res_df.rds"))

## deduplicate by pmid
clean_geoff <- clean_geoff %>%
  dplyr::filter(is.na(pmid) == FALSE)
clean_wwarn <- clean_wwarn %>%
  dplyr::rowwise() %>%
  dplyr::filter(is.na(pmid) == FALSE)

geoff_pmid <- unique(clean_geoff$pmid)
wwarn_pmid <- unique(clean_wwarn$pmid)

geoff_pmid[which(geoff_pmid %in% wwarn_pmid)]
wwarn_pmid[which(wwarn_pmid %in% geoff_pmid)]
unique(geoff_pmid[which(geoff_pmid %in% wwarn_pmid)], wwarn_pmid[which(wwarn_pmid %in% geoff_pmid)])

clean_wwarn <- clean_wwarn %>%
  dplyr::filter(pmid %in% geoff_pmid == FALSE)

## add back the geoff pre-print data
preprint <- clean_geoff %>%
  dplyr::filter(source == "preprint") 

clean_geoff <- rbind(clean_geoff, preprint)  

# start deduplication process
# identify studies that may be duplicates from WWARN
geoff_sites <- unique(clean_geoff$site)
clean_geoff %>% dplyr::filter(study_key == "s0004_fola_eth_622i_valid") %>% View()
=======
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
>>>>>>> oj_clean_pre_geoff_in
