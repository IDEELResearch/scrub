library(tidyverse)
library(here)

# read data sets
clean_geoff <- readRDS(here("analysis", "data-derived", "02_clean_geoffs_simple.RDS"))
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
