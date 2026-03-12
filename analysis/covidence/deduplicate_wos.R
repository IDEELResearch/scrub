# deduplicate_ris.R
#
# Author: Gina Cuomo-Dannenburg
# Date: 2026-03-12
#
# Inputs:
# - data-derived/pf8_STAVE.rds
# - data-derived/WHO_STAVE.rds
# - data-derived/WWARN_STAVE.rds
# - covidence/webofscience.xls
# - covidence/webofscience.ris
#
# Outputs:
# - covidence/webofscience_deduplicated.ris
#
# Purpose:
# Reads in PMID_to_k13_coverage.csv, which flags for each study (PMID) which WHO
# k13 loci were sequenced. This is used to generate the corresponding wild-type
# variant for each study, i.e. the expected AA code if no non-synonymous
# mutations were observed at any of these sequenced loci. This information is
# saved in a modified version of the original data frame in the data-derived
# folder.
#
# ------------------------------------------------------------------

# load libraries
# install.packages(c("revtools", "dplyr"))
library(revtools)
library(dplyr)
library(tidyverse)
library(STAVE)

# identify the PMIDs that are already included in other data sources to then remove these from the import
# read in the RDS objects of the databases

pf8 <- readRDS("analysis/data-derived/pf8_STAVE.rds")
who <- readRDS("analysis/data-derived/WHO_STAVE.rds")
wwarn <- readRDS("analysis/data-derived/WWARN_STAVE.rds")

# get pmids
pf8_pmid <- pf8$get_studies() %>% as.data.frame() %>% dplyr::pull(PMID) %>% unique()
who_pmid <- who$get_studies() %>% as.data.frame() %>% dplyr::pull(PMID) %>% unique()
wwarn_pmid <- wwarn$get_studies() %>% as.data.frame() %>% dplyr::pull(PMID) %>% unique()

pmids <- c(pf8_pmid,
           who_pmid,
           wwarn_pmid) %>% unique()

# identify the papers in the xls that need to be removed from the ris
wos <- readxl::read_xls("analysis/covidence/webofscience.xls")
wos_ris <- revtools::read_bibliography("analysis/covidence/webofscience.ris")

wos$PMID <- wos$`Pubmed Id`
which(wos$PMID %in% pmids) %>% length()
View(wos)

# wos$UT == wos_ris$accession -> match by this
# identify the accessions to remove
wos_accessions <- wos %>%
  dplyr::mutate(exclude = dplyr::if_else(PMID %in% pmids, TRUE, FALSE)) %>%
  dplyr::filter(exclude == TRUE) %>%
  dplyr::pull(`UT (Unique WOS ID)`)

# remove the entries that are already in other databases
wos_clean <- wos_ris %>%
  dplyr::filter((accession %in% wos_accessions) == FALSE)

# check that removed the correct number of studies
(nrow(wos_ris) - nrow(wos_clean)) == length(wos_accessions)

write_bibliography(wos_clean, file = "analysis/covidence/webofscience_deduplicated.ris")

