# remove-wwarn.R
#
# Author: Gina Cuomo-Dannenburg
# Date: 2026-03-16
#
# Inputs:
# - "analysis/data-raw/WWARN_K13_database_04-12-2023.xls"
# - "analysis/data-raw/WWARN_partnerdrug_database_04-12-2023.xls"
# - "analysis/covidence/webofscience-march26.csv" 
# - "analysis/covidence/ovid-medline-march26.ris"
# - "analysis/covidence/pubmed-march26.nbib"

#
# Outputs:
# - covidence/wos-pmids.csv
# - covidence/pubmed-pmids.csv
# - covidence/medline-pmids.csv
#
# Purpose:
# Reads in the RIS/nbib files from each of the databases for import into covidence.
# It then identifies and removes studies that are already captured in WWARN database.
#
# ------------------------------------------------------------------

# libraries
library(tidyverse)
library(revtools) # tools for working with references

# read in files
# WWARN
wwarn_k13 <- read_excel("analysis/data-raw/WWARN_K13_database_04-12-2023.xls")
wwarn_pd <-  read_excel("analysis/data-raw/WWARN_partnerdrug_database_04-12-2023.xls")

# 1. Read in files --------------------------------------------------------
# read in from csv -- these work
wos_bib <- revtools::read_bibliography("analysis/covidence/webofscience-march26.csv") 
medline_bib <- revtools::read_bibliography("analysis/covidence/ovid-medline-march26.ris")

# read in pubmed
pubmed_lines <- readLines("analysis/covidence/pubmed-march26.nbib", warn = FALSE)

# remove empty lines
pubmed_lines <- pubmed_lines[nchar(pubmed_lines) > 0]

# remove lines shorter than 3 characters (cannot contain a tag) -- checked that this doesn't remove any references
pubmed_lines <- pubmed_lines[nchar(pubmed_lines) >= 3]

writeLines(pubmed_lines, "analysis/covidence/pubmed-march26_clean.nbib")

pubmed_bib <- bibliometrix::convert2df("analysis/covidence/pubmed-march26_clean.nbib",
                                       dbsource = "pubmed",
                                       format = "pubmed")

# summarise # of papers before the date filtering
nrow(wos_bib) # 247
nrow(medline_bib) # 369
nrow(pubmed_bib) # 868

# 2. Filter for publication dates in the range so all databases ar --------
# 2a. MEDLINE -------------------------------------------------------------

# year = year; y2 is the publication date formatted a little strangely and needs fixing
medline_bib$pub_date <- sub("//", "", x = medline_bib$y2)
medline_bib$pub_date <- gsub('^(.{4})(.*)$', '\\1-\\2', medline_bib$pub_date)
medline_bib$pub_date <- gsub('^(.{7})(.*)$', '\\1-\\2', medline_bib$pub_date)
medline_bib$pub_date <- as.Date(medline_bib$pub_date)


# filter by dates
sum(is.na(medline_bib$year))
medline <- medline_bib %>% split(.$year)

# only 2014 and 2025 need looking at
# fix publication date for only study we're interested in
medline$`2014`$pub_date[which(is.na(medline$`2014`$pub_date))] <- as.Date("2014-08-14") # 25056012
medline$`2014` <- medline$`2014` %>%
  dplyr::filter(pub_date > as.Date("2014-09-25"))

# find missing publication dates for 2025
fix_pmids <- medline$`2025`$pmid[which(is.na(medline$`2025`$pub_date))]
fix_index <- which(is.na(medline$`2025`$pub_date))
medline$`2025`$pub_date[fix_index[1]] <- as.Date("2025-10-15") # 40241390 
medline$`2025`$pub_date[fix_index[2]] <- as.Date("2025-09-02") # 40219818 
medline$`2025`$pub_date[fix_index[3]] <- as.Date("2025-07-30") # 39928049 
medline$`2025`$pub_date[fix_index[4]] <- as.Date("2025-02-04") # 39367758 
medline$`2025`$pub_date[fix_index[5]] <- as.Date("2025-07-30") # 39936828 
medline$`2025`$pub_date[fix_index[6]] <- as.Date("2025-07-11") # 40138574 
medline$`2025`$pub_date[fix_index[7]] <- as.Date("2025-07-01") # 40343740

medline$`2025` <- medline$`2025` %>%
  dplyr::filter(pub_date < as.Date("2025-07-09"))

# with dates filtered out
medline_bib <- do.call(rbind, medline)


# 2b. PubMed --------------------------------------------------------------

pubmed_bib$pub_date <- as.Date(pubmed_bib$EDAT)
pubmed_bib <- pubmed_bib %>%
  dplyr::filter(pub_date %in% seq.Date(from = as.Date("2014-07-09"),
                                       to = as.Date("2025-09-25"), by = 1))

# 2c. Check WoS -----------------------------------------------------------
# from manual check need to filter WoS also
wos <- wos_bib %>% split(.$publication_year)

wos$`2014`$publication_date # all after September so fine
wos$`2025`$publication_date # one study missing and one in October

wos$`2025`$publication_date[which(wos$`2025`$publication_date == "")] <- "2024-11-21"
wos$`2025` <- wos$`2025` %>%
  dplyr::filter(publication_date != "OCT 15") # only need to remove one study

wos_bib <- do.call(rbind, wos)

nrow(medline_bib) # 326 -- removed 43
nrow(pubmed_bib) # 799 -- removed 69  
nrow(wos_bib) # 246 -- removed 1

wos_bib$pmid <- wos_bib$pubmed_id %>% as.integer()
medline_bib$pmid <- medline_bib$ID %>% as.integer()
pubmed_bib$pmid <- as.integer(pubmed_bib$PMID)

# 3. Identify WWARN PMIDs and filter dfs ----------------------------------------

wwarn_pmid <- c(wwarn_k13$pubMedId, wwarn_pd$notes) %>%
  unique()
pmid <- data.frame(pmid = wwarn_pmid, wwarn = 1) |>
  dplyr::filter(is.na(pmid) == FALSE) |>
  dplyr::filter((pmid %in% c("500000000", "99999999", "80002535", "80002424", "80002356", "80000390", "80000305")) == FALSE) |>
  dplyr::distinct()
wwarn <- pmid |> dplyr::pull(pmid) %>% as.integer

# deduplicated dataframes
wos_deduplicated <- wos_bib %>%
  dplyr::filter(!(pmid %in% wwarn)) %>%
  dplyr::filter(!(is.na(pmid))) # checked this paper -- published 2026
medline_deduplicated <- medline_bib %>%
  dplyr::filter(!(pmid %in% wwarn)) 
pubmed_deduplicated <- pubmed_bib %>%
  dplyr::filter(!(pmid %in% wwarn)) 

# before and after deduplications
nrow(wos_bib)
nrow(wos_deduplicated)
nrow(medline_bib)
nrow(medline_deduplicated)
nrow(pubmed_bib)
nrow(pubmed_deduplicated)

# 4. Output files to convert to RIS for Covidence -------------------------
wos_pmids <- wos_deduplicated$pmid
medline_pmids <- medline_deduplicated$pmid
pubmed_pmids <- pubmed_deduplicated$pmid

write.csv(wos_pmids, "analysis/covidence/wos_pmids.csv")
write.csv(medline_pmids, "analysis/covidence/medline_pmids.csv")
write.csv(pubmed_pmids, "analysis/covidence/pubmed_pmids.csv")
