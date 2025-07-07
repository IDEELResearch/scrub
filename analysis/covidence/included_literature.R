# load libraries

library(tidyverse)
library(readxl)

# read in data
# WWARN
wwarn_k13 <- read_excel("analysis/data-raw/WWARN_K13_database_04-12-2033.xls")
wwarn_pd <-  read_excel("analysis/data-raw/WWARN_partnerdrug_database_04-12-2023.xls")
wwarn_pmid <- c(wwarn_k13$pubMedId, wwarn_pd$notes) %>%
  unique()
pmid <- data.frame(pmid = wwarn_pmid, wwarn = 1) |>
  dplyr::filter(is.na(pmid) == FALSE) |>
  dplyr::filter((pmid %in% c(500000000, 99999999, 80002535, 80002424, 80002356, 80000390, 80000305)) == FALSE) |>
  dplyr::distinct()
wwarn <- pmid |> pull(pmid)

# medline
lines <- readLines("analysis/covidence/medline.txt")

# Find the lines that say "Unique Identifier"
uid_label_rows <- grep("^Unique Identifier\\s*$", lines)

# The actual PMIDs are the next line after each label
pmid_lines <- lines[uid_label_idx + 1]

# Trim leading/trailing whitespace
medline <- trimws(pmid_lines)

medline_pmids <- medline[medline %in% wwarn == FALSE]
write.csv(medline_pmids, file = "analysis/covidence/medline_pmids.csv")

# pubmed
# Read the nbib file
pubmed <- readLines("analysis/covidence/pubmed-hackathon-july.nbib")

# Extract lines that start with PMID
pubmed <- grep("^PMID- ", pubmed, value = TRUE)

# Get the numeric PMIDs
pubmed <- sub("^PMID- ", "", pubmed)

pubmed_pmids <- pubmed[pubmed %in% wwarn == FALSE]
write.csv(pubmed_pmids, file = "analysis/covidence/pubmed_pmids.csv")
