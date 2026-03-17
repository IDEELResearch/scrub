
# ------------------------------------------------------------
# get abstracts from PMID for all non-pf8 studies
# (NB, if complains because too many queries, break PMID vector into multiple
# sections and stick back together)

library(rentrez)
library(xml2)

s_geoff <- readRDS(file = here("analysis", "data-derived", "geoff_STAVE.rds"))
s_wwarn <- readRDS(file = here("analysis", "data-derived", "WWARN_STAVE.rds"))
s_who <- readRDS(file = here("analysis", "data-derived", "WHO_STAVE.rds"))

# Input PMIDs
pmids <- c(s_geoff$get_studies()$PMID,
           s_wwarn$get_studies()$PMID,
           s_who$get_studies()$PMID)

pmids <- pmids[!is.na(pmids)]

# Fetch PubMed records as XML
xml_text_raw <- entrez_fetch(
  db = "pubmed",
  id = pmids,
  rettype = "xml",
  retmode = "xml"
)

# Parse XML
doc <- read_xml(xml_text_raw)

# Find all PubMed articles
articles <- xml_find_all(doc, ".//PubmedArticle")

# Extract fields
results <- lapply(articles, function(article) {
  
  pmid_node <- xml_find_first(article, ".//PMID")
  title_node <- xml_find_first(article, ".//ArticleTitle")
  abstract_nodes <- xml_find_all(article, ".//Abstract/AbstractText")
  journal_node <- xml_find_first(article, ".//Journal/Title")
  year_node <- xml_find_first(article, ".//PubDate/Year")
  
  pmid <- if (inherits(pmid_node, "xml_missing")) NA_character_ else xml_text(pmid_node)
  title <- if (inherits(title_node, "xml_missing")) NA_character_ else xml_text(title_node)
  journal <- if (inherits(journal_node, "xml_missing")) NA_character_ else xml_text(journal_node)
  year <- if (inherits(year_node, "xml_missing")) NA_character_ else xml_text(year_node)
  
  abstract <- if (length(abstract_nodes) == 0) {
    NA_character_
  } else {
    paste(xml_text(abstract_nodes), collapse = " ")
  }
  
  data.frame(
    PMID = pmid,
    Title = title,
    Abstract = abstract,
    Journal = journal,
    Year = year,
    stringsAsFactors = FALSE
  )
})

# Combine into one data frame
pubmed_data <- do.call(rbind, results)

# ------------------------------------------------------------
# search for WGS data 

# option to read in pre-computed from file
#pubmed <- readRDS("/Users/rverity/Desktop/pubmed.rds")

# Patterns related to whole genome sequencing
wgs_patterns <- c(
  "\\bwhole genomes? sequenc\\w*\\b",
  "\\bwhole-genome sequenc\\w*\\b",
  "\\bsequencing of the whole genome\\b",
  "\\bwhole genome data\\b",
  "\\bwgs\\b",
  "\\bgenome[- ]wide sequencing\\b"
)

# match abstracts against patterns
pubmed$WGS <- FALSE
for (i in 1:nrow(pubmed)) {
  pubmed$WGS[i] <- any(str_detect(tolower(pubmed$Abstract[[i]]), wgs_patterns))
}

# manually fix papers where cannot get the abstract automatically
pubmed$WGS[346] <- FALSE

pubmed |>
  filter(WGS) |>
  View()
