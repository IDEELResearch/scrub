library(tidyverse)

# ---------------------------------------------------- o
# 1. WHO database compile ----
# ---------------------------------------------------- o

# read in WHO study info
whodf <- readxl::read_xlsx("analysis/data-raw/WHO_res_database_02-01-2024.xlsx", sheet = 2)

# read in res WHO info
whores <- readxl::read_xlsx("analysis/data-raw/WHO_res_database_02-01-2024.xlsx", sheet = 3)

# combine by the res
who <- left_join(whores, whodf, by = "ID")

# sort names as wanted
who <- who %>%
  rename(admin_0 = COUNTRY_NAME) %>%
  rename(admin_1 = ADMIN2) %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  rename(lat = LATITUDE) %>%
  rename(long = LONGITUDE) %>%
  rename(study_start_year = YEAR_START) %>%
  mutate(year = NA) %>%
  mutate(study_end_year = NA) %>%
  rename(n = SAMPLE_SIZE) %>%
  mutate(prev = as.numeric(PROPORTION)/100) %>%
  mutate(x = round(n*prev)) %>%
  rename(gene = MM_TYPE) %>%
  rename(mut = GENOTYPE) %>%
  rename(url = CITATION_URL) %>%
  mutate(publication_year = 2023) %>% # TODO: Fix this to scrape from pubmeds where exists
  # happen to know that this is the Conrad paper early
  mutate(url = replace(url, url == "99999", "http://www.ncbi.nlm.nih.gov/pubmed/37611122")) %>%
  mutate(publication_year = replace(publication_year, url == "99999", 2023)) %>%
  mutate(pmid = gsub(".*pubmed/(\\d+).*|.*gov/(\\d+).*","\\1",url)) %>%
  mutate(pmid = replace(pmid, pmid == "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7029593/", "32070355")) %>%
  mutate(pmid = replace(pmid, pmid == "", NA)) %>%
  mutate(database = "WHO") %>%
  rename(source = DATA_SOURCE) %>%
  rename(site = SITE_NAME) %>%
  select(ID, publication_year,iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

# standardise gene names
who <- who %>%
  mutate(
    gene = replace(gene, gene == "Pfkelch13", "k13"),
    gene = replace(gene, gene == "Pfcrt", "crt"),
    gene = replace(gene, gene == "Pfmdr1", "mdr1"),
    gene = replace(gene, gene == "Pfplasmepsin 2-3", "pfpm23")
  )

# split these out into each marker type
whospl <- who %>% split(who$gene)

# 1. Handle each separately for ease - crt

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$crt %>%
  mutate(mut = replace(mut, mut == "Pfcrt", "crt_76T")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = (sum(prev))) %>%
  pull(p)

whospl$crt <- whospl$crt %>%
  mutate(mut = replace(mut, mut == "Pfcrt", "crt_76T")) %>%
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  summarise(n = sum(unique(n)), x = n - sum(x[mut == "WT"]), prev = x/n) %>%
  mutate(mut = "crt:76:T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup

# 2. Handle each separately for ease - k13_valid

# sanity check that data is correctly formatted from WHO
# i.e. that all samples for a uuid grouping sum to the
# number of samples tested,
whospl$k13 %>%
  mutate(mut = gsub("&", ",", mut, fixed = TRUE)) %>%
  mutate(mut = gsub("unspecified", "-", mut, fixed = TRUE)) %>%
  rowwise() %>% 
  mutate(mut = pf7k_format_k13_for_stave(mut)) %>%
  ungroup() %>% 
  group_by(across(c(-x, -n, -prev, -mut))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  summarise(p = sum(prev)) %>%
  pull(p) %>% table()

mdf2 <- whospl$k13 %>%
  mutate(mut = gsub("&", ",", mut, fixed = TRUE)) %>%
  mutate(mut = gsub("unspecified", "-", mut, fixed = TRUE)) %>%
  rowwise() %>% 
  mutate(mut = pf7k_format_k13_for_stave(mut)) %>%
  ungroup() 

# function to split duplicated entries
split_rows_by_and <- function(data, column) {
  # Use tidyr's separate_rows to split the specified column
  data %>%
    separate_rows(!!sym(column), sep = "\\s*&&\\s*")
}

# now clean the variant_string entries with && in and bind with the rest
res <- rbind(
  mdf2 %>% 
    filter(grepl("&&", mut)) %>% 
    split_rows_by_and("mut"),
  mdf2 %>%
    filter(!grepl("&&", mut))
)

# and lastly sort out the WT calls
mutation_key_path <- here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")
mutation_key <- read.csv(mutation_key_path)
indices_to_transform <- which(res$mut == "k13:WT")
# range <- "349-726" - this is the range noted in the original pf7k file
res$mut[indices_to_transform] <- gsub("K13", "k13", collapse_k13_range("k13:349:726"))

# bring it all back together
who_res_df <- rbind(whospl$crt, res)
saveRDS(who_res_df, here::here("analysis/data-derived/who_res.rds"))

