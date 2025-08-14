library(stringr)
library(purrr)
library(dplyr)
library(here)

# k13 mdr1 and crt mutants
wwarn_data <- readRDS(here::here("analysis", "data-derived", 
                           "wwarn_res.rds"))
ref_als <- read.csv(here::here("analysis", "data-raw", 
                                 "k13_ref_protein_codon_dictionary.csv"))

# TODO: extend this so this now works for mdr1 and crt

# clean all the data to fix issues in the wwarn data structure and make compatible with stave
# TODO: check if there is anything that no longer works with pd data also included
set.seed(1)
wwarn <- wwarn_data %>%
  dplyr::group_by(pmid) %>%
  dplyr::mutate(study_uid = paste0("wwarn","_",sample(10000:99999,1, replace = FALSE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(publication_status = "peer-reviewed", 
                database = "wwarn",
                # TODO: figure out how to pull in the real publication years from the pmid -- non-urgent
                publication_year = 1000) %>% # placeholder but this is so stave doesn't error
  dplyr::rowwise() %>%
  dplyr::mutate(source = str_split(source, " ")[[1]][1]) %>%
  dplyr::mutate(site_fixed = if_else(is.na(site), # sum(is.na(wwarn$site_fixed))
                                     gsub(" ","",admin_1), 
                                     gsub(" ", "", site)),
                country_name = countrycode::countrycode(iso3c, # fix up top
                                                        origin = 'iso3c', 
                                                        destination = 'country.name')) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(study_end_year = if_else(is.na(study_end_year),
                                         as.numeric(year), study_end_year)) %>%
  dplyr::mutate(study_start_year = if_else(is.na(study_start_year),
                                         as.numeric(year), study_start_year)) %>%
  # fix the years on this study
  # TODO: find the correct values for the years
  # dplyr::mutate(study_start_year = if_else(study_start_year == 1, 2001, study_start_year))
  dplyr::filter(study_end_year != 3) %>% 
  dplyr::rename(first_author_surname = source) %>%
  dplyr::mutate(authors = stringr::word(first_author_surname,1)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(across(c(study_uid, site, study_end_year))) %>%
  dplyr::mutate(survey_id = paste0(study_uid,"_",authors,"_",site_fixed,"_",study_end_year)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(survey_id = gsub("[^a-zA-Z0-9_]", "", survey_id)) %>% # all the things throwing errors in stave
  dplyr::rowwise() %>%
  dplyr::mutate(survey_id = iconv(survey_id, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(long = as.numeric(long),
                lat = as.numeric(lat)) %>%
  dplyr::rename(study_ID = study_uid) %>%
  dplyr::mutate(study_ID = paste0(study_ID, "_", authors),
                study_name = study_ID,
                study_type = "peer_reviewed") %>%
  dplyr::mutate(study_ID = iconv(study_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(study_ID =  gsub("[^a-zA-Z0-9_]", "", study_ID)) %>% 
  dplyr::mutate(continent = countrycode::countrycode(iso3c, origin = "iso3c", destination = "continent")) %>%
  dplyr::filter(continent == "Africa") |>
  dplyr::rowwise() 

wwarn_k13 <- wwarn %>% filter(gene == "k13")
wwarn_crt <- wwarn %>% filter(gene == "crt")
wwarn_mdr1 <- wwarn %>% filter(gene == "mdr1")

studies_all_WT_k13 <- wwarn_k13 %>%
  group_by(pmid) %>%
  summarise(all_WT = all(gene_mut == "k13:WT")) %>%
  filter(all_WT) %>%
  pull(pmid)

# no all WT studies for crt or mdr1 so no need to impute
studies_all_WT_crt <- wwarn_crt %>%
  group_by(pmid) %>%
  summarise(all_WT = all(gene_mut == "crt:WT")) %>%
  filter(all_WT) %>%
  pull(pmid)

studies_all_WT_mdr1 <- wwarn_mdr1 %>%
  group_by(pmid) %>%
  summarise(all_WT = all(gene_mut == "mdr1:WT")) %>%
  filter(all_WT) %>%
  pull(pmid)

# separate out the mutant data frame

mutation_key_path <- here::here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")
mutation_key <- read.csv(mutation_key_path)

validated_path <- here::here("analysis", "data-raw", "mutation_dictionary.csv")
validated_markers <- read.csv(validated_path) |> 
  dplyr::filter(gene == "k13") |>
  dplyr::filter(annotation != "invalidated")
codons <- as.numeric(gsub("[^\\d]+", "", validated_markers$mut, perl = TRUE))
reference_validated <- mutation_key |>
  dplyr::filter(CODON %in% codons) |>
  dplyr::arrange(CODON) |>
  dplyr::mutate(gene_mut = paste0(PROTEIN, ":", CODON, ":", REF))

# to add the pmids that need manual investigation
imputed <- data.frame()
studies <- wwarn_k13 %>% split(.$pmid)

# TODO: update the function so it actually updates wt_only w/i the function

# impute wt mutatations for k13
for(i in 1:length(studies)) {
  imputed <- rbind(imputed, impute_study(studies[[i]]))
}

# this is currently missing the 28 PMIDs with only WT from that study 
# this needs manually checking 

wwarn_k13 <- imputed |>
  dplyr::select(-c(siid))

# TODO: manually check the 28 excluded PMIDs and fix and add these into the datasets

# fix column names so its stave compatible
# rename the columns so that they make sense
column_names <- get_column_names_for_clean()

# CLEAN K13
# first add in extra information needed within column_names
wwarn_k13_edit <- wwarn_k13 %>% 
  dplyr::mutate(survey_ID = paste0(study_ID,"_wwarn_", admin_1,"_",year)) |> # make a survey ID 
  ungroup() %>% 
  dplyr::mutate(database = "wwarn",
                iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c"),
                study_ID = gsub("-", "_", study_ID),
                survey_ID = gsub("-", "_", gsub(" ", "", survey_ID))) %>% 
  dplyr::mutate(survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(continent = countrycode::countrycode(sourcevar = iso3c, "iso3c", "continent"),
                spatial_notes = "lat/lon reported in wwarn_k13",
                time_notes = "Automated midpoint") %>% 
  dplyr::rename(site_name = admin_1, 
                lon = long,
                variant_string = gene_mut,
                variant_num = x,
                total_num = n)

# sort out collecton dates
convstart <- adjust_invalid_date(unique(wwarn_k13_edit$year), is_start = TRUE)
wwarn_k13_edit$collection_start <- convstart[match(wwarn_k13_edit$year, unique(wwarn_k13_edit$year))]
convend <- adjust_invalid_date(unique(wwarn_k13_edit$year), is_start = FALSE)
wwarn_k13_edit$collection_end <- convend[match(wwarn_k13_edit$year, unique(wwarn_k13_edit$year))]

# Compute collection day
wwarn_k13_edit <- wwarn_k13_edit %>% rowwise() %>% 
  mutate(collection_day = median(c(collection_start, collection_end))) %>% 
  ungroup()

# Grab just the columns we need for pairing with wwarn_k13 etc
wwarn_k13_clean <- wwarn_k13_edit %>% 
  select(all_of(column_names)) %>%
# TODO: fix variant string k13:470:X - for now, just exclude them
  dplyr::filter(variant_string != "k13:470:X")

# Last extra cleans to align with STAVE
wwarn_k13_clean$survey_ID <- gsub("'", "", wwarn_k13_clean$survey_ID)

# CLEAN PD
wwarn_pd <- rbind(wwarn_crt, wwarn_mdr1)
names(wwarn_pd)

# first add in extra information needed within column_names
wwarn_pd_edit <- wwarn_pd %>% 
  # dplyr::mutate(survey_ID = paste0(study_ID,"_wwarn_", admin_1,"_",year)) |> # make a survey ID 
  ungroup() %>% 
  dplyr::rename(survey_ID = survey_id) %>%
  dplyr::mutate(database = "wwarn",
                iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c"),
                study_ID = gsub("-", "_", study_ID),
                survey_ID = gsub("-", "_", gsub(" ", "", survey_ID))) %>%
  dplyr::mutate(survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(continent = countrycode::countrycode(sourcevar = iso3c, "iso3c", "continent"),
                spatial_notes = "lat/lon reported in wwarn_pd",
                time_notes = "Automated midpoint") %>% 
  dplyr::rename(site_name = admin_1, 
                lon = long,
                variant_string = gene_mut,
                variant_num = x,
                total_num = n) %>%
  # need to rename the WT so that stave doesn't error
  dplyr::mutate(variant_string = if_else(variant_string == "mdr1:WT", "mdr1:86:N", variant_string),
                variant_string = if_else(variant_string == "crt:WT", "crt:76:K", variant_string))

# sort out collecton dates
convstart <- adjust_invalid_date(unique(wwarn_pd_edit$year), is_start = TRUE)
wwarn_pd_edit$collection_start <- convstart[match(wwarn_pd_edit$year, unique(wwarn_pd_edit$year))]
convend <- adjust_invalid_date(unique(wwarn_pd_edit$year), is_start = FALSE)
wwarn_pd_edit$collection_end <- convend[match(wwarn_pd_edit$year, unique(wwarn_pd_edit$year))]

# Compute collection day
wwarn_pd_edit <- wwarn_pd_edit %>% rowwise() %>% 
  mutate(collection_day = median(c(collection_start, collection_end))) %>% 
  ungroup()

# Grab just the columns we need for pairing with wwarn_pd etc
wwarn_pd_clean <- wwarn_pd_edit %>% 
  select(all_of(column_names)) # all variant_string and other variables look good

# Last extra cleans to align with STAVE
wwarn_pd_clean$survey_ID <- gsub("'", "", wwarn_pd_clean$survey_ID)

wwarn_clean <- rbind(wwarn_k13_clean, wwarn_pd_clean)
saveRDS(wwarn_clean, "analysis/data-derived/wwarn_clean.rds")
