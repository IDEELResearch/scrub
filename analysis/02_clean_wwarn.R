library(stringr)
library(purrr)
library(dplyr)
library(here)

# k13 mdr1 and crt mutants
wwarn_data <- readRDS(here::here("analysis", "data-derived", 
                           "wwarn_res.rds"))
ref_als <- read.csv(here::here("analysis", "data-raw", 
                                 "k13_ref_protein_codon_dictionary.csv"))

# clean all the data to fix issues in the wwarn data structure and make compatible with stave
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

# identify studies where only WT is reported
studies_all_WT_k13 <- wwarn_k13 %>%
  group_by(pmid) %>%
  summarise(all_WT = all(gene_mut == "k13:WT")) %>%
  filter(all_WT) %>%
  pull(pmid)

# no all WT studies for crt or mdr1 so will impute more simply at the end of the script
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

# Add collection start and end for wwarn_72050_Ibrahim
wwarn_pd_edit <- wwarn_pd_edit %>%
  mutate(
    # For this one study, override start and end as proper Date objects
    collection_start = if_else(
      study_ID == "wwarn_72050_Ibrahim",
      as.Date("2003-10-01"),
      collection_start
    ),
    collection_end = if_else(
      study_ID == "wwarn_72050_Ibrahim",
      as.Date("2006-10-31"),
      collection_end
    )
  ) %>%
  mutate(
    # midpoint between start and end
    midpoint = collection_start + (collection_end - collection_start) / 2,
    # set collection_day only for this study
    collection_day = if_else(
      study_ID == "wwarn_72050_Ibrahim",
      midpoint,
      collection_day
    )
  ) %>%
  select(-midpoint)

# Grab just the columns we need for pairing with wwarn_pd etc
wwarn_pd_clean <- wwarn_pd_edit %>% 
  select(all_of(column_names)) # all variant_string and other variables look good

# Last extra cleans to align with STAVE
wwarn_pd_clean$survey_ID <- gsub("'", "", wwarn_pd_clean$survey_ID)

# 2. Impute the pd WT and mutants
# now impute the mutations in mdr1 and crt 
mdr1 <- wwarn_pd_clean %>% filter(gene == "mdr1") %>% distinct()
crt <- wwarn_pd_clean %>% filter(gene == "crt") %>% distinct()
nrow(mdr1) + nrow(crt) == nrow((wwarn_pd_clean) %>% distinct())

# make a dataset with the number of variants reported per survey
mdr1 <- mdr1 %>% 
  mutate(variant_string = standardise_mdr1_86(variant_string)) %>% 
  filter(variant_string != "mdr1:86:F") %>%
  dplyr::group_by(survey_ID, collection_day) %>%
  mutate(n = length(unique(total_num))) %>%
  filter(n == 1) %>% # exclude surveys with different Ns
  ungroup() %>%
  group_by(survey_ID, collection_day, total_num, variant_string) %>%
  mutate(variant_num = sum(variant_num),
         total_num = unique(total_num)) %>%
  ungroup() %>%
  distinct(survey_ID, collection_day, total_num, variant_string,.keep_all = TRUE) %>% 
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(num_variants = n())

# split by the number of variants in the survey
mdr1_1 <- mdr1 %>% filter(num_variants == 1) # nrow == 0
mdr1_2 <- mdr1 %>% filter(num_variants == 2)
mdr1_3 <- mdr1 %>% filter(num_variants == 3)

# now clean those with 2 mutations
# 3 possible iterations of 2 mutations
# add a column that indicates which of the three classes it falls into so that we can treat them separately
mdr1_2_groups <- mdr1_2 %>%
  group_by(survey_ID, collection_day) %>%
  mutate(
    pair = unique(variant_string) |> str_c(collapse = "+"),
    pair_simple = case_when(
      pair %in% c("mdr1:86:N+mdr1:86:Y", "mdr1:86:Y+mdr1:86:N")   ~ "N+Y",
      pair %in% c("mdr1:86:N+mdr1:86:N/Y", "mdr1:86:N/Y+mdr1:86:N") ~ "N+N/Y",
      pair %in% c("mdr1:86:Y+mdr1:86:N/Y", "mdr1:86:N/Y+mdr1:86:Y") ~ "Y+N/Y",
      TRUE                                                          ~ pair
    )
  ) %>%
  ungroup() 

mdr1_2_groups <-  mdr1_2_groups %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  mutate(correct = if_else(xs == total_num, "yes", "no"))

table(mdr1_2_groups$correct) 
# all are correct now

# split the dataframe with two unique variants into what combination of variants - M ~ mixed
mdr1_2_NY <- mdr1_2_groups %>% filter(pair_simple == "N+Y") %>% filter(xs <= total_num)
mdr1_2_NM <- mdr1_2_groups %>% filter(pair_simple == "N+N/Y") %>% filter(xs <= total_num)
mdr1_2_YM <- mdr1_2_groups %>% filter(pair_simple == "Y+N/Y") %>% filter(xs <= total_num)

# in WWARN, we only have the case of NY. this means both mutant and wt are reported
# where all are correct this also means that the data is fully imputed already
# can simply use mdr1_2

# now check when all three are reported 
(mdr1_3 %>% 
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs > total_num) %>% 
  nrow()) == 0 

# just to ensure that the names are correct but no imputation waws required
mdr1 <- rbind(mdr1_1,
              mdr1_2, 
              mdr1_3) %>%
  select(names(wwarn_k13_clean)) %>%
  mutate(prev = variant_num / total_num)

# make a dataset with the number of variants per survey
unique(crt$variant_string) # no variants we want to filter out

crt <- crt %>%
  dplyr::group_by(survey_ID, collection_day) %>%
  mutate(n = length(unique(total_num))) %>%
  filter(n == 1) %>% # exclude surveys with different Ns
  ungroup() %>%
  group_by(survey_ID, collection_day, total_num, variant_string) %>%
  mutate(variant_num = sum(variant_num),
         total_num = unique(total_num)) %>%
  ungroup() %>%
  distinct(survey_ID, collection_day, total_num, variant_string,.keep_all = TRUE) %>% 
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(num_variants = n())

table(crt$num_variants)

# split by the number of variants in the survey
crt_1 <- crt %>% filter(num_variants == 1) 
crt_2 <- crt %>% filter(num_variants == 2)
crt_3 <- crt %>% filter(num_variants == 3)

# look at the 3 surveys with only one mutant
# manually fix this ID: "wwarn_16216_Duah_Hohoe_2004" -- other paper is correct
crt_1 <- crt_1 %>%
  split(.$study_ID)

# look back at PMID and manually fix
crt_1$wwarn_16216_Duah <- crt_1$wwarn_16216_Duah %>%
  bind_rows(crt_1$wwarn_16216_Duah, crt_1$wwarn_16216_Duah)

# assume that K76T means mixed -- all values given in percentages and then converted
# K76 20 = 26/130
# T76 45 = 58 / 130
# K76T 35 = 46 / 130
crt_1$wwarn_16216_Duah$variant_string <- c("crt:76:K",
                                           "crt:76:T",
                                           "crt:76:K/T")

crt_1$wwarn_16216_Duah$variant_num <- c(26, 58, 46)
crt_1$wwarn_16216_Duah$prev <- crt_1$wwarn_16216_Duah$variant_num / crt_1$wwarn_16216_Duah$total_num

(crt_1$wwarn_16216_Duah %>% 
    mutate(xs = sum(variant_num)) %>% 
    pull(xs)) == (unique(crt_1$wwarn_16216_Duah$total_num)) # checked this is correct now

crt_1 <- do.call(rbind, crt_1)

# now clean those with 2 mutations
# 3 possible iterations of 2 mutations
# add a column that indicates which of the three classes it falls into so that we can treat them separately
crt_2_groups <- crt_2 %>%
  group_by(survey_ID, collection_day) %>%
  mutate(
    pair = unique(variant_string) |> str_c(collapse = "+"),
    pair_simple = case_when(
      pair %in% c("crt:76:K+crt:76:T", "crt:76:T+crt:76:K")   ~ "K+T",
      pair %in% c("crt:76:K+crt:76:K/T", "crt:76:K/T+crt:76:K") ~ "K+K/T",
      pair %in% c("crt:76:T+crt:76:K/T", "crt:76:K/T+crt:76:T") ~ "T+K/T",
      TRUE                                                          ~ pair
    )
  ) %>%
  ungroup() 

crt_2_groups <-  crt_2_groups %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  mutate(correct = if_else(xs == total_num, "yes", "no"))

table(crt_2_groups$correct) 
# 10 not correct -- look at these and fix

# 3 studies to fix
crt_2_fix <- crt_2_groups %>% 
  filter(correct == "no") %>%
  split(.$study_ID)

# Olefongo is correct -- the difference are other mutations that we are not interested in.
# means sum(x) != n but these are not missing mixed mutants based on the paper

# simple typo -- from looking at the whole study and checking totals 76T should have one more variant
crt_2_fix$wwarn_62503_Afoakwah$variant_num[1] <- 37
sum(crt_2_fix$wwarn_62503_Afoakwah$variant_num) == crt_2_fix$wwarn_62503_Afoakwah$total_num

# both have 100 prev -- look at paper
# 13% 76T; 87% K76; N = 180
crt_2_fix$wwarn_79899_Ngassa$variant_num <- c(23, 180 - 23)
crt_2_fix$wwarn_79899_Ngassa$prev <- crt_2_fix$wwarn_79899_Ngassa$variant_num /
  crt_2_fix$wwarn_79899_Ngassa$total_num

crt_2_fix <- do.call(rbind, crt_2_fix) 

crt_2 <- rbind(crt_2_groups %>%
  filter(correct == "yes"), crt_2_fix)

# check that we haven't lost any surveys
length(unique(crt_2_groups$survey_ID)) == length(unique(crt_2$survey_ID))

crt_3 <- crt %>% filter(num_variants == 3)
crt_3_groups <- crt_3 %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  mutate(correct = if_else(xs == total_num, "yes", "no")) %>%
  arrange(correct, survey_ID) %>%
  split(.$correct)

# fix this study: xs > n
mensah <- crt_3_groups$no %>%
  filter(study_ID == "wwarn_96528_Mensah")
  
mensah <- mensah[1:2,]
mensah$variant_string <- c("crt:76:K", "crt:76:T")
mensah$variant_num <- c(59, 15)
mensah$total_num <- c(74, 74)
mensah$prev <- mensah$variant_num / mensah$total_num

bonizzoni <- crt_3_groups$no %>%
  filter(study_ID == "wwarn_43882_Bonizzoni")
bonizzoni$variant_num[1] <- 22
bonizzoni$variant_num[6] <- 24

bonizzoni <- bonizzoni %>%
  mutate(prev = variant_num / total_num) %>%
  group_by(survey_ID) %>%
  mutate(xs = sum(variant_num)) 

# check we haven't lost any surveys
length(unique(crt_3$survey_ID)) == (length(unique(crt_3_groups$yes$survey_ID)) + 
                                      length(unique(crt_3_groups$no$survey_ID)))

crt_3 <- do.call(rbind, crt_3_groups)

# just to ensure that the names are correct but no imputation waws required
crt <- rbind(crt_1,
              crt_2, 
              crt_3) %>%
  select(names(wwarn_k13_clean)) %>%
  mutate(prev = variant_num / total_num)

wwarn_pd_clean <- rbind(mdr1, crt) %>% select(names(wwarn_k13_clean))

wwarn_clean <- rbind(wwarn_k13_clean, wwarn_pd_clean)
saveRDS(wwarn_clean, "analysis/data-derived/wwarn_clean.rds")
