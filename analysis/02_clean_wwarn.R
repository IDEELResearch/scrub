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
mdr1 <- wwarn_pd_clean %>% filter(gene == "mdr1")
crt <- wwarn_pd_clean %>% filter(gene == "crt")
nrow(mdr1) + nrow(crt) == nrow(wwarn_pd_clean)

# make a dataset with the manual fixes combined
mdr1 <- mdr1_original %>% 
  mutate(variant_string = standardise_mdr1_86(variant_string)) %>% 
  filter(variant_string != "mdr1:86:F") %>%
  dplyr::group_by(survey_ID, collection_day) %>%
  mutate(n = length(unique(total_num))) %>%
  filter(n == 1) %>% # exclude surveys with different Ns
  ungroup() %>%
  rbind(mdr1_fix) %>%
  group_by(survey_ID, collection_day, total_num, variant_string) %>%
  mutate(variant_num = sum(variant_num),
         total_num = unique(total_num)) %>%
  ungroup() %>%
  distinct(survey_ID, collection_day, total_num, variant_string,.keep_all = TRUE) %>% 
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(num_variants = n())

# split by the number of variants in the survey
mdr1_1 <- mdr1 %>% filter(num_variants == 1)
mdr1_2 <- mdr1 %>% filter(num_variants == 2)
mdr1_3 <- mdr1 %>% filter(num_variants == 3)

# start with cleaning the surveys where we only have one variant
mdr1_1 %>% group_by(variant_string) %>% summarise(n = n()) # all are mdr1:86:N or mdr1:86:Y only

# split them into the type of one variant we have
mdr1_1_n <- mdr1_1 %>% filter(variant_string == "mdr1:86:N")
mdr1_1_y <- mdr1_1 %>% filter(variant_string == "mdr1:86:Y")

# for each row in mdr1_1_n add an additional row that is the complement
complement <- NULL
for(i in 1:nrow(mdr1_1_n)) {
  df <- mdr1_1_n[i,]
  df$variant_num <- df$total_num - df$variant_num
  df$variant_string <- "mdr1:86:Y"
  df$prev <- df$variant_num / df$total_num
  complement <- rbind(complement, df)
}

# checked that sum(variants) == denom
mdr1_1_n <- rbind(mdr1_1_n, complement) %>%
  arrange(survey_ID, collection_day, total_num) %>%
  mutate(prev = variant_num/total_num)

# for each row in mdr1_1_y add an additional row that is the complement
complement <- NULL
for(i in 1:nrow(mdr1_1_y)) {
  df <- mdr1_1_y[i,]
  df$variant_num <- df$total_num - df$variant_num
  df$variant_string <- "mdr1:86:N"
  df$prev <- df$variant_num / df$total_num
  complement <- rbind(complement, df)
}

mdr1_1_y <- rbind(mdr1_1_y, complement) %>%
  arrange(survey_ID, collection_day, total_num) %>%
  mutate(prev = variant_num/total_num)

mdr1_1 <- rbind(mdr1_1_n, mdr1_1_y) %>%
  mutate(prev = variant_num / total_num) 

# checked that sum(variants) == denom
mdr1_1 %>%
  group_by(collection_day, survey_ID, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>%
  nrow() 
# [1] 0

# now clean those that report all three because this is another simple case
mdr1_3 <- mdr1 %>% 
  filter(num_variants == 3) %>%
  group_by(collection_day, survey_ID, total_num) %>%
  mutate(xs = sum(variant_num)) 

mdr1_3_correct <- mdr1_3 %>%
  filter(xs == total_num)

# TODO: one study that still needs fixing

mdr1_3 <- mdr1_3_correct

# now clean those with 2 mutations
# 3 possible iterations of 2 mutations
# add a column that indicates which of the three classes it falls into so that we can treat them separately
mdr1_2_groups <- mdr1_2 %>%
  group_by(nid) %>%
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

# these rows are already correct and not missing imputation -- leave as is 
mdr1_2_groups_correct <- mdr1_2_groups %>% filter(correct == "yes")

# these are the rows that need fixing -- fix based on the classification
mdr1_2_groups <- mdr1_2_groups %>% filter(correct == "no") 

# TODO: fix these studies where sum(x) > denom 
# I think these are data issues 
mdr1_2_fix <- mdr1_2_groups %>% filter(xs > total_num)

# split the dataframe with two unique variants into what combination of variants - M ~ mixed
mdr1_2_NY <- mdr1_2_groups %>% filter(pair_simple == "N+Y") %>% filter(xs <= total_num)
mdr1_2_NM <- mdr1_2_groups %>% filter(pair_simple == "N+N/Y") %>% filter(xs <= total_num)
mdr1_2_YM <- mdr1_2_groups %>% filter(pair_simple == "Y+N/Y") %>% filter(xs <= total_num)

# NM = 2 rows
mdr1_2_NM[3,] <- mdr1_2_NM[2,]
mdr1_2_NM[3,]$variant_num <- mdr1_2_NM[3,]$total_num - sum(mdr1_2_NM$variant_num[1:2])
mdr1_2_NM[3,]$variant_string <- "mdr1:86:Y"
# check that sum numerator = denominator now 
mdr1_2_NM %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow()

# now looking at N & Y
mdr1_2_NY <- mdr1_2_NY %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(sid = cur_group_id()) %>%
  arrange(sid)

mdr1_2_NY_impute <- NULL
for(i in 1:length(unique(mdr1_2_NY$sid))) {
  df <- mdr1_2_NY[c((2*i)-1, 2*i),]
  df[3,] <- df[2,]
  df[3,]$variant_string <- "mdr1:86:N/Y"
  df[3,]$variant_num <- df[3,]$total_num - df[3,]$xs
  df[3,]$prev <- df[3,]$variant_num / df[3,]$total_num
  mdr1_2_NY_impute <- rbind(df, mdr1_2_NY_impute)
}
mdr1_2_NY_impute %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow() # check these are all fixed

mdr1_2_NY <- mdr1_2_NY_impute

# now clean YM
mdr1_2_YM <- mdr1_2_YM %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(sid = cur_group_id()) %>%
  arrange(sid)

mdr1_2_YM_impute <- NULL
for(i in 1:length(unique(mdr1_2_YM$sid))) {
  df <- mdr1_2_YM[c((2*i)-1, 2*i),]
  df[3,] <- df[2,]
  df[3,]$variant_string <- "mdr1:86:N/Y"
  df[3,]$variant_num <- df[3,]$total_num - df[3,]$xs
  df[3,]$prev <- df[3,]$variant_num / df[3,]$total_num
  mdr1_2_YM_impute <- rbind(df, mdr1_2_YM_impute)
}

mdr1_2_YM_impute %>%
  group_by(survey_ID, collection_day, total_num) %>%
  mutate(xs = sum(variant_num)) %>%
  filter(xs != total_num) %>% nrow() # check these are all fixed

mdr1_2_YM <- mdr1_2_YM_impute

mdr1_2 <- rbind(mdr1_2_NM,
                mdr1_2_NY,
                mdr1_2_YM)

mdr1 <- rbind(mdr1_1,
              mdr1_2, 
              mdr1_3) %>%
  select(names(master_table_simplified)) %>%
  mutate(prev = variant_num / total_num)




wwarn_clean <- rbind(wwarn_k13_clean, wwarn_pd_clean)
saveRDS(wwarn_clean, "analysis/data-derived/wwarn_clean.rds")
