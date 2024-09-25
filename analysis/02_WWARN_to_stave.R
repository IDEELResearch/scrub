wwarn_data <- readRDS("analysis/data-derived/wwarn_res_df.rds")
ref_als <- read.csv("analysis/data-raw/k13_ref_protein_codon_dictionary.csv")

wwarn <- wwarn_data %>%
  dplyr::group_by(pmid) %>%
  dplyr::mutate(study_uid = paste0("wwarn","_",sample(10000:99999,1, replace = FALSE))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(publication_status = "peer-reviewed", 
                database = "wwarn",
                publication_year = 1000) %>% # placeholder but this is so stave doesnt error
  dplyr::rowwise() %>%
  dplyr::mutate(source = str_split(source, " ")[[1]][1]) %>%
  dplyr::mutate(site_fixed = if_else(is.na(site), 
                                     gsub(" ","",admin_1), 
                                     gsub(" ", "", site)),
                country_name = countrycode::countrycode(iso3c, # fix up top
                                                        origin = 'iso3c', 
                                                        destination = 'country.name')) %>%
  # manually fix the country names
  dplyr::mutate(country_name = if_else(country_name %in% c("Congo - Kinshasa",
                                                           "Congo - Brazzaville"),
                                       "Democratic Republic of the Congo", country_name)) %>%
  dplyr::mutate(country_name = if_else(country_name == "Côte d’Ivoire",
                                       "Côte d'Ivoire", country_name)) %>%
  dplyr::mutate(country_name = if_else(country_name == "Myanmar (Burma)",
                                       "Myanmar", country_name)) %>%
  dplyr::mutate(country_name = if_else(country_name == "São Tomé & Príncipe",
                                       "São Tomé and Príncipe", country_name)) %>%
  dplyr::mutate(country_name = if_else(country_name == "Cape Verde",
                                       "Cabo Verde", country_name)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(country_name != "Eswatini") %>% 
  dplyr::mutate(study_end_year = if_else(is.na(study_end_year),
                                         as.numeric(year), study_end_year)) %>%
  dplyr::mutate(study_start_year = if_else(is.na(study_start_year),
                                         as.numeric(year), study_start_year)) %>%
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
  dplyr::mutate(study_ID =  gsub("[^a-zA-Z0-9_]", "", study_ID))

# fix the wildtype mutations
wt_studies <- wwarn %>%
  dplyr::filter(gene_mut == "k13:WT") %>%
  dplyr::pull(study_ID) %>% unique()

# for each of these studies find out what the wildtype encoding should be
wt_mutations <- data.frame(study_ID = wt_studies,
                           mutations = rep("", length = length(wt_studies)))

for(i in 1:length(wt_studies)) {
  muts <- wwarn %>%
    dplyr::filter(study_ID == wt_studies[i]) %>%
    dplyr::select(study_ID, gene_mut) %>%
    dplyr::filter(str_detect(gene_mut, "k13")) %>%
    dplyr::mutate(codons = gsub("k13:","", gene_mut)) %>%
    dplyr::mutate(codons = gsub("[A-Za-z]","", codons)) %>%
    dplyr::mutate(codons = gsub(":","", codons)) %>%
    dplyr::filter(codons != "") %>%
    dplyr::pull(codons) %>% parse_number() %>% unique() %>% sort()
  # TODO: figure out how to add the reference allelles here 
  wt_mutations$mutations[i] <- paste0("k13:",paste(muts, collapse = "_"),":*")
  
}

wwarn_mut <- wwarn %>% dplyr::filter((gene_mut == "k13:WT") == FALSE) %>%
  dplyr::filter(gene_mut != "k13:470:X") %>%
  dplyr::distinct(survey_id, gene_mut, .keep_all = TRUE)
wwarn_wt <- wwarn %>% dplyr::filter((gene_mut == "k13:WT")) %>%
  dplyr::left_join(wt_mutations, by = "study_ID") %>%
  dplyr::filter(mutations != "k13::") %>%
  dplyr::mutate(gene_mut = mutations) %>%
  dplyr::select(-mutations)

# TODO: figure out why these are errors but just filter out for now so I can send an .rds
remove <- wwarn_stave$counts_dataframe[c(3714, 3717, 3727, 3750, 3752, 3756, 3762, 3767, 
                                         3776, 3785, 3789, 3801, 3813, 3815, 3486, 2490, 
                                         2483, 2930, 3088, 3107, 3111, 3202, 2672, 3642, 
                                         3644, 2772, 2982, 2715, 2698, 3417, 3405, 2471),] %>% 
  pull(survey_key)
  
  
  
# wwarn <- rbind(wwarn_mut, wwarn_wt) 
wwarn <- wwarn_mut %>% dplyr::filter((survey_id %in% remove) == FALSE)

saveRDS(wwarn, "analysis/data-derived/wwarn_stave.RDS")
