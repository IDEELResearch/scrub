wwarn_data <- wwarn_res_df # readRDS("analysis/data-derived/wwarn_res_df.rds")

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
  dplyr::mutate(long = long%%360) %>%
  dplyr::rename(study_ID = study_uid) %>%
  dplyr::mutate(study_ID = paste0(study_ID, "_", authors),
                study_name = study_ID,
                study_type = "peer_reviewed") %>%
  dplyr::mutate(study_ID = iconv(study_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
  dplyr::mutate(study_ID =  gsub("[^a-zA-Z0-9_]", "", study_ID)) 

saveRDS(wwarn, "analysis/data-derived/wwarn_stave.RDS")
