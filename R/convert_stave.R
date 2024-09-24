## functions converting the main types of data output into stave objects

wwarn_to_stave <- function(wwarn_data) {
  #TODO figure out how to extract the publication year from PMID
  wwarn <- wwarn_data %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
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
    dplyr::filter(study_end_year >= 1980) %>%
    dplyr::rename(first_author_surname = source) %>%
    # now add survey IDs here also
    dplyr::group_by(across(c(study_uid, site, study_end_year))) %>%
    dplyr::mutate(survey_id = paste0(study_uid,"_",author,"_",site_fixed,"_",study_end_year)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(survey_id = gsub("[^a-zA-Z0-9_]", "", survey_id)) %>% # all the things throwing errors in stave
    dplyr::rowwise() %>%
    dplyr::mutate(survey_id = iconv(survey_id, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(long = as.numeric(long),
                  lat = as.numeric(lat)) %>%
    dplyr::mutate(long = long%%360) %>%
    dplyr::rename(study_ID = study_uid,
                  authors = first_author_surname) %>%
    dplyr::mutate(study_ID = paste0(study_ID, "_", authors),
                  study_name = study_ID,
                  study_type = "peer_reviewed") %>%
    dplyr::mutate(study_ID = iconv(study_ID, from = "UTF-8", to = "ASCII//TRANSLIT")) %>%
    dplyr::mutate(study_ID =  gsub("[^a-zA-Z0-9_]", "", study_ID))
    
  studies <- wwarn %>% 
    dplyr::select(c("study_ID","study_name","study_type","authors","publication_year","url")) %>%
    dplyr::distinct()

  counts <- wwarn %>%
    dplyr::rename(survey_key = survey_id,
                  variant_string = gene_mut,
                  variant_num = x,
                  total_num = n) %>%
    dplyr::filter(variant_num >= 0) %>%
    # TODO: figure out why some of these are not integers and how to fix them
    dplyr::mutate(variant_num = floor(variant_num)) %>%
    dplyr::filter(variant_string != "mdr1:CNV") %>%
    dplyr::select(survey_key, variant_string, variant_num, total_num) 
  
  
  ### surveys id -- in WWARN a lot of data points are multisite 
  surveys <- wwarn %>%
    dplyr::group_by(across(c(study_ID, site_fixed))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-site) %>%
    dplyr::rename(site_name = admin_1, 
                  date_start = study_start_year,
                  date_end = study_end_year,
                  lon = long, 
                  survey_ID = survey_id) %>%
    dplyr::filter(survey_ID %in% counts$survey_key) %>%
    dplyr::mutate(collection_start = paste0(date_start,"-01-01"),
                  collection_end = paste0(date_end,"-12-31"),
                  study_key = study_ID) %>%
    dplyr::mutate(start = as.Date(collection_start),
                  end = as.Date(collection_end),
                  mid = as.Date((as.numeric(start) + as.numeric(end)) / 2)) %>%
    dplyr::mutate(collection_day = as.character(mid),
                  spatial_notes = "wwarn lat and long",
                  time_notes = "automated midpoint") %>%
    dplyr::select(study_key, survey_ID, country_name, site_name,
                  lat, lon, spatial_notes, collection_start, collection_end, 
                  collection_day, time_notes)
  studies <- studies %>%
    dplyr::filter(study_ID %in% surveys$study_key)

  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}

# repeat the same functions for pf7k and who when those datasets are ready
