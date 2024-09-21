## functions converting the main types of data output into stave objects

wwarn_to_stave <- function(wwarn_data) {
  #TODO figure out how to extract the publication year from PMID
  wwarn <- wwarn_data %>%
    dplyr::group_by(pmid) %>% # study uid is unique per article
    dplyr::mutate(study_uid = paste0("wwarn","-",sample(10000:99999,1, replace = FALSE))) %>%
    dplyr::mutate(publication_status = "peer-reviewed", 
                  database = "wwarn",
                  publication_year = "*PLACEHOLDER") %>% # placeholder until I figure this out
    dplyr::rowwise() %>%
    dplyr::mutate(source = str_split(source, " ")[[1]][1]) %>%
    dplyr::mutate(site_fixed = if_else(is.na(site), 
                                       gsub(" ","",admin_1), 
                                       gsub(" ", "", site))) %>%
    dplyr::ungroup() %>%
    dplyr::rename(first_author_surname = source) %>%
    # now add survey IDs here also
    dplyr::group_by(across(c(study_uid, site, study_end_year))) %>%
    dplyr::mutate(survey_id = paste0(study_uid,"-",author,"-",site_fixed,"-",study_end_year)) %>%
    dplyr::ungroup()
  
  studies <- wwarn %>% 
    select(c(study_uid, first_author_surname, publication_year, 
             pmid, publication_status, database))
  
  ### surveys id -- in WWARN a lot of data points are multisite 
  # TODO - figure out what to do about surveys that are multiple sites
  surveys <- wwarn %>%
    dplyr::group_by(across(c(study_uid, site_fixed))) %>%
    dplyr::mutate(survey_id = paste0(study_uid,"-",first_author_surname,"-",
                                     site_fixed,"-",study_end_year)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-site) %>%
    dplyr::rename(site = admin_1, 
                  country = admin_0,
                  location = site_fixed,
                  date_start = study_start_year,
                  date_end = study_end_year) %>%
    dplyr::select(study_uid, survey_id, iso3c, country, site, location, 
                  lat, long, date_start,date_end)
  
  counts <- wwarn %>%
    dplyr::rename(gene_mutation = gene_mut,
                  mutant_num = x,
                  total_num = n) %>%
    dplyr::select(study_uid, survey_id, gene_mutation, mutant_num, total_num)
  
  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}

# repeat the same functions for pf7k and who when those datasets are ready
