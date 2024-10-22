## functions converting the main types of data output into stave objects

wwarn_to_stave <- function(wwarn) {
  #TODO figure out how to extract the publication year from PMID
  
  studies <- wwarn %>% 
    dplyr::select(c("study_ID","study_name","study_type","authors","publication_year","url")) %>%
    dplyr::distinct(study_ID, study_name, study_type, authors, publication_year, .keep_all = TRUE)
  
  counts <- wwarn %>%
    dplyr::rename(survey_key = survey_id,
                  variant_string = gene_mut,
                  variant_num = x,
                  total_num = n) %>%
    dplyr::filter(variant_num >= 0) %>%
    # TODO: figure out why some of these are not integers and how to fix them
    dplyr::mutate(variant_num = floor(variant_num)) %>%
    dplyr::filter(variant_string != "mdr1:CNV") %>%
    dplyr::select(survey_key, variant_string, variant_num, total_num) %>%
    dplyr::distinct(survey_key, variant_string, .keep_all = TRUE)
  
  
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
                  collection_day, time_notes) %>%
    dplyr::distinct(study_key, survey_ID, .keep_all = TRUE)
  
  studies <- studies %>%
    dplyr::filter(study_ID %in% surveys$study_key)  

  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}

# repeat the same functions for pf7k and who when those datasets are ready
