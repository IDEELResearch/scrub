#' Convert data into a stave compatible list of dataframes
#' 
#' @param data A wide dataset to be converted into a form compatible with {STAVE}
#' 
#' @return A list of three dataframes easily convertible to be compatible with {STAVE} using stave$append_data
#' 
#' @export
#' 

convert_stave <- function(data) {
  
  #TODO figure out how to extract the publication year from PMID (is this still needed)
  
  # grab our distinct studies data frame 
  studies <- data %>% 
    dplyr::select(c("study_ID","study_name","study_type","authors","publication_year","url")) %>%
    dplyr::mutate(publication_year = as.integer(publication_year)) %>% 
    dplyr::distinct(study_ID, study_name, study_type, authors, publication_year, .keep_all = TRUE)
  
  # grab our distinct counts data frame
  counts <- data %>%
    dplyr::rename(survey_key = survey_ID) %>% 
    dplyr::filter(variant_string != "mdr1:CNV") %>%
    dplyr::filter(variant_string != "pfpm23:CNV") %>%
    dplyr::select(survey_key, variant_string, variant_num, total_num) %>%
    dplyr::mutate(variant_num = as.integer(variant_num)) %>% 
    dplyr::mutate(total_num = as.integer(total_num)) %>% 
    dplyr::distinct(survey_key, variant_string, .keep_all = TRUE)
  
  # grab our distinct surveys data frame
  surveys <- data %>%
    dplyr::rename(study_key = study_ID) %>%
    dplyr::select(study_key, survey_ID, country_name, site_name,
                  lat, lon, spatial_notes, collection_start, collection_end, 
                  collection_day, time_notes) %>%
    dplyr::mutate(lat = as.integer(lat)) %>% 
    dplyr::mutate(lon = as.integer(lon)) %>% 
    dplyr::distinct(study_key, survey_ID, .keep_all = TRUE)
  
  studies <- studies %>%
    dplyr::filter(study_ID %in% surveys$study_key)  

  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}
