#' Convert data into a stave compatible list of dataframes
#' 
#' @param data A wide dataset to be converted into a form compatible with `STAVE`
#' 
#' @return A list of three dataframes easily convertible to be compatible with `STAVE` using stave$append_data
#' 
#' @export
#' 
convert_stave <- function(data) {
  
  #TODO figure out how to extract the publication year from PMID (is this still needed?)
  
  # grab our distinct studies data frame 
  studies <- data %>% 
    dplyr::select(c("study_ID","study_name","study_type","authors","publication_year","url")) %>%
    dplyr::mutate(publication_year = as.integer(.data$publication_year)) %>% 
    dplyr::distinct(.data$study_ID, .data$study_name, .data$study_type, .data$authors, .data$publication_year, .keep_all = TRUE) %>% 
    dplyr::rename(study_id = study_ID) %>% 
    dplyr::mutate(study_type = replace(.data$study_type, .data$study_type == "published", "peer_reviewed"))
  
  # grab our distinct counts data frame
  counts <- data %>%
    dplyr::rename(survey_key = .data$survey_ID) %>% 
    dplyr::filter(.data$variant_string != "mdr1:CNV") %>%
    dplyr::filter(.data$variant_string != "pfpm23:CNV") %>%
    dplyr::select(study_key = .data$study_ID, .data$survey_key, .data$variant_string, .data$variant_num, .data$total_num) %>%
    dplyr::mutate(variant_num = as.integer(.data$variant_num)) %>% 
    dplyr::mutate(total_num = as.integer(.data$total_num)) %>% 
    dplyr::distinct(.data$survey_key, .data$variant_string, .keep_all = TRUE)
  
  # grab our distinct surveys data frame
  surveys <- data %>%
    dplyr::rename(study_key = .data$study_ID) %>%
    dplyr::select(.data$study_key, .data$survey_ID, .data$country_name, .data$site_name,
                  latitude = .data$lat, longitude = .data$lon, .data$spatial_notes, 
                  .data$collection_start, .data$collection_end, 
                  .data$collection_day, .data$time_notes) %>%
    dplyr::mutate(latitude = as.numeric(.data$latitude)) %>% 
    dplyr::mutate(longitude = as.numeric(.data$longitude)) %>% 
    dplyr::distinct(.data$study_key, .data$survey_ID, .keep_all = TRUE) %>% 
    dplyr::rename(survey_id = survey_ID)

  
  studies <- studies %>%
    dplyr::filter(.data$study_id %in% surveys$study_key)  

  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}
