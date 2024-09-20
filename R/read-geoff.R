## reading in the tsv's routinely fixing the formatting
# read in all three tsv's from the filepaths passed and then reformat into a list of dataframes for stave

read_geoff <- function(study_path, site_path, prev_path) {
  study <- readr::read_tsv(file = study_path, show_col_types = FALSE)
  site <- readr::read_tsv(file = site_path, show_col_types = FALSE)
  prev <- readr::read_tsv(file = prev_path, show_col_types = FALSE)
  # pull out study ID for the rest of the data frame
  study_id <- study$DATA[which(study$FIELDS == "study_uid")]
  author <- study$DATA[which(study$FIELDS == "first_author_surname")]
  
  # studies dataframe
  study_names <- study$FIELDS # safe column names
  study <- t(study[,-1]) %>% as.data.frame() # convert to df
  colnames(study) <- study_names
  rownames(study) <- NULL
  studies <- study %>%
    # select the columns we want to keep [PRELIMINARY] -- this is determined by {stave} structure
    select(c(study_uid, first_author_surname, publication_year, pmid, publication_status)) 
  
  # surveys dataframe is a combination of sites and prevalence
  prev <- prev %>%
    dplyr::group_by(across(c(site_uid, date_start, date_end))) %>%
    dplyr::mutate(substudy = cur_group_id()) %>% # fix substudy numbering
    dplyr::mutate(survey_id = paste0(study_id,"-",author,"-",site_uid,"-",date_end)) # check what happens if date_end is a date and not a character string
  
  site_prev <- full_join(site, prev, by = "site_uid")
  surveys <- site_prev %>%
    dplyr::mutate(study_id = study_id) %>%
    dplyr::select(c(study_id, survey_id, iso3c, country, site_name, collection_location,
                      lat_n, lon_e, date_start, date_end)) %>%
    dplyr::distinct() %>%
    dplyr::rename(site = site_name,
                  location = collection_location,
                  lat = lat_n,
                  long = lon_e)
  
  counts <- site_prev %>%
    dplyr::mutate(study_id = study_id) %>%
    dplyr::select(survey_id, gene_mutation, mutant_num, total_num)
  
  return(list(studies_dataframe = studies, 
              surveys_dataframe = surveys, 
              counts_dataframe = counts))
}

# test this object
stave_obj <- read_geoff(study_path = "analysis/data-example/ideel_study_data_validated.tsv",
                        site_path =  "analysis/data-example/ideel_site_data_validated.tsv",
                        prev_path = "analysis/data-example/ideel_prevalence_data_long_validated.tsv")
