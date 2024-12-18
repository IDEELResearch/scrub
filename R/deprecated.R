# this is a collection of deprecated code from other files not currently needed
# this has been created to make it easier to get the package to pass checks. 


# read_geoff <- function(study_path, site_path, prev_path) {
#   study <- readr::read_tsv(file = study_path, show_col_types = FALSE)
#   site <- readr::read_tsv(file = site_path, show_col_types = FALSE)
#   prev <- readr::read_tsv(file = prev_path, show_col_types = FALSE)
#   # pull out study ID for the rest of the data frame
#   study_id <- study$DATA[which(study$FIELDS == "study_uid")]
#   author <- study$DATA[which(study$FIELDS == "first_author_surname")]
#   
#   # studies dataframe
#   study_names <- study$FIELDS # safe column names
#   study <- t(study[,-1]) %>% as.data.frame() # convert to df
#   colnames(study) <- study_names
#   rownames(study) <- NULL
#   studies <- study %>%
#     # select the columns we want to keep [PRELIMINARY] -- this is determined by {stave} structure
#     dplyr::select(c(study_uid, first_author_surname, publication_year, pmid, publication_status))  %>%
#     dplyr::mutate(database = "geoff")
#   
#   # surveys dataframe is a combination of sites and prevalence
#   prev <- prev %>%
#     dplyr::group_by(across(c(site_uid, date_start, date_end))) %>%
#     dplyr::mutate(substudy = cur_group_id()) %>% # TODO: change this
#     dplyr::mutate(survey_id = paste0(study_id,"-",author,"-",site_uid,"-",date_end)) # check what happens if date_end is a date and not a character string
#   
#   # TODO: make a unique start end numbering and use that for the survey id 
#   site_prev <- full_join(site, prev, by = "site_uid")
#   surveys <- site_prev %>%
#     dplyr::mutate(study_id = study_id) %>%
#     dplyr::select(c(study_id, survey_id, iso3c, country, site_name, collection_location,
#                     lat_n, lon_e, date_start, date_end)) %>%
#     dplyr::distinct() %>%
#     dplyr::rename(site = site_name,
#                   location = collection_location,
#                   lat = lat_n,
#                   long = lon_e)
#   
#   counts <- site_prev %>%
#     dplyr::mutate(study_uid = study_id) %>%
#     dplyr::select(study_uid, survey_id, gene_mutation, mutant_num, total_num)
#   
#   df <- full_join(study, surveys) %>%
#     dplyr::full_join(counts)
#   
#   
# }



# placeholder function to get a clean geoff 

# TODO: make this compatible with geoff data. at the moment this is just a placeholder copied from wwarn
# geoff_to_stave <- function(geoff) {
#   #TODO figure out how to extract the publication year from PMID
#   
#   studies <- geoff %>% 
#     dplyr::select(c("study_ID","study_name","study_type","authors","publication_year","url")) %>%
#     dplyr::distinct(study_ID, study_name, study_type, authors, publication_year, .keep_all = TRUE)
#   
#   counts <- geoff %>%
#     dplyr::rename(survey_key = survey_id,
#                   variant_string = gene_mut,
#                   variant_num = x,
#                   total_num = n) %>%
#     dplyr::filter(variant_num >= 0) %>%
#     # TODO: figure out why some of these are not integers and how to fix them
#     dplyr::mutate(variant_num = floor(variant_num)) %>%
#     dplyr::filter(variant_string != "mdr1:CNV") %>%
#     dplyr::select(survey_key, variant_string, variant_num, total_num) %>%
#     dplyr::distinct(survey_key, variant_string, .keep_all = TRUE)
#   
#   
#   ### surveys id -- in geoff a lot of data points are multisite 
#   surveys <- geoff %>%
#     dplyr::group_by(across(c(study_ID, site_fixed))) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(-site) %>%
#     dplyr::rename(site_name = admin_1, 
#                   date_start = study_start_year,
#                   date_end = study_end_year,
#                   lon = long, 
#                   survey_ID = survey_id) %>%
#     dplyr::filter(survey_ID %in% counts$survey_key) %>%
#     dplyr::mutate(collection_start = paste0(date_start,"-01-01"),
#                   collection_end = paste0(date_end,"-12-31"),
#                   study_key = study_ID) %>%
#     dplyr::mutate(start = as.Date(collection_start),
#                   end = as.Date(collection_end),
#                   mid = as.Date((as.numeric(start) + as.numeric(end)) / 2)) %>%
#     dplyr::mutate(collection_day = as.character(mid),
#                   spatial_notes = "geoff lat and long",
#                   time_notes = "automated midpoint") %>%
#     dplyr::select(study_key, survey_ID, country_name, site_name,
#                   lat, lon, spatial_notes, collection_start, collection_end, 
#                   collection_day, time_notes) %>%
#     dplyr::distinct(study_key, survey_ID, .keep_all = TRUE)
#   
#   studies <- studies %>%
#     dplyr::filter(study_ID %in% surveys$study_key)  
#   
#   return(list(studies_dataframe = studies, 
#               surveys_dataframe = surveys, 
#               counts_dataframe = counts))
# }
