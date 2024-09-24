## example workflow to read in all the data files and convert them to a stave object
devtools::load_all() # make sure all functions & packages are loaded that we need in this library

# load stave - currently STAVE is empty 
devtools::install_github("mrc-ide/STAVE")
library(STAVE)

# TODO: this will all need fixing once Bob fixes stave
# TODO: all of the data should be combined into one and then cleaned and deduplicated before conversion to stave
# read in data from databases (excluding geoff)
wwarn_data <- readRDS("analysis/data-derived/wwarn_res_df.rds")

# make an empty stave object
stave <- STAVE::STAVE_object$new()
# append WWARN data
wwarn_stave <- wwarn_to_stave(wwarn_data = wwarn_data)
stave <- stave$append_data(studies_dataframe = wwarn_stave$studies_dataframe,
                           surveys_dataframe = wwarn_stave$surveys_dataframe,
                           counts_dataframe = wwarn_stave$counts_dataframe)



# now convert all of the geoff objects into a list
geoff_path <- "analysis/data-geoff/" # filepath to all of the geoff tsvs which are in separate folders
folders <- list.dirs(path = geoff_path)

geoff <- list(studies_dataframe = data.frame(),
              surveys_dataframe = data.frame(),
              counts_dataframe = data.frame())

for(i in 2:length(folders)) { # the first path is the root which we want to omit
  files <- list.files(path = folders[i],
                      all.files = TRUE, full.names = TRUE) %>%
    str_subset(".tsv")
  # relies on a standardised naming convention of the .tsv files
  study_path <- str_subset(files, "study")
  site_path <- str_subset(files, "site")
  prev_path <- str_subset(files, "prevalence")
  
  geoff_obj <- read_geoff(study_path = study_path,
                          site_path = site_path,
                          prev_path = prev_path)
  
  # this is now a list of data frames ready for STAVE
  geoff <- list(studies_dataframe = rbind(geoff$studies_dataframe, geoff_obj$studies_dataframe),
                surveys_dataframe = rbind(geoff$surveys_dataframe, geoff_obj$surveys_dataframe),
                counts_dataframe = rbind(geoff$counts_dataframe, geoff_obj$counts_dataframe))
  
}

stave <- STAVE::append_data(stave, geoff)

