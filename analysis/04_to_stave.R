# load stave  
# devtools::install_github("mrc-ide/STAVE")
library(STAVE)

# TODO: all of the data should be combined into one and then cleaned and deduplicated before conversion to stave
# read in data from databases (excluding geoff)
data <- readRDS("analysis/data-derived/final_data.rds")

# make an empty stave object
stave <- STAVE::STAVE_object$new()
# append WWARN data
data_stave <- convert_stave(data)

# TODO: this throws errors in specific rows -- have to debug in console to work
stave$append_data(studies_dataframe = data_stave$studies_dataframe,
                  surveys_dataframe = data_stave$surveys_dataframe,
                  counts_dataframe = data_stave$counts_dataframe)

dir.create("analysis/data-out", showWarnings = FALSE)
saveRDS(stave, "analysis/data-out/stave_data.rds")
