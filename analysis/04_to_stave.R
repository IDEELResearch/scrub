# load stave  
# devtools::install_github("mrc-ide/STAVE")
library(STAVE)

# read in data from databases (excluding geoff)
data <- readRDS("analysis/data-derived/final_data.rds")

# hot fixes till upstream data entry and STAVE is fixed
data <- data %>% filter(!grepl("-|TK|BND", variant_string))
data <- data %>% filter(!(survey_ID %in% (distinct_variants %>% filter(distinct_total_num>1) %>% pull(survey_key))))

# make an empty stave object
stave <- STAVE::STAVE_object$new()

# append WWARN data
data_stave <- convert_stave(data)

# Convert into a STAVE object
stave$append_data(studies_dataframe = data_stave$studies_dataframe,
                  surveys_dataframe = data_stave$surveys_dataframe,
                  counts_dataframe = data_stave$counts_dataframe[-c(11975, 14542, 14558, 14561, 14550, 14566, 14574, 14582),]) # another hot fix

# Save the output in data-out ready for downstream analysis
dir.create("analysis/data-out", showWarnings = FALSE)
saveRDS(stave, "analysis/data-out/stave_data.rds")
