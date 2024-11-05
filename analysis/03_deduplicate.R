library(here)

# read in the data
# TODO: read in the other data sources
wwarn_data <- readRDS(here("analysis", "data-derived", "wwarn_data.rds"))

## placeholder for reading in all the data sources and deduplicating
# TODO: deduplicate the studies
all_data <- wwarn_data

# output final data
saveRDS(all_data, "analysis/data-derived/final_data.rds")
