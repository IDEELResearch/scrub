# Load required libraries
library(dplyr)
library(lubridate)
library(here)

################################################################################
#
# Script: 01_read_geoffs.R
# Purpose: Combine all geoff studies into one dataframe
#
################################################################################

# Define root directory and initialize empty dataframe for combined all studies wide table
root_dir <- here("analysis", "data-geoff")
master_table <- data.frame()

# Iterate over each project directory within root_dir
study_folders <- list.dirs(root_dir, recursive = FALSE)
for (project_dir in study_folders) {
  
  # Extract the project identifier
  study_uid <- basename(project_dir)
  print(paste("Processing project:", study_uid))
  
  # Identify file paths
  study_overview_path <- list.files(
    path = project_dir,
    pattern = "*_study_data_validated\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  site_overview_path <- list.files(
    path = project_dir,
    pattern = "*_site_data_validated\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  prev_table_path <- list.files(
    path = project_dir,
    pattern = "*_prevalence_data_LONG_validated\\.tsv\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Skip if any files are missing
  if (length(study_overview_path) == 0 || length(site_overview_path) == 0 || length(prev_table_path) == 0) {
    print(paste("Project missing all three input tables:", study_uid))
    next
  }
  
  # Read in the TSV files
  study_overview <- read.table(gzfile(study_overview_path), sep = "\t", quote = "", header = TRUE, fill = TRUE)
  site_overview <- read.table(gzfile(site_overview_path), sep = "\t", quote = "", header = TRUE, fill = TRUE)
  prev_table <- read.table(gzfile(prev_table_path), sep = "\t", quote = "", header = TRUE, fill = TRUE)
  
  # Transpose and merge data
  
  study_overview_t <- as.data.frame(t(study_overview$DATA), stringsAsFactors = FALSE)
  colnames(study_overview_t) <- (study_overview$FIELDS)
  study_overview_t <- janitor::clean_names(study_overview_t)
  
  # TODO: fix the duplicate columns in the outputs created
  wide_data <- prev_table %>%
    left_join(site_overview, by = "site_uid") %>%
    mutate(study_key = study_uid) %>%
    cbind(study_overview_t) %>%
    as.data.frame(check.names = TRUE)
  
  # Convert numeric columns to character
  wide_data <- wide_data %>%
    mutate(across(where(is.numeric), as.character))
  
  # Append data to master_table
  master_table <- bind_rows(master_table, wide_data)
}

# Save the combined master table as an RDS file in the data-derived directory
saveRDS(master_table, here("analysis", "data-derived", "01_read_geoffs_output_table.rds"))
print("Master table saved. Ready for cleaning with 02_clean_geoffs.R.")
