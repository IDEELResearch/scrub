# Load required packages --------------------------------------------
library(dplyr)
library(lubridate)

#' Clean general data by trimming whitespace, replacing spaces, and converting empty strings to NA
#'
#' This function takes a data frame and applies cleaning operations to all character columns.
#' It trims leading/trailing whitespace, replaces spaces with underscores, and converts empty strings to NA.
#'
#' @param df A data frame to be cleaned.
#' @return A cleaned data frame with modified character columns.
#' @export
clean_data <- function(df) {
  df[] <- lapply(df, function(x) if (is.character(x)) trimws(gsub(" ", "_", x)) else x)
  df[] <- lapply(df, function(x) if (is.character(x)) ifelse(x == "", NA, x) else x)
  return(df)
}

#' Adjust invalid dates by setting date_start to the first day and date_end to the last day of the month
#'
#' This function adjusts invalid dates by converting them to the first day of the month (for start dates) 
#' and the last day of the month (for end dates).
#'
#' @param date_str A character vector of date strings.
#' @param mode A character string specifying whether it's "start" or "end" date adjustment. 
#'             "start" sets the date to the first day of the month, "end" to the last day.
#' @return A Date object with invalid dates adjusted accordingly.
#' @export
adjust_invalid_date <- function(date_str, mode = c("start", "end")) {
  mode <- match.arg(mode)
  date_fixed <- suppressWarnings(ymd(date_str))
  is_invalid <- is.na(date_fixed)
  
  if (any(is_invalid)) {
    year_month <- substr(date_str[is_invalid], 1, 7)
    
    if (mode == "start") {
      # Set the date to the first day of the month
      date_fixed[is_invalid] <- ymd(paste0(year_month, "-01"))
    } else if (mode == "end") {
      # Set the date to the last day of the month
      last_day_of_month <- sapply(year_month, function(ym) days_in_month(as.Date(paste0(ym, "-01"))))
      date_fixed[is_invalid] <- ymd(paste0(year_month, "-", last_day_of_month))
    }
  }
  
  return(date_fixed)
}

# Set project directory and file paths ------------------------------
project_dir <- "/Users/george/Bailey_lab/GEM_hackathon/scrubbing_geoff_dev/s0022_molinadelafuente_2023_v01"
study_uid <- basename(project_dir)

study_overview_path <- paste0(study_uid, "_study_data_validated.tsv")
site_overview_path <- paste0(study_uid, "_site_data_validated.tsv")
prev_table_path <- paste0(study_uid, "_prevalence_data_LONG_validated.tsv")

#' Read and transpose study overview data ---------------------------
#' 
#' Read the study overview TSV and transpose the data for further processing.
#' The first row is removed and the column names are set to the first row of the original data.
#'
#' @param project_dir Directory where the TSV file is located.
#' @param study_overview_path File path to the study overview data.
#' @return A transposed data frame of the study overview.
#' @export
study_overview <- read.table(file.path(project_dir, study_overview_path), sep = "\t", header = TRUE)
site_overview <- read.table(file.path(project_dir, site_overview_path), sep = "\t", header = TRUE)
prev_table <- read.table(file.path(project_dir, prev_table_path), sep = "\t", header = TRUE)

# Transpose study_overview
study_overview_t <- as.data.frame(t(study_overview))
colnames(study_overview_t) <- study_overview$FIELDS
study_overview_t <- study_overview_t[-1, , drop = FALSE]  # Remove the first row

# Clean and merge the data ------------------------------------------
study_overview_t <- clean_data(study_overview_t)
site_overview <- clean_data(site_overview)
prev_table <- clean_data(prev_table)

wide_data <- prev_table %>%
  left_join(site_overview, by = "site_uid") %>%
  mutate(study_key = study_overview_t$study_uid) %>%
  {
    # Use withCallingHandlers to catch and suppress the specific warning
    withCallingHandlers(
      cbind(., study_overview_t),  # The operation potentially generating the warning
      warning = function(w) {
        if (grepl("row names were found from a short variable and have been discarded", w$message)) {
          invokeRestart("muffleWarning")  # Suppress only this specific warning
        }
      }
    )
  }

#' Clean, convert dates, and calculate midpoints ---------------------
#' 
#' This step cleans start and end date columns, adjusts invalid dates, and calculates
#' the midpoint of collection periods. It also converts date columns to characters.
#'
#' @return A data frame with cleaned and adjusted dates and midpoints.
wide_data <- wide_data %>%
  mutate(
    # Adjust start date to the first day of the month
    collection_start = adjust_invalid_date(date_start, mode = "start"),
    
    # Adjust end date to the last day of the month
    collection_end = adjust_invalid_date(date_end, mode = "end"),
    
    collection_start = ymd(collection_start),
    collection_end = ymd(collection_end),
    
    # Calculate midpoint
    mid = case_when(
      !is.na(collection_start) & !is.na(collection_end) ~ as.Date((as.numeric(collection_start) + as.numeric(collection_end)) / 2, origin = "1970-01-01"),
      !is.na(collection_start) & is.na(collection_end) ~ collection_start,
      is.na(collection_start) & !is.na(collection_end) ~ collection_end,
      TRUE ~ NA_Date_
    )
  ) %>%
  # Convert dates to character strings for output
  mutate(
    collection_start = as.character(collection_start),
    collection_end = as.character(collection_end),
    mid = as.character(mid)
  )


# Generate unique survey IDs -----------------------------------------
set.seed(1234)
random_ids <- sample(1:1000000, size = nrow(wide_data), replace = FALSE)

wide_data <- wide_data %>%
  mutate(
    survey_ID = paste0(study_uid, "_", first_author_surname, "_", site_name, "_", publication_year, "_", random_ids),
    survey_ID = gsub("[^a-zA-Z0-9_]", "", survey_ID),
    survey_ID = iconv(survey_ID, from = "UTF-8", to = "ASCII//TRANSLIT")
  )

# Create data frames for output --------------------------------------
studies <- data.frame(
  study_ID = study_overview_t$study_uid,
  study_name = paste(study_overview_t$first_author_surname, "_et al_", study_overview_t$publication_year, sep = ""),
  study_type = "other",
  authors = study_overview_t$first_author_surname,
  publication_year = as.numeric(study_overview_t$publication_year),
  url = study_overview_t$study_url
)

surveys <- wide_data %>%
  dplyr::mutate(
    spatial_notes = "lat_and_long",
    time_notes = "automated_midpoint"
  ) %>%
  dplyr::select(
    study_key,
    survey_ID,
    country_name = country,
    site_name,
    lat = lat_n,
    lon = lon_e,
    spatial_notes,
    collection_start,
    collection_end,
    collection_day = mid,
    time_notes
  )

#' Clean amino acid sequence -----------------------------------------
#' 
#' Ensures that the amino acid sequence is valid and in uppercase.
#' 
#' @param sequence A character vector of amino acid sequences.
#' @return A cleaned character vector of amino acid sequences.
#' @export
clean_amino_acid_sequence <- function(sequence) {
  valid_chars <- allowed_amino_acids()$amino_acid
  cleaned_sequence <- gsub("[^ACDEFGHIKLMNPQRSTVWY|/]", "", toupper(sequence))
  cleaned_sequence <- gsub("_", "", cleaned_sequence)
  return(cleaned_sequence)
}

# Process counts dataframe -------------------------------------------
counts <- wide_data %>%
  dplyr::mutate(
    survey_key = survey_ID,
    gene_mutation = tolower(gene_mutation),
    variant_string = sapply(gene_mutation, function(gene) {
      parts <- unlist(strsplit(gene, ":"))
      if (length(parts) == 3) {
        gene_name <- parts[1]
        amino_acid_position <- parts[2]
        amino_acid_sequence <- clean_amino_acid_sequence(parts[3])
        return(paste(gene_name, amino_acid_position, amino_acid_sequence, sep = ":"))
      } else {
        return(gene)
      }
    })
  ) %>%
  dplyr::select(
    survey_key,
    variant_string,
    variant_num = mutant_num,
    total_num
  )

# Output the processed dataframes ------------------------------------
str(studies)
str(surveys)
str(counts)

# Append to STAVE object ---------------------------------------------
stave <- STAVE::STAVE_object$new()
geoff_stave <- list(studies_dataframe = studies, surveys_dataframe = surveys, counts_dataframe = counts)

stave$append_data(studies_dataframe = geoff_stave$studies_dataframe,
                  surveys_dataframe = geoff_stave$surveys_dataframe,
                  counts_dataframe = geoff_stave$counts_dataframe)
