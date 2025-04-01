#' Deduplicate Data in a Principled Way
#'
#' This function deduplicates a dataset, ensuring each entry is unique based on 
#' specified criteria, and returns a cleaned data frame suitable for input into STAVE.
#'
#' @param df A data.frame to be deduplicated. It must contain the following columns:
#'   - `"iso3c"`: ISO3 country code
#'   - `"admin_0"`: Country-level administrative division
#'   - `"admin_1"`: First-level administrative division
#'   - `"site"`: Site name or identifier
#'   - `"lat"`: Latitude of the site
#'   - `"long"`: Longitude of the site
#'   - `"year"`: Year of the data entry
#'   - `"study_start_year"`: Starting year of the study
#'   - `"study_end_year"`: Ending year of the study
#'   - `"variant_num"`: Numeric variable
#'   - `"total_num"`: Sample size
#'   - `"prev"`: Prevalence value
#'   - `"gene"`: Gene identifier
#'   - `"mut"`: Mutation identifier
#'   - `"variant_string"`: Gene-mutation combination identifier
#'   - `"annotation"`: Annotation details
#'   - `"database"`: Source database
#'   - `"pmid"`: PubMed identifier for source
#'   - `"url"`: URL for additional source information
#'   - `"source"`: Source of the data
#'
#' @return A cleaned data.frame ready for use in STAVE 
#'   (see \url{https://github.com/mrc-ide/STAVE}).
#'
#' @export
deduplicate <- function(df) {
  ### Scenario 1: Identify potential duplicates reported across different study_IDs
  # Group by administrative region, mutation, and collection timeframe
  # Keep groups where more than one unique study_ID reports the same data
  same_data_diff_study <- df %>%
    dplyr::group_by(variant_string, prev, total_num, collection_year_start) %>%
    dplyr::filter(n_distinct(study_ID) > 1) %>%
    dplyr::ungroup()
  
  # Split the dataframe into a list where each list element represents a potential duplicate group
  duplicate_diff_study_list <- same_data_diff_study %>%
    dplyr::group_by(variant_string, prev,  total_num, collection_year_start) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::group_split()
  
  # Check if list of potential duplicates are within 50km radius of one another
  duplicate_diff_study_50km_list <- lapply(duplicate_diff_study_list, filter_duplicates_within_radius)
  duplicate_diff_study_50km_list <- duplicate_diff_study_50km_list[!sapply(duplicate_diff_study_50km_list, is.null)] # Remove NULL elements dataframes that didn't meet the condition)
  
  # Apply tagging function to flag which studies to keep or remove
  tagged_duplicate_diff_study_list <- purrr::map(duplicate_diff_study_50km_list, add_tags_diff_studyID)
  
  # Get summary of duplicates
  summary_diff_study <- summarize_duplicates(tagged_duplicate_diff_study_list)
  
  ### Scenario 2: Identify potential duplicate entries within the same study_ID
  # Group by administrative region, site, mutation, collection timeframe, and study_ID
  # Keep groups with more than one entry (i.e., internal duplication within the study)
  same_data_same_study <- df %>%
    dplyr::group_by(variant_string, prev,  total_num, collection_year_start, study_ID) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
    dplyr::ungroup()
  
  # Split into a list where each list element is a group of internal duplicates
  duplicate_same_study_list <- same_data_same_study %>%
    dplyr::group_by(variant_string, prev,  total_num, collection_year_start, study_ID) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::group_split()
  
  # Check if list of potential duplicates are within 50km radius of one another
  duplicate_same_study_50km_list <- lapply(duplicate_same_study_list, filter_duplicates_within_radius)
  duplicate_same_study_50km_list <- duplicate_same_study_50km_list[!sapply(duplicate_same_study_50km_list, is.null)] # Remove NULL elements dataframes that didn't meet the condition)
  
  # Apply logic to determine which row to keep for internal duplicates
  tagged_duplicate_same_study_list <- lapply(duplicate_same_study_50km_list, add_tags_same_studyID)
  ### TO-DO CECILE: REMOVE EMPTY LIST ITEMS ONCE GEOFF IS FIXED
  tagged_duplicate_same_study_list <- tagged_duplicate_same_study_list[!sapply(tagged_duplicate_same_study_list, is.null)] # Remove null elements from the list (these represent cases where no valid tag was applied)
  length(tagged_duplicate_same_study_list)
  # Get summary of duplicates
  summary_same_study <- summarize_duplicates(tagged_duplicate_same_study_list)
  
  ### Combine duplicates
  # Combine duplicates identified across studies and within studies
  duplicate_list <- c(tagged_duplicate_diff_study_list, tagged_duplicate_same_study_list)
  duplicates_df <- dplyr::bind_rows(duplicate_list) # Convert list of duplicates into one dataframe
  
  # Get all rows from df that are not in the duplicates_df
  df <- df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  duplicates_df <- duplicates_df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  non_duplicate_rows <- df %>%
    dplyr::anti_join(duplicates_df, by = colnames(df)) %>%  # assumes all columns in `df` are relevant for identifying uniqueness
    dplyr::mutate(keep_row = "keep")
  
  # Combine duplicates with non-duplicates into a single dataframe
  deduplicate_df_with_tags <- bind_rows(duplicates_df, non_duplicate_rows)
  
  return(list(df = deduplicate_df_with_tags, 
              summary_same = summary_same_study, 
              summary_diff = summary_diff_study))
}

#' Tag rows in a data.frame for duplicates across different study_IDs
#' 
#' This function tags rows in a data.frame with a "keep_row" column to indicate which 
#' records should be kept or removed based on duplicate study IDs across different databases.
#' The function handles four distinct cases:
#' 
#' - **Case 1:** All records are from the "WHO" or "Pf7k" databases:
#'   - Keep the first study_ID alphabetically and remove the rest.
#' 
#' - **Case 2:** All records are from the "GEOFF" database:
#'   - If duplicates exist between study_IDs starting with "S00" and "S01", keep the study_ID starting with "S01".
#'   - If all study_IDs start with "S00" or "S01", keep the most recent publication based on the `publication_year`.
#' 
#' - **Case 3:** Mixed records from "GEOFF" and either "WHO" or "Pf7k":
#'   - Keep records from "GEOFF" and remove records from "WHO" or "Pf7k".
#' 
#' - **Case 4:** Mixed records from "Pf7k" and "WHO" only:
#'   - Keep records from "Pf7k" and remove records from "WHO".
#' 
#' @param df A `data.frame` containing the data to be processed. It should include the following columns:
#'   - `study_ID`: Unique identifiers for each study.
#'   - `database`: The source database for each study (e.g., "WHO", "Pf7k", "GEOFF").
#'   - `publication_year`: The publication year of each study (used in Case 2).
#' 
#' @return A `data.frame` with an added `keep_row` column, which contains the values "keep" or "remove".
#'   The column indicates whether the corresponding row should be kept or removed based on the criteria for each case.
#'
#' @export
add_tags_diff_studyID <- function(df) {
  # Case 1: All records are from WHO or Pf7k
  if (all(df$database %in% c("WHO", "Pf7k"))) {
    # Identify first study_ID
    first_study_id <- df %>%
      dplyr::arrange(study_ID) %>%
      dplyr::slice(1) %>%
      dplyr::pull(study_ID)
    # Keep the first study_ID (alphabetically) and remove the others
    df <- df %>%
      dplyr::mutate(keep_row = case_when(
        study_ID == first_study_id ~ "keep",
        TRUE ~ "remove"))
  } 
  # Case 2: All records are from GEOFF
  else if (all(df$database == "GEOFF")) {
    df <- df %>%
      dplyr::mutate(keep_row = case_when(
        # If duplicates are between study_IDs begin with S00 and S01, pick S01
        any(grepl("^S00", study_ID)) & any(grepl("^S01", study_ID)) ~ ifelse(grepl("^S01", study_ID), "keep", "remove"),
        # If all study_IDs begin with S00 or S01, pick most recent publication
        all(grepl("^S00|^S01", study_ID)) ~ {
          latest_year <- max(df$publication_year, na.rm = TRUE) # Get the latest publication year
          latest_studies <- df %>%
            filter(publication_year == latest_year) # Subset to those studies
          selected_study <- latest_studies %>%
            slice(1) %>%
            pull(study_ID) # Choose one study arbitrarily (first one)
          ifelse(study_ID == selected_study, "keep", "remove")
        },
        # Default fallback: tag as NA
        TRUE ~ NA_character_
      ))
  }
  # Case 3: Mixture of GEOFF and WHO or Pf7k
  else if (any(df$database == "GEOFF") & any(df$database %in% c("WHO", "Pf7k"))) {
    # Prefer GEOFF data over WHO/Pf7k
    df$keep_row <- ifelse(df$database == "GEOFF", "keep", "remove")
  } 
  # Case 4: Pf7k and WHO only (again), use Pf7k as priority
  else if (any(df$database %in% c("Pf7k", "WHO"))) {
    df <- df %>%
      dplyr::mutate(keep_row = case_when(
        database == "Pf7k" ~ "keep",
        database == "WHO" ~ "remove",
        TRUE ~ NA_character_  # For any other unexpected entries
      ))
  }
  return(df)
}

#' Add tag to rows in data.frame for duplicates within the same studyID
#' 
#' This function tags rows for removal based on duplicate geographic coordinates 
#' within the same study ID. If the geographic coordinates (longitude and latitude) 
#' are the same when rounded to the first decimal place for GEOFF studies, it keeps the 
#' first row and removes the others. If the coordinates are not the same, it returns `NULL`. 
#' For all other studies, the first row is kept, and the others are removed.
#'
#' @param df A data.frame where rows should be tagged as "keep" or "remove".
#' 
#' @return A data.frame with the `keep_row` column indicating if the study should be kept or removed. 
#'   Returns `NULL` if the rounded coordinates for GEOFF studies are not the same.
#'
#' @export
add_tags_same_studyID <- function(df) { 
  # Round longitude and latitude to 1 decimal place
  df$lon_rounded <- round(df$lon, 2)
  df$lat_rounded <- round(df$lat, 2)
  
  # Check for GEOFF studies
  if (all(df$database == "GEOFF")) {
    # Check if the rounded coordinates are the same across all rows
    if (length(unique(paste(df$lon_rounded, df$lat_rounded))) != 1 || length(unique(df$collection_day)) != 1) {
      return(NULL)  # If coordinates are not the same (after rounding), return NULL
    } else {
      # For GEOFF studies with the same coordinates, keep only the first row
      df <- df %>%
        dplyr::mutate(keep_row = ifelse(dplyr::row_number() == 1, "keep", "remove"))
    }
  } else {
    # For all other cases, keep only the first row
    df <- df %>%
      dplyr::mutate(keep_row = ifelse(dplyr::row_number() == 1, "keep", "remove"))
  }
  
  # Remove the rounded columns before returning
  df <- df %>% dplyr::select(-lon_rounded, -lat_rounded)
  
  return(df)
}

#' Filter Dataframe for Duplicates Within a Specified Radius
#' 
#' This function checks whether all geographic coordinates (longitude and latitude) in the given 
#' dataframe are within a specified radius (in kilometers). If any point in the dataframe is outside 
#' the radius, the function returns `NULL`, effectively excluding the dataframe. If all points are within
#' the specified radius, the original dataframe is returned unchanged.
#'
#' @param df A `data.frame` containing at least two columns: 
#'   - `lon`: Longitude values of the geographic points.
#'   - `lat`: Latitude values of the geographic points.
#' 
#' @param radius_km A numeric value specifying the maximum allowed distance (in kilometers) between 
#'   all points. The default value is 50 km. This defines the "radius" within which all points must 
#'   fall in order to keep the dataframe.
#'
#' @return The original `data.frame` if all points are within the specified radius of each other. 
#'   If any point is outside the specified radius, the function returns `NULL`.
#'
#' @export
filter_duplicates_within_radius <- function(df, radius_km = 50) {
  df$lon <- as.numeric(df$lon)
  df$lat <- as.numeric(df$lat)
  coords <- as.matrix(df[, c("lon", "lat")])

  for (i in seq_len(nrow(coords))) {
    # Calculate distances from the current point to all other points
    distances <- sp::spDistsN1(pts = coords, pt = coords[i, ], longlat = TRUE)
    
    # Check if all distances are within the radius (excluding self)
    if (!all(distances <= radius_km | distances == 0)) {
      return(NULL)  # Exclude this dataframe if any point is outside the radius
    }
  }
  
  return(df)  # Keep the dataframe if all points are within the radius
}

#' Summarize Duplicates Across Multiple Dataframes
#'
#' This function takes a list of dataframes containing potential duplicate records and summarizes their similarities and differences.
#' It identifies matching columns, checks for unique values in specific fields, rounds and compares latitude/longitude values, 
#' and determines which study should be kept based on the "keep_row" column.
#'
#' @param duplicate_dfs A list of dataframes, where each dataframe contains potential duplicate entries.
#'                      Each dataframe must have a `study_ID` column and optionally other relevant columns.
#' @param round_digits A numeric value indicating how many decimal places longitude and latitude values should be rounded to
#'                     before checking for duplicates. Default is 1.
#'
#' @return A dataframe summarizing duplicate pairs with the following columns:
#' \item{study_ID_1}{The study ID of the first entry in the duplicate dataframe.}
#' \item{study_ID_2}{The study ID of the second entry in the duplicate dataframe.}
#' \item{matching_columns}{A comma-separated list of columns that have identical values across duplicates.}
#' \item{variant_string, total_num, variant_num, prev, site_name, iso3c}{Columns that maintain their unique value if consistent across duplicates, otherwise marked as "not duplicate".}
#' \item{lon, lat}{Rounded longitude and latitude values if consistent across duplicates, otherwise marked as "not duplicate".}
#' \item{keep}{The study ID that should be retained based on the "keep_row" column. If no entry has "keep", it is set to NA.}
#'
#' @export
summarize_duplicates <- function(duplicate_dfs, round_digits = 1) {
  summary_list <- lapply(duplicate_dfs, function(df) {
    # Extract study IDs (assuming first two columns contain the IDs)
    study_ID_1 <- df$study_ID[1]
    study_ID_2 <- df$study_ID[2]
    
    # Find matching columns where all values are identical
    matching_columns <- names(df)[sapply(df, function(col) length(unique(col)) == 1)]
    
    # Columns to check for uniqueness
    cols_to_check <- c("variant_string", "total_num", "variant_num", "prev", "site_name", "iso3c")
    
    # Extract unique values, set "not duplicate" if multiple unique values exist
    additional_info <- sapply(cols_to_check, function(col) {
      unique_values <- unique(df[[col]])
      if (length(unique_values) == 1) {
        return(unique_values)
      } else {
        return("not duplicate")  # Mark non-duplicate values explicitly
      }
    })
    
    # Special handling for long/lat: round and check uniqueness
    long_lat_info <- sapply(c("lon", "lat"), function(col) {
      rounded_values <- unique(round(df[[col]], round_digits))
      if (length(rounded_values) == 1) {
        return(rounded_values)
      } else {
        return("not duplicate")  # Mark non-duplicate values explicitly
      }
    })
    
    # Determine which study_ID should be marked for "keep"
    keep_value <- ifelse(any(df$keep_row == "keep"), df$study_ID[df$keep_row == "keep"][1], NA)
    
    # Return as a data frame
    data.frame(
      study_ID_1 = study_ID_1,
      study_ID_2 = study_ID_2,
      matching_columns = paste(matching_columns, collapse = ", "),
      as.list(additional_info),  # Add extracted unique values
      as.list(long_lat_info),    # Add rounded long/lat values
      keep = keep_value, 
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all results into a single data frame
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

#' Helper function to identify rows removed by distinct and corresponding kept rows
#'
#' This function identifies duplicate rows in a data frame based on a specified set of columns,
#' returning a list where each element shows the rows that were removed due to duplication and
#' the corresponding row that was kept.
#'
#' @param df A data frame from which duplicates will be identified.
#' @param ... Columns to use for identifying duplicates. These columns are used to define
#'   the criteria for uniqueness.
#'
#' @return A list of data frames, where each element corresponds to a unique combination of the specified columns.
#'   Each data frame in the list contains:
#'   - Rows that were removed as duplicates (`type = "removed"`)
#'   - The row that was kept (`type = "kept"`) for that combination.
#'
#' @details The function groups the data by the specified columns, flags the first occurrence of each unique combination
#'   as the row to keep, and identifies subsequent rows as duplicates. For each duplicate, it finds the corresponding row
#'   that caused the duplication and combines the kept and removed rows into a single output.
#'
#' @examples
#' # Example usage with a data frame `df` having columns `id` and `value`
#' # rows_removed_by_distinct(df, id, value)
#'
#' @export
rows_removed_by_distinct <- function(df, ...) {

  # Capture columns specified in ...
  grouping_cols <- dplyr::ensyms(...)
  
# Identify and keep rows to show removed and retained rows
df_with_flags <- df %>%
  dplyr::group_by(...) %>%
  dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# Rows being kept (distinct rows)
kept_rows <- df_with_flags %>%
  dplyr::filter(.data$keep_row) %>%
  dplyr::select(-.data$keep_row)

# Rows being removed (duplicates)
removed_rows <- df_with_flags %>%
  dplyr::filter(!.data$keep_row) %>%
  dplyr::select(-.data$keep_row)

# Filter `kept_rows` to match unique combinations in `removed_rows`
matching_kept_rows <- kept_rows %>%
  dplyr::semi_join(removed_rows, by = as.character(grouping_cols))  # Only keeps matching rows

# Bind the matching kept rows with removed rows
paired_rows <- dplyr::bind_rows(removed_rows, matching_kept_rows, .id = "type") %>%
  dplyr::mutate(type = ifelse(.data$type == "1", "removed", "kept"))


# Split into a list by the grouping columns to show each set of duplicates separately
paired_list <- paired_rows %>%
  dplyr::group_by(...) %>%
  dplyr::group_split()

return(paired_list)

}
