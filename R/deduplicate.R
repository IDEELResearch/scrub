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
    dplyr::group_by(name_2, variant_string, prev, total_num, collection_year) %>%
    dplyr::filter(n_distinct(study_ID) > 1) %>%
    dplyr::ungroup()
  
  # Split the dataframe into a list where each list element represents a potential duplicate group
  duplicate_diff_study_list <- same_data_diff_study %>%
    dplyr::group_by(name_2, variant_string, prev,  total_num, collection_year) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::group_split()
  
  # Apply tagging function to flag which studies to keep or remove
  tagged_duplicate_diff_study_list <- purrr::map(duplicate_diff_study_list, add_tags_diff_studyID)
  
  ### Scenario 2: Identify potential duplicate entries within the same study_ID
  # Group by administrative region, site, mutation, collection timeframe, and study_ID
  # Keep groups with more than one entry (i.e., internal duplication within the study)
  same_data_same_study <- df %>%
    dplyr::group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
    dplyr::ungroup()
  
  # Split into a list where each list element is a group of internal duplicates
  duplicate_same_study_list <- same_data_same_study %>%
    dplyr::group_by(name_2, site_name, variant_string, collection_start, collection_end, study_ID) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::group_split()
  
  # Apply logic to determine which row to keep for internal duplicates
  tagged_duplicate_same_study_list <- lapply(duplicate_same_study_list, handle_same_studyID_duplicates)
  
  ### TO-DO CECILE: REMOVE EMPTY LIST ITEMS ONCE GEOFF IS FIXED
  # Remove null elements from the list (these represent cases where no valid tag was applied)
  tagged_duplicate_same_study_list <- tagged_duplicate_same_study_list[!sapply(tagged_duplicate_same_study_list, is.null)]
  
  # Combine duplicates identified across studies and within studies
  duplicate_list <- c(tagged_duplicate_diff_study_list, tagged_duplicate_same_study_list)
  
  # Convert list of duplicates into one dataframe
  duplicates_df <- dplyr::bind_rows(duplicate_list)
  
  # Get all rows from df that are not in the duplicates_df
  non_duplicate_rows <- df %>%
    dplyr::anti_join(duplicates_df, by = colnames(df)) %>%  # assumes all columns in `df` are relevant for identifying uniqueness
    dplyr::mutate(keep_row = "keep")
  
  # Combine duplicates with non-duplicates into a single dataframe
  deduplicate_df_with_tags <- bind_rows(duplicates_df, non_duplicate_rows)
  
  return(deduplicate_df_with_tags)
}

#' Spatial join for admin2 regions
#'
#' This function spatially joins the each survey_IDs (row) reported longitude 
#' and latitude with the admin2 region reported by geodata.
#'
#' @param df A data.frame to be . It must contain the following columns:
#'   - `"iso3c"`: ISO3 country code
#'   - `"lat"`: Latitude of the site
#'   - `"long"`: Longitude of the site
#'
#' @return A spatially joined data.frame with the corresponding admin2 region
#' per row.
#'
#' @export
sf_join_admn2 <- function(df, admin2_sf) {
  df_sf <- df %>%
    filter(!is.na(lon) & !is.na(lat)) %>% # Delete rows where lat/long is NA - NOTE: study S0147Jeang2024 is missing lat
    mutate(lon_keep = lon, lat_keep = lat) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # Ensure lat lon coordinates use CRS 4326
    sf::st_transform(crs = sf::st_crs(admin2_sf)) # Ensure both dataframes use the same CRS (WGS 84)
  
  # Perform spatial join: Match points (df_sf) to Admin 2 polygons (admin2_sf)
  sf::sf_use_s2(FALSE)
  # 1. Perform initial spatial join
  df_sf_joined <- sf::st_join(df_sf, admin2_sf, left = TRUE)
  
  # 2. Identify unmatched points
  unmatched <- df_sf_joined %>% filter(is.na(name_2))
  
  # 3. Get nearest polygon for unmatched points
  if (nrow(unmatched) > 0) {
    nearest_idx <- sf::st_nearest_feature(unmatched, admin2_sf)
    nearest_matches <- admin2_sf[nearest_idx, ]
    
    # Compare iso codes and keep name_2 only if they match
    name_2_corrected <- ifelse(
      unmatched$iso == nearest_matches$admin2_iso3c,
      nearest_matches$name_2,
      NA_character_
    )
    
    # Update admin2_iso3c accordingly (optional, NA if iso doesn't match)
    admin2_iso3c_corrected <- ifelse(
      unmatched$iso == nearest_matches$admin2_iso3c,
      nearest_matches$admin2_iso3c,
      NA_character_
    )
    
    # Bind corrected columns and retain original geometry
    unmatched_filled <- unmatched %>%
      select(-name_2, -admin2_iso3c) %>%
      mutate(
        name_2 = name_2_corrected,
        admin2_iso3c = admin2_iso3c_corrected
      )
    
    # 4. Combine matched and filled unmatched results
    matched <- df_sf_joined %>% filter(!is.na(name_2))
    df_sf_joined_final <- bind_rows(matched, unmatched_filled)
    
  } else {
    df_sf_joined_final <- df_sf_joined
  }
  
  # Convert back to a regular dataframe
  df_final <- df_sf_joined %>% 
    sf::st_drop_geometry() %>% 
    mutate(lon = lon_keep, lat = lat_keep) %>%  # Restore lat/lon
    select(-lon_keep, -lat_keep)  # Remove temporary columns
  
  return(df_final)
}

#' Add tag to rows in data.frame for duplicates between different studies
#' 
#' This function handles cases where duplicates exist across different study_IDs
#' by adding a column keep_row that determines which rows should be keep and 
#' removed. There are 4 cases:
#'    - Case 1: All records are from WHO or Pf7k: keep the first study_ID 
#'      (alphabetically) and remove the others
#'    - Case 2: All records are from GEOFF: If duplicates are between study_IDs begin 
#'      with S00 and S01, pick S01 and if all study_IDs begin with S00 or S01, 
#'      pick most recent publication
#'    - Case 3: Mixture of GEOFF and WHO or Pf7k: Keep GEOFF data over WHO/Pf7k
#'    - Case 4: Mixture of Pf7k and WHO: keep Pf7k and remove WHO
#'
#' @param df A data.frame where rows should be tagged as "keep" or "remove".
#' 
#' @return A data.frame with the keep_row column indicating if the study should 
#' be kept or removed. 
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

#' Add tag to rows in data.frame for duplicates within the same study
#' 
#' This function adds column keep_row and marks the first the first row in the group
#' to keep and remove the other rows which are duplicates.
#'
#' @param df A data.frame where rows should be tagged as "keep" or "remove".
#' 
#' @return A data.frame with the keep_row column indicating if the study should 
#' be kept or removed. 
#'
#' @export
handle_same_studyID_duplicates <- function(df) { 
  # Check if all rows in the group come from the GEOFF database
  if (all(df$database == "GEOFF")) {
    # If all rows are GEOFF, we currently exclude these from processing
    # This is a temporary rule until GEOFF duplicates are resolved
    return(NULL)
  } else {
    # For all other databases, keep only the first row in the group
    # This helps retain a representative entry and remove internal duplicates
    df %>%
      dplyr::mutate(keep_row = ifelse(dplyr::row_number() == 1, "keep", "remove"))
  }
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
