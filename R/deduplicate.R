#' Identify and summarize potential duplicate studies across multiple databases
#'
#' This function identifies potential duplicates among WHO, Pf7k, and GEOFF datasets,
#' groups them by relevant variables, filters by spatial proximity, tags which to keep/remove,
#' and summarizes the results.
#'
#' @param df A data.frame to be deduplicated. It must contain the following columns:
#'   - `"study_ID"`: Study ID
#'   - `"study_name"`: Name of the study
#'   - `"study_type"`: Type of study conducted
#'   - `"authors"`: Authors of the study
#'   - `"publication_year"`: Year the study was published
#'   - `"url"`: URL for additional source information
#'   - `"survey_ID"`: Survey ID
#'   - `"country_name"`: Name of the country
#'   - `"site_name"`: Site name or identifier
#'   - `"lat"`: Latitude of the site
#'   - `"lon"`: Longitude of the site
#'   - `"spatial_notes"`: Notes on spatial data
#'   - `"collection_start"`: Start date of sample collection
#'   - `"collection_end"`: End date of sample collection
#'   - `"collection_day"`: Day of sample collection
#'   - `"time_notes"`: Notes on time or collection period
#'   - `"variant_string"`: Gene-mutation combination identifier
#'   - `"variant_num"`: Numeric representation of the variant
#'   - `"total_num"`: Sample size
#'   - `"iso3c"`: ISO3 country code
#'   - `"continent"`: Continent of the study site
#'   - `"pmid"`: PubMed ID for the source
#'   - `"prev"`: Prevalence value
#'   - `"gene"`: Gene identifier
#'   - `"database"`: Source database
#'   - `"collection_year_start"`: Starting year of the sample collection
#'   - `"collection_year_end"`: Ending year of the sample collection
#'
#' @return A cleaned data.frame ready for use in STAVE 
#'   (see \url{https://github.com/mrc-ide/STAVE}).
#'
#' @export
deduplicate <- function(df) {
  ### Scenario 1: Identify potential duplicates reported across different study_IDs
  # Group by administrative region, mutation, and collection timeframe
  # Keep groups where more than one unique study_ID reports the same data
  
  # ---- Subgroup A: Duplicates among WHO and Pf7k ----
  df_who_pf7k <- df %>% filter(database %in% c("WHO", "Pf7k"))
  
  # Find duplicates by collection_year_start
  dup_who_pf7k_start <- df_who_pf7k %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_start) %>%
    filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # Find duplicates by collection_year_end
  dup_who_pf7k_end <- df_who_pf7k %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_end) %>%
    filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # ---- Subgroup B: Duplicates between GEOFF and WHO/Pf7k ----
  df_geoff <- df %>% filter(database == "GEOFF")
  
  # Find GEOFF/WHO-Pf7k duplicates by collection_year_start
  geoff_matches_start <- df_geoff %>%
    semi_join(df_who_pf7k, by = c("iso3c", "gene", "collection_year_start"))
  who_pf7k_matches_start <- df_who_pf7k %>%
    semi_join(df_geoff, by = c("iso3c", "gene", "collection_year_start"))
  dup_geoff_who_pf7k_start <- bind_rows(geoff_matches_start, who_pf7k_matches_start) %>%
    group_by(iso3c, gene, collection_year_start) %>%
    filter(n() > 1) %>%
    group_split()
  
  # Find GEOFF/WHO-Pf7k duplicates by collection_year_end
  geoff_matches_end <- df_geoff %>%
    semi_join(df_who_pf7k, by = c("iso3c", "gene", "collection_year_end"))
  who_pf7k_matches_end <- df_who_pf7k %>%
    semi_join(df_geoff, by = c("iso3c", "gene", "collection_year_end"))
  dup_geoff_who_pf7k_end <- bind_rows(geoff_matches_end, who_pf7k_matches_end) %>%
    group_by(iso3c, gene, collection_year_end) %>%
    filter(n() > 1) %>%
    group_split()
  
  # ---- Subgroup C: Duplicates within GEOFF ----
  # Find GEOFF duplicates by collection_year_start
  dup_geoff_start <- df_geoff %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_start) %>%
    filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # Find GEOFF duplicates by collection_year_end
  dup_geoff_end <- df_geoff %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_end) %>%
    filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # ---- Generate signatures for each group ----
  sigs_who_pf7k_start <- map_chr(dup_who_pf7k_start, group_signature_who_pf7k)
  sigs_who_pf7k_end  <- map_chr(dup_who_pf7k_end, group_signature_who_pf7k)
  sigs_geoff_who_pf7k_start <- map_chr(dup_geoff_who_pf7k_start, group_signature_geoff_who_pf7k)
  sigs_geoff_who_pf7k_end <- map_chr(dup_geoff_who_pf7k_end, group_signature_geoff_who_pf7k)
  sigs_geoff_start <- map_chr(dup_geoff_start, group_signature_geoff)
  sigs_geoff_end <- map_chr(dup_geoff_end, group_signature_geoff)
  
  # ---- Combine all groups and remove duplicates based on signature ----
  all_groups_diff_study <- c(dup_who_pf7k_start, dup_who_pf7k_end, dup_geoff_who_pf7k_start, dup_geoff_who_pf7k_end, dup_geoff_start, dup_geoff_end)
  all_sigs_diff_study <- c(sigs_who_pf7k_start, sigs_who_pf7k_end, sigs_geoff_who_pf7k_start, sigs_geoff_who_pf7k_end, sigs_geoff_start, sigs_geoff_end)
  
  # Keep only unique groups
  diff_data_diff_study_list <- all_groups_diff_study[!duplicated(all_sigs_diff_study)]
  
  # ---- Filter by spatial proximity (e.g., within 50km) ----
  duplicate_diff_study_50km_list <- lapply(diff_data_diff_study_list, filter_duplicates_within_radius)
  duplicate_diff_study_50km_list <- duplicate_diff_study_50km_list[!sapply(duplicate_diff_study_50km_list, is.null)] # Remove NULL elements dataframes that didn't meet the condition)
  
  # ---- Tag which studies to keep or remove ----
  tagged_duplicate_diff_study_list <- lapply(duplicate_diff_study_50km_list, add_tags_diff_studyID)
  tagged_duplicate_diff_study_list <- tagged_duplicate_diff_study_list[!sapply(tagged_duplicate_diff_study_list, is.null)]
  
  ### Scenario 2: Identify potential duplicate entries within the same study_ID
  # Group by administrative region, site, mutation, collection timeframe, and study_ID
  # Keep groups with more than one entry (i.e., internal duplication within the study)
  same_data_same_study_start <- df %>%
    dplyr::group_by(variant_string, prev, total_num, variant_num, study_ID, collection_year_start) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
    dplyr::group_split()
  
  same_data_same_study_end <- df %>%
    dplyr::group_by(variant_string, prev, total_num, variant_num, study_ID, collection_year_end) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
    dplyr::group_split()
  
  # ---- Generate signatures for each group ----
  sigs_same_study_start <- map_chr(same_data_same_study_start, group_signature_same_study)
  sigs_same_study_end  <- map_chr(same_data_same_study_end, group_signature_same_study)
  
  # ---- Combine all groups and remove duplicates based on signature ----
  all_groups_same_study <- c(same_data_same_study_start, same_data_same_study_end)
  all_sigs_same_study   <- c(sigs_same_study_start, sigs_same_study_end)
  
  duplicate_same_study_list <- all_groups_same_study[!duplicated(all_sigs_same_study)]
  
  # ---- Filter by spatial proximity (e.g., within 50km) ----
  duplicate_same_study_50km_list <- lapply(duplicate_same_study_list, filter_duplicates_within_radius)
  duplicate_same_study_50km_list <- duplicate_same_study_50km_list[!sapply(duplicate_same_study_50km_list, is.null)] # Remove NULL elements dataframes that didn't meet the condition)
  
  # ---- Tag which studies to keep or remove ----
  tagged_duplicate_same_study_list <- lapply(duplicate_same_study_50km_list, add_tags_same_studyID)
  tagged_duplicate_same_study_list <- tagged_duplicate_same_study_list[!sapply(tagged_duplicate_same_study_list, is.null)]
  
  # ---- Combine duplicates between different study IDs and within the same study ID ----
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
              list_same = duplicate_same_study_50km_list, 
              list_diff = duplicate_diff_study_50km_list))
}

#' Generate a Group Signature for WHO/Pf7k Data
#'
#' Creates a unique signature string for a WHO or Pf7k dataset by concatenating key fields 
#' (variant information, prevalence, sample sizes, and collection years). This is useful 
#' for identifying and comparing groups of observations in deduplication tasks.
#'
#' @param df A data frame containing the following columns:
#'   - `variant_string`
#'   - `prev`
#'   - `total_num`
#'   - `variant_num`
#'   - `collection_year_start`
#'   - `collection_year_end`
#'
#' @return A single string that uniquely identifies the group, with values concatenated 
#' using underscores and records separated by `|`.
#' @export
group_signature_who_pf7k <- function(df) {
  df %>%
    select(variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end) %>%
    distinct() %>%
    arrange(across(everything())) %>%
    mutate_all(as.character) %>%
    unite("sig", everything(), sep = "_", remove = TRUE) %>%
    pull(sig) %>%
    paste(collapse = "|")
}

#' Generate a Group Signature for GEOFF vs WHO/Pf7k Comparison
#'
#' Creates a unique signature string to compare a GEOFF dataset to WHO or Pf7k datasets.
#' Uses either `collection_year_start` or `collection_year_end`, depending on which is 
#' available. Useful in deduplication and data harmonization.
#'
#' @param df A data frame containing the following columns:
#'   - `variant_string`
#'   - `collection_year_start`
#'   - `collection_year_end`
#'
#' @return A single string representing the group signature, with values concatenated 
#' using underscores and records separated by `|`.
#' @export
group_signature_geoff_who_pf7k <- function(df) {
  # Use the year that is not all NA
  if (all(is.na(df$collection_year_start))) {
    df %>%
      select(variant_string, collection_year_end) %>%
      distinct() %>%
      arrange(across(everything())) %>%
      mutate_all(as.character) %>%
      unite("sig", everything(), sep = "_", remove = TRUE) %>%
      pull(sig) %>%
      paste(collapse = "|")
  } else {
    df %>%
      select(variant_string, collection_year_start) %>%
      distinct() %>%
      arrange(across(everything())) %>%
      mutate_all(as.character) %>%
      unite("sig", everything(), sep = "_", remove = TRUE) %>%
      pull(sig) %>%
      paste(collapse = "|")
  }
}

#' Generate a Group Signature for GEOFF Data
#'
#' Creates a unique signature string for GEOFF datasets based on variant and prevalence data,
#' as well as sample sizes and collection years. This aids in identifying and distinguishing 
#' unique groups in deduplication workflows.
#'
#' @param df A data frame containing the following columns:
#'   - `variant_string`
#'   - `prev`
#'   - `total_num`
#'   - `variant_num`
#'   - `collection_year_start`
#'   - `collection_year_end`
#'
#' @return A single concatenated string signature for the group.
#' @export
group_signature_geoff <- function(df) {
  df %>%
    select(variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end) %>%
    distinct() %>%
    arrange(across(everything())) %>%
    mutate_all(as.character) %>%
    unite("sig", everything(), sep = "_", remove = TRUE) %>%
    pull(sig) %>%
    paste(collapse = "|")
}

#' Generate Study-Specific Signature for Deduplication
#'
#' Creates a unique string signature for a group of entries originating from the same study.
#' This is used for identifying and deduplicating records that share the same `study_ID` and
#' other key identifying variables such as variant information, prevalence, and collection years.
#'
#' @param df A data frame containing at least the following columns:
#'   - `"study_ID"`: Unique study identifier
#'   - `"variant_string"`: Identifier combining gene and mutation
#'   - `"prev"`: Prevalence estimate
#'   - `"total_num"`: Total sample size
#'   - `"variant_num"`: Number of samples with the variant
#'   - `"collection_year_start"`: Start year of data collection
#'   - `"collection_year_end"`: End year of data collection
#'
#' @return A single character string that concatenates all unique combinations of
#' the input fields, sorted and joined by underscores, and collapsed using `"|"`.
#'
#' @export
group_signature_same_study <- function(df){
  df %>%
    select(study_ID, variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end) %>%
    distinct() %>%
    arrange(across(everything())) %>%
    mutate_all(as.character) %>%
    unite("sig", everything(), sep = "_", remove = TRUE) %>%
    pull(sig) %>%
    paste(collapse = "|")
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
  # Check for GEOFF studies
  if (all(df$database == "GEOFF")) {
    # Check if the rounded coordinates are the same across all rows
    return (NULL)
    }
  
  # Remove the rounded columns before returning
  df <- df %>% 
    dplyr::mutate(keep_row = if_else(row_number() == 1, "keep", "remove"))
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
