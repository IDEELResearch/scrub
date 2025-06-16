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
  
  # ---- Subgroup A: Duplicates among WHO, WWARN and Pf7k ----
  df_who_pf7k_wwarn <- df %>% dplyr::filter(database %in% c("WHO", "Pf7k", "WWARN"))
  
  # Find duplicates by variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end
  dup_who_pf7k_wwarn_list <- df_who_pf7k_wwarn %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end) %>%
    dplyr::filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # Filter by spatial proximity (e.g., within 50km)
  dup_who_pf7k_wwarn_filtered <- lapply(dup_who_pf7k_wwarn_list, filter_duplicates_within_radius)
  dup_who_pf7k_wwarn_filtered <- dup_who_pf7k_wwarn_filtered[!sapply(dup_who_pf7k_wwarn_filtered, is.null)]
  
  # ---- Subgroup B: Duplicates between GEOFF and WHO/Pf7k/WWARN ----
  df_geoff <- df %>% filter(database == "GEOFF")
  
  # Find GEOFF and WHO, WWARN or Pf7k duplicates by collection_year_start
  geoff_matches_start <- df_geoff %>%
    semi_join(df_who_pf7k_wwarn, by = c("iso3c", "gene", "collection_year_start"))
  who_pf7k_wwarn_matches_start <- df_who_pf7k_wwarn %>%
    semi_join(df_geoff, by = c("iso3c", "gene", "collection_year_start"))
  dup_geoff_who_pf7k_wwarn_start <- bind_rows(geoff_matches_start, who_pf7k_wwarn_matches_start) %>%
    group_by(iso3c, gene, collection_year_start) %>%
    dplyr::filter(n() > 1) %>%
    group_split()
  
  # Find GEOFF and WHO, WWARN or Pf7k duplicates by collection_year_end
  geoff_matches_end <- df_geoff %>%
    semi_join(df_who_pf7k_wwarn, by = c("iso3c", "gene", "collection_year_end"))
  who_pf7kk_wwarn_matches_end <- df_who_pf7k_wwarn %>%
    semi_join(df_geoff, by = c("iso3c", "gene", "collection_year_end"))
  dup_geoff_who_pf7k_wwarn_end <- bind_rows(geoff_matches_end, who_pf7kk_wwarn_matches_end) %>%
    group_by(iso3c, gene, collection_year_end) %>%
    dplyr::filter(n() > 1) %>%
    group_split()
  
  # Filter by spatial proximity (e.g., within 50km)
  dup_geoff_who_pf7k_wwarn_filtered_start <- lapply(dup_geoff_who_pf7k_wwarn_start, filter_duplicates_within_radius)
  dup_geoff_who_pf7k_wwarn_filtered_start <- dup_geoff_who_pf7k_wwarn_filtered_start[!sapply(dup_geoff_who_pf7k_wwarn_filtered_start, is.null)]
  dup_geoff_who_pf7k_wwarn_filtered_end <- lapply(dup_geoff_who_pf7k_wwarn_end, filter_duplicates_within_radius)
  dup_geoff_who_pf7k_wwarn_filtered_end <- dup_geoff_who_pf7k_wwarn_filtered_end [!sapply(dup_geoff_who_pf7k_wwarn_filtered_end, is.null)]
  
  # Keep only unique groups
  dup_geoff_who_pf7k_wwarn_filtered <- deduplicate_list(
    c(dup_geoff_who_pf7k_wwarn_filtered_start, dup_geoff_who_pf7k_wwarn_filtered_end),
    group_signature_geoff_who_pf7k_wwarn
  )
  
  # ---- Subgroup C: Duplicates within GEOFF ----
  # Find GEOFF duplicates based on variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end
  dup_geoff <- df_geoff %>%
    group_by(variant_string, prev, total_num, variant_num, collection_year_start, collection_year_end) %>%
    dplyr::filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # Filter by spatial proximity (e.g., within 50km)
  dup_geoff_filtered <- lapply(dup_geoff, filter_duplicates_within_radius)
  dup_geoff_filtered <- dup_geoff_filtered[!sapply(dup_geoff_filtered, is.null)]
  
  # Combine all subgroup duplicates
  duplicate_diff_studies <- c(dup_who_pf7k_wwarn_filtered, dup_geoff_who_pf7k_wwarn_filtered, dup_geoff_filtered)
  
  # Tag which studies to keep or remove
  tagged_duplicate_diff_study_list <- lapply(duplicate_diff_studies, add_tags_diff_studyID)
  tagged_duplicate_diff_study_list <- tagged_duplicate_diff_study_list[!sapply(tagged_duplicate_diff_study_list, is.null)]
  
  ### Scenario 2: Identify potential duplicate entries within the same study_ID
  # Group by administrative region, site, mutation, collection timeframe, and study_ID
  # Keep groups with more than one entry (i.e., internal duplication within the study)
  duplicate_same_study_list <- df %>%
    dplyr::group_by(variant_string, prev, total_num, variant_num, study_ID, collection_year_start, collection_year_end) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::mutate(keep_row = dplyr::row_number() == 1) %>%
    dplyr::group_split()
  
  # Filter by spatial proximity (e.g., within 50km)
  duplicate_same_study_50km_list <- lapply(duplicate_same_study_list, filter_duplicates_within_radius)
  duplicate_same_study_50km_list <- duplicate_same_study_50km_list[!sapply(duplicate_same_study_50km_list, is.null)] # Remove NULL elements dataframes that didn't meet the condition)
  
  # Tag which studies to keep or remove
  tagged_duplicate_same_study_list <- lapply(duplicate_same_study_50km_list, add_tags_same_studyID)
  tagged_duplicate_same_study_list <- tagged_duplicate_same_study_list[!sapply(tagged_duplicate_same_study_list, is.null)]
  
  # ---- Combine duplicates between different study IDs and within the same study ID ----
  # Combine duplicates identified across studies and within studies
  duplicate_list <- c(tagged_duplicate_diff_study_list, tagged_duplicate_same_study_list)
  duplicates_df <- dplyr::bind_rows(duplicate_list) # Convert list of duplicates into one dataframe
  
  # Get all rows from df that are not in the duplicates_df
  df <- df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  duplicates_df <- duplicates_df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  cols_for_join <- setdiff(colnames(df), "keep_row")
  
  # Identify datapoints that are no duplicates
  non_duplicate_rows <- df %>%
    dplyr::anti_join(duplicates_df %>% select(all_of(cols_for_join)), by = cols_for_join) %>%
    dplyr::mutate(keep_row = "keep")
  
  # Combine duplicates with non-duplicates into a single dataframe
  deduplicate_df_with_tags_pre <- bind_rows(duplicates_df, non_duplicate_rows)
  deduplicate_df_with_tags <- deduplicate_df_with_tags_pre %>% distinct()
  
  return(list(df = deduplicate_df_with_tags,
              list_same = tagged_duplicate_same_study_list, 
              list_diff = tagged_duplicate_diff_study_list))
}


#' Deduplicate a List of Data Frames Based on Group Signatures
#'
#' Iterates through a list of data frames and retains only the first occurrence of each
#' unique group based on a provided signature-generating function. This helps eliminate 
#' duplicate groupings across overlapping criteria (e.g., collection year start and end).
#'
#' @param list_of_dfs A list of data frames to deduplicate. Each data frame represents a group.
#' @param sig_func A function that takes a data frame and returns a unique string signature
#'        representing the group (e.g., via concatenation of key identifying columns).
#'
#' @return A deduplicated list of data frames, where each entry corresponds to a unique signature.
#' 
#' @export
deduplicate_list <- function(list_of_dfs, sig_func) {
  seen_sigs <- character()
  out <- list()
  for (df in list_of_dfs) {
    sig <- sig_func(df)
    if (!(sig %in% seen_sigs)) {
      seen_sigs <<- c(seen_sigs, sig)
      out <- append(out, list(df))
    }
  }
  return(out)
}

#' Deduplicate Studies by PMID with Conflicting Study Identifiers
#'
#' Identifies groups of records that share the same PubMed ID (`pmid`) but have different 
#' `study_ID`s, and applies a deduplication strategy to tag each row as `"keep"` or `"remove"`. 
#' The deduplication logic is defined in the helper function `add_tags_diff_studyID()`, 
#' which applies a set of prioritized rules based on the `database` source.
#'
#' @param df A data frame containing at minimum the following columns:
#'   - `pmid`: PubMed identifier
#'   - `study_ID`: Study identifier
#'   - `database`: Data source (e.g., WHO, GEOFF, Pf7k, WWARN)
#'   - Any additional fields required by `add_tags_diff_studyID()` (e.g., `publication_year`)
#'
#' @return A list of data frames, each corresponding to a set of records with the same 
#' `pmid` but different `study_ID`s. Each data frame includes a `keep_row` column 
#' indicating which rows to retain (`"keep"`) or discard (`"remove"`).
#'
#' @export
deduplicate_pmid <- function(df){
  # Identify groups with duplicate pmids but differing study_IDs
  duplicate_pmid_list <- df %>%
    dplyr::filter(!is.na(pmid)) %>%
    group_by(pmid) %>%
    dplyr::filter(n_distinct(study_ID) > 1) %>%
    group_split()
  
  # Apply tagging logic to each group
  tagged_deduplicate_pmid_list <- lapply(duplicate_pmid_list, add_tags_diff_studyID)
  
  # Remove the duplicates with Pf7k since they all have the same pmid
  deduplicate_pmid_list_filtered <- tagged_deduplicate_pmid_list[
    !sapply(tagged_deduplicate_pmid_list, function(df) all(df$database == "Pf7k"))]
  
  # Combine tagged duplicates into a single dataframe
  tagged_duplicates_df <- dplyr::bind_rows(deduplicate_pmid_list_filtered)
  
  # Ensure numeric comparison compatibility for coordinates
  df <- df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  tagged_duplicates_df <- tagged_duplicates_df %>% dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))
  
  # Identify all rows not involved in deduplication
  non_duplicate_rows <- df %>%
    dplyr::anti_join(tagged_duplicates_df, by = colnames(df)) %>%
    dplyr::mutate(keep_row = NA)
  
  # Combine tagged duplicates and untouched rows into one complete dataframe
  deduplicate_df_with_tags <- bind_rows(tagged_duplicates_df, non_duplicate_rows)
  
  return(list(
    tagged_list = deduplicate_pmid_list_filtered,
    deduplicated_df = deduplicate_df_with_tags
  ))
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
group_signature_geoff_who_pf7k_wwarn <- function(df) {
  df %>%
    select(iso3c, gene, collection_year_start, collection_year_end) %>%
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
  # Case 1: Mixture of GEOFF and WHO, WWARN or Pf7k
  if (any(df$database == "GEOFF") & any(df$database %in% c("WHO", "WWARN", "Pf7k"))) {
    df <- df %>% mutate(keep_row = if_else(database == "GEOFF", "keep", "remove"))
  } 
  # Case 2: Pf7k, WWARN and WHO only, keep Pf7k as priority
  else if (setequal(unique(df$database), c("Pf7k", "WHO", "WWARN"))) {
    df <- df %>%
      mutate(keep_row = if_else(database == "Pf7k", "keep", "remove"))
  } 
  # Case 3: WWARN and WHO only, keep WWARN as priority
  else if (setequal(unique(df$database), c("WHO", "WWARN"))) {
    df <- df %>%
      mutate(keep_row = if_else(database == "WWARN", "keep", "remove"))
  }
  # Case 4: All WHO, WWARN, or Pf7k pick one random study_ID to keep
  else if (all(df$database == "WHO") || all(df$database == "WWARN") || all(df$database == "Pf7k")) {
    set.seed(42)  # for reproducibility
    one_study <- sample(unique(df$study_ID), 1)
    df <- df %>%
      mutate(keep_row = if_else(study_ID == one_study, "keep", "remove"))
  } 
  # Case 5: All records are from GEOFF, pick study_ID with geoff_S01 over geoff_S00, or keep most recent publication
  else if (all(df$database == "GEOFF")) {
    study_ids <- unique(df$study_ID)
    # If mixture of geoff_S00 and geoff_S01 → keep all S01s
    if (any(grepl("^geoff_S00", study_ids)) && any(grepl("^geoff_S01", study_ids))) {
      df$keep_row <- ifelse(grepl("^geoff_S01", df$study_ID), "keep", "remove")
      # If all geoff_S00 or all geoff_S01 → keep one study with latest year (randomly if tied)
    } else if (all(grepl("^geoff_S00|^geoff_S01", study_ids))) {
      latest_year <- max(df$publication_year, na.rm = TRUE)
      candidates <- df %>%
        filter(publication_year == latest_year) %>%
        distinct(study_ID)
      set.seed(42)  # for reproducibility
      selected_study <- sample(candidates$study_ID, 1)
      df$keep_row <- ifelse(df$study_ID == selected_study, "keep", "remove")
      # Fallback if no geoff_S00/geoff_S01 pattern
    } else {
      df$keep_row <- NA_character_
    }
  }
  # Case 6: Fallback assign NA
  else {
    df <- df %>% mutate(keep_row = NA_character_)
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
