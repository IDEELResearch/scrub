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
#'   - `"x"`: Numeric variable
#'   - `"n"`: Sample size
#'   - `"prev"`: Prevalence value
#'   - `"gene"`: Gene identifier
#'   - `"mut"`: Mutation identifier
#'   - `"gene_mut"`: Gene-mutation combination identifier
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
#'
deduplicate <- function(df) {
  
  df %>% 
    # fully remove complete duplicates at all columns
    dplyr::distinct() %>% 
    # fully remove complete duplicates at all columns minus database
    dplyr::distinct(dplyr::across(-database))
  
  return(df)
  
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
  grouping_cols <- ensyms(...)
  
# Identify and keep rows to show removed and retained rows
df_with_flags <- df %>%
  group_by(...) %>%
  mutate(keep_row = row_number() == 1) %>%
  ungroup()

# Rows being kept (distinct rows)
kept_rows <- df_with_flags %>%
  filter(keep_row) %>%
  select(-keep_row)

# Rows being removed (duplicates)
removed_rows <- df_with_flags %>%
  filter(!keep_row) %>%
  select(-keep_row)

# Filter `kept_rows` to match unique combinations in `removed_rows`
matching_kept_rows <- kept_rows %>%
  semi_join(removed_rows, by = as.character(grouping_cols))  # Only keeps matching rows

# Bind the matching kept rows with removed rows
paired_rows <- bind_rows(removed_rows, matching_kept_rows, .id = "type") %>%
  mutate(type = ifelse(type == "1", "removed", "kept"))


# Split into a list by the grouping columns to show each set of duplicates separately
paired_list <- paired_rows %>%
  group_by(...) %>%
  group_split()

return(paired_list)

}
