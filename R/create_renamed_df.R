#' Create Dataframe with All Columns from the Original Data with Renamed Cols for WWARN rbind
#'
#' This function creates a dataframe compatible with wwarn cleaned df rbind that includes all columns from the original 
#' `master_table` and adds new columns as specified by the mapping. It preserves the original 
#' columns and appends or replaces columns as needed.
#'
#' @param df A dataframe, typically the `master_table`, containing the original data.
#' @param mapping A named list that maps new column names to existing columns in the dataframe. 
#' If the value is `NA`, the column will be filled with `NA`. If the column is "database", it will 
#' be filled with the `default_database` value.
#' @param default_database A character string specifying the default value for the "database" column. 
#' Defaults to "GEOFF".
#' @return A dataframe that includes all original columns from the input dataframe as well as any additional 
#' columns specified in the mapping.
#' @examples
#' master_table <- data.frame(study_uid = c("s001", "s002"), gene_mutation = c("K13:580Y", "CRT:76T"))
#' mapping <- list(
#'   "gene" = "gene_mutation",
#'   "study" = "study_uid",
#'   "database" = NA
#' )
#' renamed_df <- create_renamed_df(master_table, mapping)
#' print(renamed_df)
#' 
#' 
create_renamed_df <- function(df, mapping, default_database = "GEOFF") {
  result_df <- df  # Start with all columns from the original dataframe
  for (new_col in names(mapping)) {
    if (!is.na(mapping[[new_col]]) && mapping[[new_col]] %in% colnames(df)) {
      result_df[[new_col]] <- df[[mapping[[new_col]]]]
    } else if (new_col == "database") {
      result_df[[new_col]] <- default_database
    } else {
      result_df[[new_col]] <- NA
    }
  }
  return(result_df)
}
