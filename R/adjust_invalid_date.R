#' Adjust Invalid Dates
#'
#' This function adjusts date strings that may be incomplete or invalid by
#' converting them to valid `Date` objects. It supports year-only, year-month,
#' and full year-month-day formats. For incomplete dates, it returns the first
#' or last day of the period, depending on whether the date is a start or end date.
#'
#' @param date_str A string representing a date, which could be in one of the following formats:
#'        "YYYY", "YYYY-MM", or "YYYY-MM-DD".
#' @param is_start A logical value. If `TRUE`, incomplete dates will be adjusted to the
#'        start of the period (e.g., "YYYY" becomes "YYYY-01-01"). If `FALSE`, they will
#'        be adjusted to the end of the period (e.g., "YYYY" becomes "YYYY-12-31").
#' @return A valid `Date` object. If the date is invalid or unrecognized, it returns `NA`.
#' @examples
#' adjust_invalid_date("2021", is_start = TRUE)
#' # Returns: "2021-01-01"
#'
#' adjust_invalid_date("2021-05", is_start = FALSE)
#' # Returns: "2021-05-31"
#'
#' adjust_invalid_date("2021-05-15")
#' # Returns: "2021-05-15"
#'
#' @export
adjust_invalid_date <- function(date_str, is_start = TRUE) {
  date_fixed <- suppressWarnings(dplyr::case_when(
    grepl("^[0-9]{4}$", date_str) ~ {
      if (is_start)
        lubridate::ymd(paste0(date_str, "-01-01"))
      else
        lubridate::ymd(paste0(date_str, "-12-31"))
    },
    grepl("^[0-9]{4}-[0-9]{2}$", date_str) ~ {
      if (is_start)
        lubridate::ymd(paste0(date_str, "-01"))
      else
        lubridate::ceiling_date(lubridate::ymd(paste0(date_str, "-01")), "month") - lubridate::days(1)
    },
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_str) ~ lubridate::ymd(date_str),
    TRUE ~ as.Date(NA)
  ))
  return(date_fixed)
}


#' @noRd
clean_excel_formulas <- function(str) {
  # growing list of known function conversions
  excel_funcs <- list("AVERAGE" = mean, "average" = mean)
  
  # are these patterns in our strings
  fail <- grep(paste0(names(excel_funcs), collapse = "|"), str)
  if (length(fail) > 0) {
    # if so apply the function to the numbers inside the string
    funcs_to_match <- gsub("(=)([A-z]*)(\\(.*)", "\\2", str[fail])
    numbers_to_use <- stringr::str_extract_all(str[fail], "\\d+\\.\\d+")
    results <- vector("numeric", length(numbers_to_use))
    for (i in seq_along(fail)) {
      results[i] <- excel_funcs[[funcs_to_match[i]]](as.numeric(numbers_to_use[[i]]))
    }
    str[fail] <- as.character(results)
  }
  
  return(str)
  
}

#' @noRd
get_column_names_for_clean <- function() {
  c(
    "study_ID",
    "study_name",
    "study_type",
    "authors",
    "publication_year",
    "url",
    "survey_ID",
    "country_name",
    "site_name",
    "lat",
    "lon",
    "spatial_notes",
    "collection_start",
    "collection_end",
    "collection_day",
    "time_notes",
    "variant_string",
    "variant_num",
    "total_num",
    "iso3c",
    "continent",
    "pmid",
    "prev",
    "gene",
    "database"
  )
  
  
}
