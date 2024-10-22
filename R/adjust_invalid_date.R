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
  date_fixed <- suppressWarnings(
    case_when(
      grepl("^[0-9]{4}$", date_str) ~ {
        if (is_start) ymd(paste0(date_str, "-01-01")) else ymd(paste0(date_str, "-12-31"))
      },
      grepl("^[0-9]{4}-[0-9]{2}$", date_str) ~ {
        if (is_start) ymd(paste0(date_str, "-01")) else ceiling_date(ymd(paste0(date_str, "-01")), "month") - days(1)
      },
      grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_str) ~ ymd(date_str),
      TRUE ~ as.Date(NA)
    )
  )
  return(date_fixed)
}
