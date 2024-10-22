#' Check and Validate Substudy Entries
#'
#' This function checks the `substudy` column for invalid entries by first correcting known bad strings 
#' and then ensuring all entries match the allowed categories or follow the "day X extracted/calculated" pattern. 
#' It prints any invalid entries found, or a message confirming that all entries are valid after corrections.
#'
#' @param substudy_column A character vector containing substudy entries that may need validation and corrections.
#' @return A character vector with corrected and validated substudy entries. Prints invalid entries if found.
#' @examples
#' substudy_data <- c("day 0", "untreated_extracted", "untreadextracted", "day3calculated", "day 5 extracted")
#' validated_substudy <- check_substudy_entries(substudy_data)
#' # Output:
#' # Corrected 1 instances of day0extracted
#' # Corrected 1 instances of untreated_extracted
#' # Corrected 1 instances of untreadextracted
#' # Invalid entries found in 'substudy':
#' # [1] "day 5 extracted"
#'
#' @export
check_substudy_entries <- function(substudy_column) {
  # Correct bad strings first and track corrections
  substudy_column <- correct_substudy_entries(substudy_column)
  
  # Define allowed categories
  allowed_categories <- c(
    "untreatedextracted", "untreatedcalculated", "day0extracted", "day0calculated",
    "day3calculated", "day24calculated", "treatedextracted", "treatedcalculated"
  )
  
  # Regular expression for "day X extracted" or "day X calculated"
  day_extracted_calculated_pattern <- "^day[0-9]+(extracted|calculated)$"
  
  # Find any entries not matching the allowed categories or the day X extracted/calculated pattern
  invalid_entries <- substudy_column[!substudy_column %in% allowed_categories & 
                                       !grepl(day_extracted_calculated_pattern, substudy_column)]
  
  # If there are invalid entries, print them
  if (length(invalid_entries) > 0) {
    message("Invalid entries found in 'substudy':")
    print(invalid_entries)
  } else {
    message("All 'substudy' entries are valid after corrections.")
  }
  
  return(substudy_column)
}
