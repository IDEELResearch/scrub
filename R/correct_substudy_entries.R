#' Correct Known Bad Strings in the Substudy Column
#'
#' This function corrects known issues in the `substudy` column, such as converting 
#' "day 0" to "day0extracted", "untreated_extracted" to "untreatedextracted", 
#' and "untreadextracted" to "untreatedextracted". It counts and prints the number 
#' of corrections made for each known issue.
#'
#' @param substudy_column A character vector containing substudy entries that may 
#'        contain known typos or bad strings.
#' @return A character vector with corrected substudy entries. Prints the number of 
#'         corrections made for each issue.
#' @examples
#' substudy_data <- c("day 0", "untreated_extracted", "untreadextracted", "day0extracted")
#' corrected_substudy <- scrub:::correct_substudy_entries(substudy_data)
#' # Output:
#' # Corrected 1 instances of day0extracted
#' # Corrected 1 instances of untreated_extracted
#' # Corrected 1 instances of untreadextracted
#'
#' 
correct_substudy_entries <- function(substudy_column) {
  # Correct entries and count the corrections
  corrections <- list(
    day0extracted = sum(grepl("^day 0$", substudy_column)),  # Count "day 0"
    untreated_extracted = sum(grepl("^untreated_extracted$", substudy_column)),  # Count "untreated_extracted"
    untreadextracted = sum(grepl("^untreadextracted$", substudy_column))  # Count "untreadextracted"
  )
  
  # Perform the corrections
  substudy_column <- gsub("^day 0$", "day0extracted", substudy_column)
  substudy_column <- gsub("^untreated_extracted$", "untreatedextracted", substudy_column)
  substudy_column <- gsub("^untreadextracted$", "untreatedextracted", substudy_column)
  
  # Print the corrections made and their counts
  for (correction in names(corrections)) {
    if (corrections[[correction]] > 0) {
      message(paste("Corrected", corrections[[correction]], "instances of", correction))
    }
  }
  
  return(substudy_column)
}
