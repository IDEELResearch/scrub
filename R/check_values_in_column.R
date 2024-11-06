#' Check whether columns in the dataframe satisfy allowed values
#'
#' This function looks at all values in a column and returns a warning if it is not within the set of allowed values
#' This function is generalisable to any column in the dataset and a defined suitable list
#'
#' @param data The preliminarily cleaned data where the format has been standardised, "nas" fixed etc.
#' @param column_name A string in " " with the name of the column of interest
#' @param allowed_values A vector of the values allowed within column_name
#' 
#' @return either a warning or a message detailing whether the test passes or fails 
#' 
#' @export
#' 
#' @examples
#' 
#' data <- data.frame(country = c("Angola", "Burundi", "Eritrea", "Vietnam"))
#' countries <- data.frame(country = countrycode::codelist$country.name.en,
#'   continent = countrycode::codelist$continent)
#' african_countries <- countries %>%
#'   dplyr::filter(continent == "Africa") %>%
#'   dplyr::pull(country) %>% unique()
#' check_values_in_column(data = data, column_name = "country", allowed_values = african_countries)
#' pub <- data.frame(publication_status = c("peer_reviewed", "pre-print", "na"))
#' pub_status <- c("peer_reviewed", "preprint", NA)
#' check_values_in_column(data = pub, column_name = "publication_status", allowed_values = pub_status)

check_values_in_column <- function(data, column_name, allowed_values) {
  # Check if any values in the column are not in the allowed_values vector
  if (any(!data[[column_name]] %in% allowed_values)) {
    index <- which((data[[column_name]] %in% allowed_values) == FALSE)
    vals <- data[[column_name]][index] %>% unique() %>% c()
    warning(paste("Warning: Column", column_name, "contains values not in the allowed list"))
    print(vals)
  } else {
    message(paste("All values in column", column_name, "are within the allowed list"))
  }
}
