#' Impute the WT prevalence for non-mutant codons in WWARN dataset
#'
#' @details
#' 
#' This function takes the mut information and extracts the codon number. 
#' This is applied over all of the values in the column to determine the markers found
#' in that study
#' @name extract_from_entry
#' @param entry Character string. A string representing k13 SNP data in various formats.
#' This is programmed to handle double mutants, mixed infections etc.
#' @return Integer. Processed mutation information into the codon position.
#' @export
#' @examples
#' extract_from_entry("C580Y") = 580
#' extract_from_entry("C580Y_A481G") = 580   481
#' extract_from_entry("C580Y/G") = 580
#' mut <- c("C580Y", "A481G", "R622I")
#' unlist(lapply(mut, extract_from_entry)) = 580 481 622

# Function to extract codons from a single mutation entry
extract_from_entry <- function(entry) {
  # Handle double mutant (underscore)
  if (grepl("_", entry)) {
    parts <- unlist(strsplit(entry, "_"))
    codons <- as.integer(gsub(".*?([0-9]+).*", "\\1", parts))
  } else {
    # Mixed or single mutant: extract first number only, remove leading letters
    codon <- as.integer(gsub(".*?([0-9]+).*", "\\1", entry))
    codons <- codon
  }
  return(codons)
}

#' @details
#' 
#' This function adds additional rows to the dataframe. It keeps all column variables 
#' the same except for x and mut.
#' @name add_a_row
#' @export
#' 
#' @param df WWARN dataframe that you want to add a row to. 
#' Note: this will probably be a specific survey
#' @param x_new The value of x in the new row. For imputing wild type, x_new == n
#' @param mut_new The value of mut in the new row. For imputing wild type, these are the 
#' WT encodings for codons not detected in the survey

add_a_row <- function(df, x_new, mut_new) {
  rows <- nrow(df) 
  df <- df %>%
    bind_rows(slice(., 1))
  df$x[rows+1] <- x_new
  df$gene_mut[rows+1] <- mut_new
  df$prev <- df$x/df$n
  return(df)
}

#' @details
#' 
#' This imputes the WT in the survey.
#' @name impute_survey
#' @export
#' 
#' @param survey_df WWARN dataframe for this survey (note: not the study)
#' @param study_markers Validated markers, including the gene-mut of the WT,
#' within the min and max codons of that study
#' study_markers == study_validated in {impute_study}
#' @return df Returns the df of a survey with the imputed WT rows added
#' 
impute_survey <- function(survey_df, study_markers) {
  
  mutations <- survey_df$mut[survey_df$mut != "WT"]
  survey_mut <- unlist(lapply(mutations, extract_from_entry))
  survey_imputation <- study_markers |>
    dplyr::filter((codon %in% survey_mut) == FALSE) |>
    dplyr::pull(gene_mut) |>
    unique()
  
  if(length(survey_imputation) == 0) {
    if(length(mutations) != 0) {
      # impute the codons with reported mutants, from the dataframe with all of them
      non_val_wt <- mutation_key |>
        dplyr::filter(CODON %in% survey_mut) |>
        dplyr::mutate(gene_mut = paste0("k13:",CODON,":",REF))
      
      # sometimes there may only be a non validated marker
      # I just need to filter out 
      for(i in 1:length(non_val_wt)) {
        add_a_row(survey_df)
      }
      
    }
  }
  
  
  for(i in 1:length(survey_imputation)) {
    survey_df <- add_a_row(survey_df, x_new = survey_df$n[1], mut_new = survey_imputation[i])
  }
    # need to fix mut
  survey_df <- survey_df |>
    dplyr::filter(gene_mut != "k13:WT") |>
    dplyr::mutate(mut = gsub("k13:", "", gene_mut)) |> # does this need further refining?
    dplyr::mutate(mut = gsub(":", "", mut)) |> 
    dplyr::distinct() # remove the duplicates at codon w/ two mutants
  
  return(survey_df)
}

#' @details
#' 
#' This imputes the WT in the study. 
#' Comprises of several parts:
#' - Identifies the mutations in the study and extracts codon #
#' - Finds min and max codons
#' - Identifies validated markers in this range
#' - Imputes each survey in the study and rbinds them together
#' 
#' @name impute_study
#' @export
#' 
#' @param study_df WWARN dataframe for this survey 
#' @return df Returns the df of a study with the imputed WT rows added
impute_study <- function(study_df) {
  # Remove "WT" entries
  mutations <- study_df$mut[study_df$mut != "WT"] |>
    unique()
  
  # Apply to each entry and flatten the result
  codon_numbers <- unlist(lapply(mutations, extract_from_entry)) |>
    unique() |> sort()
  
  if(length(codon_numbers) == 0) {
    wt_only <- c(wt_only, unique(study_df$pmid))
    surveys <- NULL
  } else {
    # add these steps into the above function rather than separate
    min <- min(codon_numbers)
    max <- max(codon_numbers)
    
    # pull out all of the validated markers within the codon range of the study
    study_validated <- validated_markers |>
      dplyr::rowwise() |>
      dplyr::mutate(codon = extract_from_entry(mut)) |>
      dplyr::filter(codon >= min) |>
      dplyr::filter(codon <= max) |>
      dplyr::left_join(reference_validated, join_by("codon" == "CODON"))
    
    if(nrow(study_validated) == 0) {
      surveys <- study_df 
      wt_only <- c(wt_only, unique(study_df$pmid)) # add to this list
    } else {
      
      # now for each survey, filter out the rows with mutants
      surveys <- study_df |>
        dplyr::group_by(site, year) |>
        dplyr::mutate(siid = cur_group_id()) %>%
        split(.$siid)
      
      surveys <- lapply(surveys, function(df) impute_survey(df, study_validated)) %>% do.call(rbind,.)
      
    }
    
  }
  return(surveys)
}
