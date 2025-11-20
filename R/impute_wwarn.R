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
#' extract_from_entry("C580Y") == 580
#' extract_from_entry("C580Y_A481G") == c(580, 481)
#' extract_from_entry("C580Y/G") == 580
#' mut <- c("C580Y", "A481G", "R622I")

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

#' Add a row to dataframe with all the same info except x and mut
#' 
#' @details
#' 
#' This function adds additional rows to the dataframe. It keeps all column variables 
#' the same except for x and mut.
#' @name add_a_row_k13
#' @export
#' 
#' @param df WWARN dataframe that you want to add a row to. 
#' Note: this will probably be a specific survey
#' @param x_new The value of x in the new row. For imputing wild type, x_new == n
#' @param gene_mut_new The value of mut in the new row. For imputing wild type, these are the 
#' WT encodings for codons not detected in the survey

add_a_row_k13 <- function(df, x_new, gene_mut_new) {
  rows <- nrow(df) 
  df <- df %>%
    dplyr::bind_rows(df[1, ])
  df$x[rows+1] <- x_new
  df$gene_mut[rows+1] <- gene_mut_new
  df$prev <- df$x/df$n
  df <- df |>
    dplyr::mutate(mut = gsub("k13:", "", gene_mut)) |> 
    dplyr::mutate(mut = gsub(":", "", mut)) 
  return(df)
}

#' Add a row to dataframe with all the same info except x and mut
#' 
#' @details
#' 
#' This function adds additional rows to the dataframe. 
#' It keeps all column variables 
#' the same except for x and mut.
#' @name add_a_row_pd
#' @export
#' 
#' @param df WWARN partner drug dataframe that you want to add a row to. 
#' Note: this will probably be a specific survey
#' @param x_new The value of x in the new row. For imputing wild type, x_new == n
#' @param mut_new The value of mut in the new row. For imputing wild type, these are the 
#' WT encodings for codons not detected in the survey

# write a function to add a row
add_a_row_pd <- function(df, x_new, mut_new) {
  rows <- nrow(df) 
  df <- df %>%
    bind_rows(slice(., 1))
  df$x[rows+1] <- x_new
  df$mut[rows+1] <- mut_new
  df$prev <- df$x/df$n
  if(sum(df$x) != df$n[1]) {
    print("x does not sum to n still")
  }
  return(df)
}


#' Impute the WT in surveys -- designed for cleaning WWARN
#' 
#' @details
#' 
#' This imputes the WT in the survey.
#' @name impute_survey
#' @export
#' 
#' @param survey_df WWARN dataframe for this survey (note: not the study)
#' @param impute_markers The df of markers to impute for the survey
#' When there are validated markers in the study, this is all of the WT markers
#' in the min-max codon range
#' In the case where only non-validated markers were found, this is the WT 
#' study_markers == study_validated in (impute_study)
#' @return df Returns the df of a survey with the imputed WT rows added
#' 
impute_survey <- function(survey_df, impute_markers) {
  ## now split by the classification
  classification <- impute_markers$classification[1]
  if((classification %in% c("validated", "non-validated"))==FALSE) {
    stop("Classification not eligible for imputation")
  }
  
  if(classification == "validated") {
    mutations <- survey_df$mut[survey_df$mut != "WT"]
    survey_mut <- unlist(lapply(mutations, extract_from_entry))
    survey_imputation <- impute_markers |>
      dplyr::filter((codon %in% survey_mut) == FALSE) |>
      dplyr::pull(gene_mut) |>
      unique()
    if(length(survey_imputation) > 0) {
      for(i in 1:length(survey_imputation)) {
        survey_df <- add_a_row_k13(survey_df, x_new = survey_df$n[1], gene_mut_new = survey_imputation[i])
      }
    }
  } else {
    # this is where impute_markers is only WT at the codons with mutants
    # WT and one mutant in survey -- remove WT
    mutations <- survey_df$mut[survey_df$mut != "WT"]
    survey_mut <- unlist(lapply(mutations, extract_from_entry))
    survey_imputation <- impute_markers |>
      dplyr::filter((codon %in% survey_mut) == FALSE) |>
      dplyr::pull(gene_mut) |>
      unique()
    if(length(survey_imputation) > 0) {
      for(i in 1:length(survey_imputation)) {
        survey_df <- add_a_row_k13(survey_df, x_new = survey_df$n[1], gene_mut_new = survey_imputation[i])
      }
    }
  }
  # now fix the formatting
  survey_df <- survey_df |>
    dplyr::filter(gene_mut != "k13:WT") |>
    dplyr::mutate(mut = gsub("k13:", "", gene_mut)) |> 
    dplyr::mutate(mut = gsub(":", "", mut)) 
  return(survey_df)
}


#' Impute the WT in studies by applying impute_survey -- designed for cleaning WWARN
#' 
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
      dplyr::left_join(reference_validated, join_by("codon" == "CODON")) |>
      dplyr::mutate(classification = "validated") |>
      dplyr::select(gene, mut, gene_mut, codon, classification) |>
      dplyr::arrange(codon) |>
      dplyr::distinct(across(-mut), .keep_all = TRUE)  |>
      dplyr::select(gene, mut, gene_mut, codon, classification)
    study_impute <- study_validated
    
    # if there are only non-validated codons with mutants reports
    if(nrow(study_validated) == 0) {
      # need to create a df with at least codon and gene_mut for the study imputation WT
      wt_impute <- mutation_key |>
        dplyr::filter(CODON %in% codon_numbers) |> 
        dplyr::mutate(gene_mut = paste0("k13:",CODON,":",REF)) |>
        dplyr::mutate(classification = "non-validated") |>
        dplyr::rename(gene = PROTEIN,
                      codon = CODON) |>
        dplyr::mutate(mut = gsub("k13:", "", gene_mut)) |> 
        dplyr::mutate(mut = gsub(":", "", mut)) |>
        dplyr::select(gene, mut, gene_mut, codon, classification) |>
        dplyr::arrange(codon)
      
      study_impute <- wt_impute
    }
    
    # now for each survey, filter out the rows with mutants
    surveys <- study_df |>
      dplyr::group_by(site, year) |>
      dplyr::mutate(siid = cur_group_id()) %>%
      split(.$siid)
    
    # # to allow debugging
    # study_impute
    # for(j in 1:length(surveys)) {
    #   impute_survey(surveys[[j]], study_impute)
    # }
    
    surveys <- lapply(surveys, function(df) impute_survey(df, study_impute)) %>% do.call(rbind,.)
    
  }
  return(surveys)
}

