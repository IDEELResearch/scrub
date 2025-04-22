# load stave  
# devtools::install_github("mrc-ide/STAVE")
library(STAVE)
library(tidyverse)

# read in data from databases (excluding geoff)
data <- readRDS("analysis/data-derived/final_data.rds")

# hot fixes till upstream data entry and STAVE is fixed
data <- data %>% filter(!grepl("-|TK|BND", variant_string))

# make an empty stave object
stave <- STAVE::STAVE_object$new()

# append WWARN data
data_stave <- convert_stave(data)

# set URL for NA in WHO to WHO MTM
data_stave$studies_dataframe <- data_stave$studies_dataframe %>% 
  mutate(url = replace(url, is.na(url) & grepl("^who", study_id), "WHO MTM"))


# Erroring:

# Error 1: number of amino acid loci (5) must equal the number of codon positions (3)
# e.g. crt:72_73_74_75_76:CVMNK/CVIET/SVMNK which should be crt:72_73_74_75_76:C/S_V_M/I_N/E_K/T

# Error 2: entry 72945: amino acid sequence contains invalid characters. See ?allowed_amino_acids()
# e.g. k13:446:*

# For now grab these errors (but pass back to data entry most likely)
errors <- c(
  71776, 71777, 71778, 71779, 71780, 71781, 71782,
  71787, 71788, 71789, 71790, 71791, 71792, 71793,
  72699, 72701, 72703,
  72945, 72946, 72947, 72948, 72949, 72950, 72951,
  72952, 72953, 72956, 72957, 72958, 72959, 72960,
  72961, 72962, 72963, 72964,
  73102
)

# Error 3: NA variant num - 1,692 rows
# For now add these errors (but pass back to data entry most likely)
errors <- c(errors, which(is.na(data_stave$counts_dataframe$variant_num)))

# Error 4: Not so much an error but a clean that should probably go into convert_stave
# Issue is, in the below where there the same variant string, once sorted, is a 
# duplicate in a particular survey. So in convert_stave (likely) we should first 
# order our variant strings, then group by study_key, surveykey and variant_string 
# and sum variant_num and leave total_num the same across s
# pf7k_1001_PF_ML_DJIMDE pf7k_1001_PF_ML_DJIMDE_MalariaGen_Bamako_2007 crt:76:K|T  5 47
# pf7k_1001_PF_ML_DJIMDE pf7k_1001_PF_ML_DJIMDE_MalariaGen_Bamako_2007 crt:76:T|K  7 47

# For now just filter these out (n.b. we pass errors in here as these will fail the checks in order_variant_string)
variant_order <- data_stave$counts_dataframe[-errors,] |>
  mutate(combined_ID = paste(study_key, survey_key, sep = ":")) %>% 
  mutate(variant_string_order = variantstring::order_variant_string(variant_string))

counts_df <- variant_order %>% group_by(combined_ID, variant_string_order) %>% 
  mutate(tn = any(duplicated(variant_string_order))) %>% 
  filter(!tn) %>% 
  ungroup()

# Error 5: The exact same variant cannot be present more than once in the same survey.
# This includes the same variant with the genes listed in different order. 
# For example, in the below STAVE assumes that these should have been mcr1:111:L/M (presumably)
# But perhaps the correct fix is just to remove these - why do we care about recording zero?
# study_key      survey_key                    variant_string variant_num total_num
# S0013MsmtTza21 geoff_S0013MsmtTza21_Bkg_2021 mdr1:111:L               0       153
# S0013MsmtTza21 geoff_S0013MsmtTza21_Bkg_2021 mdr1:111:M               0       110

# For now just filter these out
counts_df <- counts_df %>% mutate(position_string = variantstring::position_from_variant_string(variant_string)) %>% 
  group_by(combined_ID, position_string) %>% 
  mutate(dtn = length(unique(total_num))) %>% 
  filter(dtn == 1) %>%
  ungroup() %>% 
  select(1:5)

# Error 6: There are some studies where the total number of variants do not add
# up to total_num. For example, in the below it likely should be that 
# 1 is A/F and another one is just F
# But probably should go back to data entry to fix these

# in the counts_dataframe: for a given study, survey and variant, the sum of variant_num 
# cannot exceed the total_num. Problem rows in the counts_dataframe:
#   68281, 70913, 70929, 70932, 70921, 70937, 70945, 70953, 71745
# Error in private$check_structure(studies_dataframe, surveys_dataframe,  :

# S0026Verity geoff_S0026Verity_NordUbangiH_2013 dhps:436:A  1 2
# S0026Verity geoff_S0026Verity_NordUbangiH_2013 dhps:436:F  2 2

# For now just filter these out (n.b. we pass errors in here as these will fail the checks in order_variant_string)
counts_df <- dplyr::anti_join(
  counts_df %>% mutate(gene = gsub("(^[[:alnum:]]*)(:.*)", "\\1", variant_string), .before = 1), 
  counts_df[c(68281, 70913, 70929, 70932, 70921, 70937, 70945, 70953, 71745), ] %>%
  mutate(gene = gsub("(^[[:alnum:]]*)(:.*)", "\\1", variant_string), .before = 1) %>% 
  dplyr::select(1:3)
) %>% select(-1)

# Convert into a STAVE object
stave$append_data(studies_dataframe = data_stave$studies_dataframe,
                  surveys_dataframe = data_stave$surveys_dataframe,
                  counts_dataframe = counts_df) 

# Save the output in data-out ready for downstream analysis
dir.create("analysis/data-out", showWarnings = FALSE)
saveRDS(stave, "analysis/data-out/stave_data.rds")
