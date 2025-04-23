# load stave  
# devtools::install_github("mrc-ide/STAVE")
library(STAVE)
library(tidyverse)

# read in data from databases (excluding geoff)
data <- readRDS("analysis/data-derived/final_data.rds")

# make an empty stave object
stave <- STAVE::STAVE_object$new()

# append data
data_stave <- convert_stave(data)

# STAVE errors related to Geoff that we can fix now: ------------------------

# Error 1 - Fixed Below: Not so much an error but a clean that should probably go into convert_stave
# Issue is, in the below where there the same variant string, once sorted, is a 
# duplicate in a particular survey. So in convert_stave (likely) we should first 
# order our variant strings, then group by study_key, surveykey and variant_string 
# and sum variant_num and leave total_num the same across s
# pf7k_1001_PF_ML_DJIMDE pf7k_1001_PF_ML_DJIMDE_MalariaGen_Bamako_2007 crt:76:K|T  5 47
# pf7k_1001_PF_ML_DJIMDE pf7k_1001_PF_ML_DJIMDE_MalariaGen_Bamako_2007 crt:76:T|K  7 47

# For now just fix here and group and sum
variant_order_pre <- data_stave$counts_dataframe |>
  mutate(combined_ID = paste(study_key, survey_key, sep = ":")) %>% 
  mutate(variant_string_order = variantstring::order_variant_string(variant_string))

variant_order <- variant_order_pre %>% 
  group_by(study_key, survey_key, combined_ID, variant_string_order) %>% 
  summarise(variant_string = unique(variant_string_order), 
            variant_num = sum(variant_num), 
            total_num = unique(total_num))

# should be TRUE if this has worked correctly
all(variant_order$variant_num<=variant_order$total_num)

# similarly this should be TRUE
((variant_order %>% 
  group_by(combined_ID, variant_string_order) %>% 
  mutate(tn = any(duplicated(variant_string_order))) %>% 
  filter(tn) %>% nrow()) == 0)

counts_df <- variant_order

# Error 2 - Fixed Below: The exact same variant cannot be present more than once in the same survey.
# This includes the same variant with the genes listed in different order. 
# For example, in the below STAVE does not understand how there can be two different
# denominators. 

# This is a MIP tool bug and we decided to take the higher as the denominator and
# proportion the variant num accordingly and make a note of it
# study_key      survey_key                    variant_string variant_num total_num
# S0013MsmtTza21 geoff_S0013MsmtTza21_Bkg_2021 mdr1:111:L               0       153
# S0013MsmtTza21 geoff_S0013MsmtTza21_Bkg_2021 mdr1:111:M               0       110

# Start by creating the position string
pre <- counts_df %>% mutate(position_string = variantstring::position_from_variant_string(variant_string)) %>% 
  group_by(combined_ID, position_string) %>% 
  mutate(dtn = length(unique(total_num))) 

# check that this is TRUE
pre %>% 
  filter(dtn == 1) %>%
  group_by(combined_ID, position_string) %>%
  mutate(variant_num = round(variant_num/total_num * max(total_num))) %>% 
  mutate(total_num = max(total_num)) %>% 
  pull(total_num) %>% all.equal(pre %>% filter(dtn == 1) %>% pull(total_num))

# update counts_df  
counts_df <- pre %>% 
  group_by(combined_ID, position_string) %>%
  mutate(variant_num = round(variant_num/total_num * max(total_num))) %>% 
  mutate(total_num = max(total_num)) 

# Error 3 - Fixed Below: There are some studies where the total number of variants do not add
# up to total_num. These Alfred has confirmed happen when the mip read stops half
# way through a codon so it is not an artefact. So similar trick as before except here
# we set the tootal num to be the sum of the variant_num

# in the counts_dataframe: for a given study, survey and variant, the sum of variant_num 
# cannot exceed the total_num. Problem rows in the counts_dataframe:
#   68281, 70913, 70929, 70932, 70921, 70937, 70945, 70953, 71745
# Error in private$check_structure(studies_dataframe, surveys_dataframe,  :

# S0026Verity geoff_S0026Verity_NordUbangiH_2013 dhps:436:A  1 2
# S0026Verity geoff_S0026Verity_NordUbangiH_2013 dhps:436:F  2 2

# this code below identifies where this happens
total_num_check <- dplyr::mutate(
  dplyr::summarise(
    dplyr::group_by(
      dplyr::mutate(
        counts_df %>% ungroup,
        row_number = dplyr::row_number()
      ),
      combined_ID,
      position_string
    ),
    total_num_check = sum(variant_num),
    total_num = total_num[1],
    row_number = row_number[1],
    .groups = "drop"
  ),
  overcount = total_num_check > total_num
)

# Below we filter to those locations and correct the total size
counts_df <- rbind(
  semi_join(counts_df, filter(total_num_check, overcount), by = c("combined_ID", "position_string")) %>% 
    group_by(combined_ID, position_string) %>% 
    mutate(total_num = sum(variant_num)),
  # and join with the antijoin, i.e. all the locations where there was not an overcount
  anti_join(counts_df, filter(total_num_check, overcount), by = c("combined_ID", "position_string"))
  ) %>% 
  ungroup() %>% 
  select(study_key, survey_key, variant_string, variant_num, total_num)

# Convert to STAVE and save-------------------

# Convert into a STAVE object
stave$append_data(studies_dataframe = data_stave$studies_dataframe,
                  surveys_dataframe = data_stave$surveys_dataframe,
                  counts_dataframe = counts_df) 

# Save the output in data-out ready for downstream analysis
dir.create("analysis/data-out", showWarnings = FALSE)
saveRDS(stave, "analysis/data-out/stave_data.rds")
