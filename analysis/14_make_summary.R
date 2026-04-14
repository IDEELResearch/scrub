# 14_make_summary.R
#
# Author: Bob Verity
# Date: 2026-03-19
#
# Purpose:
# Makes a .xlsx version of the STAVE prevalence data for the supplementary materials.
#
# ------------------------------------------------------------------

fix_encoding <- function(df) {
  is_chr <- vapply(df, is.character, logical(1))
  df[is_chr] <- lapply(df[is_chr], function(x) {
    iconv(x, from = "", to = "UTF-8")
  })
  names(df) <- iconv(names(df), from = "", to = "UTF-8")
  df
}

insert_name <- function(x, name, pos) {
  x <- x[x != name]              # remove if already present
  append(x, name, after = pos - 1)
}

# ------------------------------------------------------------------

library(openxlsx)
library(readxl)
library(writexl)
library(stringi)

# read in combined STAVE object and existing summary xlsx
s <- readRDS(here("analysis", "data-out", "stave_data_2026.03.17.rds"))
summary_xl_path <- here("analysis", "data-out", "marker_prevalence_2026.03.17.xlsx")

s

# open connection to summary xlsx
wb <- openxlsx::loadWorkbook(summary_xl_path)

# replace studies sheet
removeWorksheet(wb, "studies")
addWorksheet(wb, "studies")
writeData(wb, "studies", fix_encoding(s$get_studies()))

# replace surveys sheet
removeWorksheet(wb, "surveys")
addWorksheet(wb, "surveys")
writeData(wb, "surveys", fix_encoding(s$get_surveys()))

# get data for all mutations
mut_names <- c("crt_76_T", "mdr1_86_Y", "k13_446_I", "k13_469_Y", "k13_476_I",
               "k13_493_H", "k13_539_T", "k13_553_L", "k13_561_H", "k13_574_L",
               "k13_580_Y", "k13_622_I", "k13_675_V", "k13_441_L", "k13_449_A",
               "k13_469_F", "k13_515_K", "k13_527_H", "k13_538_V", "k13_568_G")

for (i in seq_along(mut_names)) {
  message(sprintf("%s of %s", i, length(mut_names)))
  
  df_mut <- s$get_prevalence(target_variant = gsub("_", ":", mut_names[i])) |>
    select(study_id, survey_id, numerator, denominator, prevalence, prevalence_lower, prevalence_upper)
  
  removeWorksheet(wb, mut_names[i])
  addWorksheet(wb, mut_names[i])
  writeData(wb, mut_names[i], df_mut)
}

# save (overwrite file)
saveWorkbook(wb, summary_xl_path, overwrite = TRUE)

