# 10_read_who.R
#
# Author: OJ Watson, Bob Verity
# Date: 2025-12-17
#
# Inputs:
# - data-raw/pf7k_raw.csv
# - data-raw/mdr.csv
# - data-raw/pf7k_samples.csv
# - data-raw/Pf7_study_details.csv
#
# Outputs:
# - data-derived/pf7_STAVE.rds
#
# Purpose:
# Reads in the Pf7 data. Wrangles into STAVE format for crt, mdr1, and k13. Saves to file.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)
library(readxl)
#remotes::install_github("mrc-ide/variantstring@v1.8.6")
library(variantstring)
#remotes::install_github("mrc-ide/STAVE@v2.0.2")
library(STAVE)
library(data.table)

# ---------------------------------------------------- o
# 1. Malariagen wrangle
# ---------------------------------------------------- o

# mdf <- data.table::fread("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_drug_resistance_marker_genotypes.txt") |>
#   rename(sample = V1)
# writing in case they remove it later
# write.csv(mdf, "analysis/data-raw/pf7k_raw.csv")
mdf <- data.table::fread(here("analysis", "data-raw", "pf7k_raw.csv")) |>
  as.data.frame()

# get the vcf grabbed GTs
mdrex <- read.csv(here::here("analysis/data-raw/mdr.csv")) |>
  rename(sample = Sample)

# wrangle these to be similar to mdf styling
# One sample (SPT35516) has a mixed "other" mutation. Setting to NA for now.
aa_annot <- function(ref, alt, ref_aa = "N", alt_aa = "Y") {
  hom_ref <- which(ref == 0 & alt == 0)
  het <- which(ref == 0 & alt == 1)
  hom_alt <- which(ref == 1 & alt == 1)
  aa <- rep(NA, length(ref))
  aa[hom_ref] <- ref_aa
  aa[het] <- paste(ref_aa, alt_aa, sep = ",")
  aa[hom_alt] <- alt_aa
  return(aa)
}
mdrex_df <- mdrex |>
  mutate(`mdr1_86[Y]` = aa_annot(mdr1_N86Y, mdr1_N86Y_oth, "N", "Y")) |>
  mutate(`mdr1_184[F]` = aa_annot(mdr1_Y184F, mdr1_Y184F_oth, "Y", "F")) |>
  select(sample, `mdr1_86[Y]`, `mdr1_184[F]`)

# bind together
mdf <- left_join(mdf, mdrex_df, by = "sample")

# grab the meta information for getting lat long and sample year
# meta <- read.csv("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_samples.txt", sep = "\t")

# writing in case they remove it later
# write.csv(meta, "analysis/data-raw/pf7k_samples.csv")
meta <- read.csv(here("analysis", "data-raw", "pf7k_samples.csv"))

# and join
mdf <- left_join(mdf, meta |> rename(sample = Sample), by = "sample")

# sort this into a consistent format that is the same as the WHO/WWWARN information

# 1. Rename annoying markers
mdf <- mdf |> janitor::clean_names()

# select all that are needed
mdf <- mdf |> select(
  sample:mdr1_184_f,
  study = study,
  admin_0 = country,
  admin_1 = admin_level_1,
  lat = admin_level_1_latitude,
  long = admin_level_1_longitude,
  year = year
) |> 
  rename(k13_markers = kelch13_349_726_ns_changes)

# convert study into valid identifier
mdf$study <- sprintf("pf7_%s", gsub("-", "_", mdf$study))

# count studies
length(unique(mdf$study))

# ----------------------------
# drop non-African samples

pf7_Africa <- c("Mauritania", "Gambia", "Guinea", "Kenya", "Tanzania", "Ghana", 
                "Burkina Faso", "Mali", "Malawi", "Uganda", "Democratic Republic of the Congo", 
                "Nigeria", "Madagascar", "Cameroon", "Côte d'Ivoire", "Ethiopia", 
                "Benin", "Senegal", "Gabon", "Sudan", "Mozambique")

mdf <- mdf |>
  filter(admin_0 %in% pf7_Africa)

# count studies
length(unique(mdf$study))

# ----------------------------
# crt

# drop bad values (actually replace these with NA to retain constant df size, as
# bad entries here might be fine for e.g. kelch)
df_crt <- mdf |>
  # drop missing
  mutate(crt_72_76_cvmnk = ifelse(crt_72_76_cvmnk == "-", NA, crt_72_76_cvmnk)) |>
  # replace frame-shifts with NA
  mutate(crt_72_76_cvmnk = ifelse(grepl(pattern = "!", x = crt_72_76_cvmnk), NA, crt_72_76_cvmnk)) |>
  # ignore unphased flags. These are phasing issues at the nucleotide level - the AA is still resolved
  mutate(crt_72_76_cvmnk = gsub("\\*", "", crt_72_76_cvmnk)) |>
  filter(!is.na(crt_72_76_cvmnk))

# get all unique crt haplotypes in data
crt_old <- unique(df_crt$crt_72_76_cvmnk)

# find the equivalent variantstring for each
crt_new <- rep(NA, length(crt_old))
for (i in seq_along(crt_new)) {
  x_split <- strsplit(crt_old[i], ",")[[1]]
  x_mat <- mapply(function (s) c(s, rep(NA, 5 - length(s))), strsplit(x_split, ""))
  x_consensus <- apply(x_mat, 1, function(y) {
    if (is.na(y[1]) & is.na(y[2])) {
      NA
    } else if (is.na(y[1])) {
      y[2]
    } else if (is.na(y[2])) {
      y[1]
    } else if (y[1] == y[2]) {
      y[1]
    } else {
      sprintf("%s|%s", y[1], y[2])
    }
  })
  if (!all(is.na(x_consensus))) {
    crt_new[i] <- paste0(x_consensus, collapse = "_")
  }
}

#cbind(crt_old, crt_new)

# replace crt
df_crt <- df_crt |>
  mutate(crt_72_76_cvmnk = crt_new[match(crt_72_76_cvmnk, crt_old)],
         crt_72_76_cvmnk = sprintf("crt:72-76:%s", crt_72_76_cvmnk),
         crt_72_76_cvmnk = variantstring::order_variant_string(crt_72_76_cvmnk))

# ----------------------------
# mdr

df_mdr1 <- mdf |>
  filter(!is.na(mdr1_86_y)) |>
  mutate(mdr1_86_y = gsub(",", "/", mdr1_86_y),
         mdr1_86_y = sprintf("mdr1:86:%s", mdr1_86_y))

# ----------------------------
# k13

# when we filter to African samples and WHO-validated loci, there are only a
# very small number of Pf7 k13 variants. These can be manually converted to
# variantstring format.
k13_stem <- "k13:441_446_449_458_469_476_481_493_515_527_537_538_539_543_553_561_568_574_580_622_675:"
k13_convert <- data.frame(k13_markers = c("g538v",
                                          "n458d",
                                          "N458D",
                                          "p441s"),
                          k13_variantstring = sprintf("%s%s", k13_stem, c("PFGNCMAYRPNVRIPRVPCRA",
                                                                          "PFGDCMAYRPNGRIPRVPCRA",
                                                                          "PFGDCMAYRPNGRIPRVPCRA",
                                                                          "SFGNCMAYRPNGRIPRVPCRA")))

df_k13 <- mdf |>
  left_join(k13_convert, by = join_by(k13_markers))

# all other samples should be marked as WT at the loci of interest
k13_WT <- sprintf("%s%s", k13_stem, "PFGNCMAYRPNGRIPRVPCRA")
df_k13 <- df_k13 |>
  mutate(k13_variantstring = ifelse(is.na(k13_variantstring), k13_WT, k13_variantstring))

# ----------------------------
# tidy up and make STAVE format

# sum up variants
crt_long <- df_crt |>
  group_by(study, lat, long, year, admin_0, admin_1, variant_string = crt_72_76_cvmnk) |>
  summarise(variant_num = n(),
            .groups = "drop") |>
  ungroup() |>
  group_by(study, lat, long, year, admin_0, admin_1) |>
  mutate(total_num = sum(variant_num)) |>
  ungroup()

mdr_long <- df_mdr1 |>
  group_by(study, lat, long, year, admin_0, admin_1, variant_string = mdr1_86_y) |>
  summarise(variant_num = n(),
            .groups = "drop") |>
  ungroup() |>
  group_by(study, lat, long, year, admin_0, admin_1) |>
  mutate(total_num = sum(variant_num)) |>
  ungroup()

k13_long <- df_k13 |>
  group_by(study, lat, long, year, admin_0, admin_1, variant_string = k13_variantstring) |>
  summarise(variant_num = n(),
            .groups = "drop") |>
  group_by(study, lat, long, year, admin_0, admin_1) |>
  mutate(total_num = sum(variant_num)) |>
  ungroup()

mdf_long <- bind_rows(crt_long, mdr_long, k13_long)

# make STAVE study-level data frame
df_study <- data.frame(study_id = unique(mdf_long$study),
                       study_label = NA,
                       description = NA,
                       access_level = "public",
                       contributors = NA,
                       reference = "https://pubmed.ncbi.nlm.nih.gov/36864926/",
                       reference_year = 2023,
                       PMID = 36864926)

# merge in details of individual Pf7 studies
pf7_details <- read.csv(here("analysis", "data-raw", "Pf7_study_details.csv"))

df_study <- df_study |>
  select(-c(study_label, contributors)) |>
  left_join(pf7_details) |>
  select(study_id, study_label, description, access_level, contributors, reference, reference_year, PMID)

# make STAVE survey-level data frame
df_survey <- mdf_long |>
  mutate(survey_id = sprintf("%s_%s_%s", study, admin_1, year),
         survey_id = gsub(" ", "_", survey_id),
         survey_id = gsub("-", "_", survey_id)) |>
  group_by(study, survey_id) |>
  summarise(country_name = admin_0[1],
            site_name = admin_1[1],
            latitude = lat[1],
            longitude = long[1],
            location_method = NA,
            location_notes = NA,
            collection_start = as.Date(sprintf("%s-01-01", year[1])),
            collection_end = as.Date(sprintf("%s-12-31", year[1])),
            collection_day = as.Date(sprintf("%s-07-01", year[1])),
            time_method = "only year of collection given. Start and end dates set at first and last days of that year, and collection day taken as midpoint",
            time_notes = NA,
            .groups = "drop") |>
  mutate(study_id = study) |>
  select(-study) |>
  relocate(study_id)

# make STAVE count-level data frame
df_counts <- mdf_long |>
  mutate(survey_id = sprintf("%s_%s_%s", study, admin_1, year),
         survey_id = gsub(" ", "_", survey_id),
         survey_id = gsub("-", "_", survey_id),
         notes = NA) |>
  select(study_id = study, survey_id, variant_string, variant_num, total_num, notes)

# make STAVE object and append data
s <- STAVE::STAVE_object$new()
s$append_data(studies_dataframe = df_study,
              surveys_dataframe = df_survey,
              counts_dataframe = df_counts)
s


# write to file
saveRDS(s, file = here("analysis", "data-derived", "pf7_STAVE.rds"))

