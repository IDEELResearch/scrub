library(tidyverse)
devtools::load_all()

# ---------------------------------------------------- o
# 1. Malariagen wrangle
# ---------------------------------------------------- o

# mdf <- data.table::fread("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_drug_resistance_marker_genotypes.txt") |>
#   rename(sample = V1)
# writing in case they remove it later
# write.csv(mdf, "analysis/data-raw/pf7k_raw.csv")
mdf <- data.table::fread("analysis/data-raw/pf7k_raw.csv") |> as.data.frame()

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
meta <- read.csv("analysis/data-raw/pf7k_samples.csv")

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

# now start collapsing into helpful formats

# ----------------------------
# crt

# drop bad values (actually replace these with NA to retain constant df size, as
# bad entries here might be fine for e.g. kelch)
mdf <- mdf |>
  # drop missing
  mutate(crt_72_76_cvmnk = ifelse(crt_72_76_cvmnk == "-", NA, crt_72_76_cvmnk)) |>
  # replace frame-shifts with NA
  mutate(crt_72_76_cvmnk = ifelse(grepl(pattern = "!", x = crt_72_76_cvmnk), NA, crt_72_76_cvmnk)) |>
  # ignore unphased flags. These are phasing issues at the nucleotide level - the AA is still resolved
  mutate(crt_72_76_cvmnk = gsub("\\*", "", crt_72_76_cvmnk))

# get all unique crt haplotypes in data
crt_old <- unique(mdf$crt_72_76_cvmnk)

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
mdf <- mdf |>
  mutate(crt_72_76_cvmnk = crt_new[match(crt_72_76_cvmnk, crt_old)],
         crt_72_76_cvmnk = sprintf("crt:72-76:%s", crt_72_76_cvmnk))

# ----------------------------
# mdr

mdf <- mdf |>
  mutate(mdr1_86_y = gsub(",", "|", mdr1_86_y),
         mdr1_184_f = gsub(",", "|", mdr1_184_f),
         mdr1_86_184 = case_when(
           is.na(mdr1_86_y) & is.na(mdr1_184_f) ~ NA,
           is.na(mdr1_86_y) ~ sprintf("mdr1:184:%s", mdr1_184_f),
           is.na(mdr1_184_f) ~ sprintf("mdr1:86:%s", mdr1_86_y),
           TRUE ~ sprintf("mdr1:86_184:%s_%s", mdr1_86_y, mdr1_184_f)
         ))

# ----------------------------
# k13

table(mdf$k13_markers)

WHO_pos <- c(441, 446, 449, 458, 469, 469, 476, 481, 493, 515, 527, 537, 
             537, 538, 539, 543, 553, 561, 568, 574, 580, 622, 675)

# ----------------------------


# A. ART

# First clean the missing values
mdf <- mdf |>
  # 1. All blanks in k13 are WT
  mutate(k13_markers = replace(k13_markers, k13_markers == "", "WT")) |>
  # 2. All "-" in k13 are missing (N = 843)
  mutate(k13_markers = replace(k13_markers, k13_markers == "-", NA)) |>
  # 3. All "!" in k13 are frame-shift in the haplotype, consider missing (N = 4)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!", NA)) |>
  # 4. All "!*" in k13 are frameshift and unphased het followed by het. consider missing (N = 1)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!*", NA)) |>
  # 5. All "*" in k13 are unphased het followed by het. consider missing (N = 8)
  mutate(k13_markers = replace(k13_markers, k13_markers == "*", NA))

# function to split duplicated entries
split_rows_by_and <- function(data, column) {
  # Use tidyr's separate_rows to split the specified column
  data |>
    separate_rows(!!sym(column), sep = "\\s*&&\\s*")
}

# clean the biallelic and k13 markers
mdf2 <- mdf |> 
  mutate(across(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$"), pf7k_format_bia_for_stave)) |> 
  rowwise() |> 
  mutate(k13_markers = pf7k_format_k13_for_stave(k13_markers))

# now clean the variant_string entries with && in and bind with the rest
res <- rbind(
  mdf2 |> pivot_longer(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$|k13_markers")) |> 
    select(sample, study, admin_0, admin_1, lat, long, year, name, value) |> 
    filter(grepl("&&", value)) |> 
    split_rows_by_and("value"),
  mdf2 |> pivot_longer(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$|k13_markers")) |> 
    select(sample, study, admin_0, admin_1, lat, long, year, name, value) |> 
    filter(!grepl("&&", value))
)

# FINALLY:
# and work out totals etc
pf7k_res_df <- res |> 
  na.omit() |> 
  group_by(across(study:name)) |>
  mutate(total_num = length(unique(sample))) |>
  ungroup() |>
  group_by(across(study:total_num)) |>
  summarise(variant_string = unique(value), 
            variant_num = sum(!is.na(value))) |>
  rename(study_ID = study) |>
  mutate(url = "https://pubmed.ncbi.nlm.nih.gov/36864926/",
         pmid = "36864926",
         study_name = study_ID,
         study_type = "peer_reviewed",
         database = "Pf7k",
         site = paste("Pf7k Study:", study_ID),
         source = NA,
         authors = "MalariaGen",
         publication_year = 2023,
         survey_ID = paste0(study_ID, "_", authors, "_", admin_1, "_", year),
  ) |> 
  ungroup()

pf7k_res_df |>
  head()

# and lastly sort out the WT calls
mutation_key_path <- here::here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")
mutation_key <- read.csv(mutation_key_path)
indices_to_transform <- which(pf7k_res_df$variant_string == "k13:WT")

# range <- "349-726" - this is the range noted in the original pf7k file
pf7k_res_df$variant_string[indices_to_transform] <- gsub("K13", "k13", collapse_k13_range("k13:349:726", mutation_key))

saveRDS(pf7k_res_df, here::here("analysis/data-derived/pf7k_res.rds"))

