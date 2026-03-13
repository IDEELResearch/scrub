# plot_sample_collection.R
#
# Author: Bob Verity
# Date: 2025-12-22
#
# Inputs:
# data-out/stave_data_ ... .rds
#
# Outputs:
# plots/sample_collection.pdf 
#
# Purpose:
# (this is an example header)
#
# ------------------------------------------------------------------
library(tidyverse)
library(here)
library(STAVE)

# read in combined STAVE object
s <- readRDS(here("analysis", "data-out", "stave_data_2026.03.13.rds"))

# get WHO positions
k13_dictionary <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET)) |>
  mutate(WT_variant = sprintf("k13:%s:%s", CODON, REF))

# get prevalence (of WT) over all positions
l <- list()
for (i in 1:nrow(k13_dictionary)) {
  message(sprintf("%s of %s", i, nrow(k13_dictionary)))
  l[[i]] <- s$get_prevalence(k13_dictionary$WT_variant[i]) |>
    select(study_id, survey_id, collection_day, denominator) |>
    mutate(target_variant = k13_dictionary$WT_variant[i],
           collection_year = year(collection_day))
}

# also get prevalence of PD
l[[22]] <- s$get_prevalence("crt:76:T") |>
  select(study_id, survey_id, collection_day, denominator) |>
  mutate(target_variant = "crt:76:T",
         collection_year = year(collection_day))
l[[23]] <- s$get_prevalence("mdr1:86:Y") |>
  select(study_id, survey_id, collection_day, denominator) |>
  mutate(target_variant = "crt:86:Y",
         collection_year = year(collection_day))

df_l <- bind_rows(l)

# find max denominator over all markers
df_comb <- df_l |>
  group_by(study_id, survey_id, collection_year) |>
  summarise(denom_max = max(denominator))

sum(df_comb$denom_max)

# ------------------------------------------------------------------------
# Plot samples collected over time

augmented_list <- c("s0002_ayelaw_eth2023",
                    "s0003_gidey_eth2022",
                    "s0006_WRAIR_kenread_ken",
                    "s0007_connelly_2024_zim",
                    "s0008_MSMT_TZA_22_23",
                    "s0010_karamoko_mgem_unpub",
                    "s0011_jacques-mari_gmms_unpub",
                    "s0021_roh_2023",
                    "s0063_ngasala_2024",
                    "s0074_lepiscopia_2025",
                    "s0075_eloff_2025")

mutate(bind_rows(l[1:21]), target_variant = "k13") |>
  bind_rows(mutate(l[[22]], target_variant = "crt 76")) |>
  bind_rows(mutate(l[[23]], target_variant = "mdr1 86")) |>
  group_by(study_id, survey_id, collection_year, target_variant) |>
  summarise(denom_max = max(denominator)) |>
  mutate(source = ifelse(grepl("^WWARN_", study_id), "WWARN",
                  ifelse(grepl("^WHO_", study_id), "WHO MTM",
                  ifelse(grepl("^MalariaGEN_", study_id), "MalariaGen Pf8",
                  ifelse(study_id %in% augmented_list, "Augmented datasets", "PRISMA review")))),
         source = factor(source, levels = c("PRISMA review", "Augmented datasets",
                                            "WWARN", "WHO MTM", "MalariaGen Pf8"))) |>
  group_by(collection_year, source, target_variant) |>
  summarise(n = sum(denom_max)) |>
  mutate(target_variant = factor(target_variant, levels = c("crt 76", "mdr1 86", "k13"))) |>
  ggplot() + theme_bw() +
  geom_col(aes(x = collection_year, y = n, fill = source)) +
  labs(x = "Collection Year", y = "Samples Sequenced",
       fill = "Data Source") +
  facet_wrap(~target_variant) +
  theme(legend.position = "bottom")

# save to file
ggsave(filename = here("analysis", "plots", "sample_collection.pdf"),
       width = 8, height = 4)
       

