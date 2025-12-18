# 01_make_PMID_wildtype.R
#
# Author: Bob Verity
# Date: 2025-12-09
#
# Inputs:
# - data-raw/PMID_to_k13_coverage.csv
# - data-raw/k13_ref_protein_codon_dictionary.csv
#
# Outputs:
# - data-derived/PMID_to_k13_coverage.csv
#
# Purpose:
# Reads in PMID_to_k13_coverage.csv, which flags for each study (PMID) which WHO
# k13 loci were sequenced. This is used to generate the corresponding wild-type
# variant for each study, i.e. the expected AA code if no non-synonymous
# mutations were observed at any of these sequenced loci. This information is
# saved in a modified version of the original data frame in the data-derived
# folder.
#
# ------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(here)

# read in coverage data
k13_coverage <- read.csv(here("analysis", "data-raw", "PMID_to_k13_coverage.csv"))

# get WT allele for all WHO candidate and validated positions
WHO_targets <- read.csv(here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")) |>
  filter(!is.na(WHO_TARGET))

# work out WT in variantstring format from coverage info
k13_coverage$WT_variant <- NA
cov_mat <- k13_coverage |>
  select(starts_with("X")) |>
  as.matrix()
for (i in 1:nrow(k13_coverage)) {
  w <- which(cov_mat[i,] == 1)
  if (any(w)) {
    k13_coverage$WT_variant[i] <- sprintf("k13:%s:%s",
                                          paste(WHO_targets$CODON[w], collapse = "_"),
                                          paste(WHO_targets$REF[w], collapse = ""))
  }
}

# save back to file
write.csv(k13_coverage, file = here("analysis", "data-derived", "PMID_to_k13_coverage.csv"), row.names = FALSE)
