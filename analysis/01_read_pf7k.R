library(tidyverse)

# ---------------------------------------------------- o
# 1. Malariagen wrangle
# ---------------------------------------------------- o

mdf <- data.table::fread("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_drug_resistance_marker_genotypes.txt") %>%
  rename(sample = V1)
# writing in case they remove it later
write.csv(mdf, "analysis/data-raw/pf7k_raw.csv")
mdf <- data.table::fread("analysis/data-raw/pf7k_raw.csv") %>% as.data.frame()

# get the vcf grabbed GTs
mdrex <- read.csv(here::here("analysis/data-raw/mdr.csv")) %>%
  rename(sample = Sample)

# wrangle these to be similar to mdf styling
# For note, there is one sample with the second allele for mdr1_86
# but I can't find what this could be. Sample SPT35516. Setting to NA for now.
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
mdrex_df <- mdrex %>%
  mutate(`mdr1_86[Y]` = aa_annot(mdr1_N86Y, mdr1_N86Y_oth, "N", "Y")) %>%
  mutate(`mdr1_184[F]` = aa_annot(mdr1_Y184F, mdr1_Y184F_oth, "Y", "F")) %>%
  select(sample, `mdr1_86[Y]`, `mdr1_184[F]`)

# bind together
mdf <- left_join(mdf, mdrex_df, by = "sample")

# grab the meta information for getting lat long and sample year
meta <- read.csv("https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_samples.txt", sep = "\t")
# writing in case they remove it later
write.csv(meta, "analysis/data-raw/pf7k_samples.csv")
meta <- read.csv("analysis/data-raw/pf7k_samples.csv")

# and join
mdf <- left_join(mdf, meta %>% rename(sample = Sample), by = "sample")

# sort this into a consistent format that is the same as the WHO/WWWARN information

# 1. Rename annoying markers
mdf <- mdf %>% janitor::clean_names()

# select all that are needed
mdf <- mdf %>% select(
  sample:mdr1_184_f,
  study = study,
  admin_0 = country,
  admin_1 = admin_level_1,
  lat = admin_level_1_latitude,
  long = admin_level_1_longitude,
  year = year
) %>% 
  rename(k13_markers = kelch13_349_726_ns_changes)

# now start collapsing into helpful formats

# A. ART

# First clean the misings
mdf <- mdf %>%
  # 1. All blanks in k13 are WT
  mutate(k = replace(k13_markers, k13_markers == "", "WT")) %>%
  # 2. All "-" in k13 are missing (N = 843)
  mutate(k13_markers = replace(k13_markers, k13_markers == "-", NA)) %>%
  # 3. All "!" in k13 are frame-shift in the haplotype, consider missing (N = 4)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!", NA)) %>%
  # 4. All "!*" in k13 are frameshift and unphased het followed by het. consider missing (N = 1)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!*", NA)) %>%
  # 5. All "*" in k13 are unphased het followed by het. consider missing (N = 8)
  mutate(k13_markers = replace(k13_markers, k13_markers == "*", NA))

# FOR UNDERSTANDING
# Upper case = homozygous mutations
# lower case = heterozygous
# lower case without second mutation noted is heterozygous at SNP level but same AA
# , dictates separate and distinct haplotypes (different clones)
# / dictates additional NS changes in the same clone
# * in the context of mutations, these indicate that the sample could not be phased

convert_k13 <- function(marker) {
  
  # Case 1: NA - Missing
  if (is.na(marker)) {
    return(NA)
  }
  
  # Case 2: "" - Wild Type
  if (marker == "") {
    return("WT")
  }
  
  # Case 3:  No / or , therefore simple
  if (!grepl("\\/|\\,", marker)) {
    return(clean_mutations(paste0("k13-", marker)))
  }
  
  # Final Case
  spl <- toupper(strsplit(marker, "\\,|\\/")[[1]])
  vapply(paste0("k13-", spl), clean_mutations, character(1))
  
  
  spls <- strsplit(marker, ",", fixed = TRUE)
  
  lapply(spls, function(x){
    #message(x)
    if (length(x) == 2) {
      return(sum(grepl(validated, x, ignore.case = TRUE))/2)
    } else {
      if (is.na(x)){
        return(NA)
      } else if (x == "WT") {
        return(0)
      } else {
        ret <- as.integer(grepl(validated, x, ignore.case = TRUE))
        if(grepl("[[:lower:]]", x)) {
          ret <- ret/2
        }
        return(ret)
      }
    }
  }) %>% unlist
  
}

# send example to Bob to comment on
mdf %>% select(crt_K76T, mdr1_N86Y, mdr1_Y184F, mdr1_dup_call, k13_markers,pm2_dup_call) %>% 
  unique() %>% 
  write.csv("analysis/data-derived/pf7k_phase_examples.csv")

## NB FOR STAVE: "|" reflects phased 
## NB FOR STAVE: "/" reflects unphased 

mdf %>% 
  mutate(crt_K76T = case_when(
    crt_K76T == "T,K" ~ "crt:76:T|K",
    crt_K76T == "T" ~ "crt:76:T",
    crt_K76T == "K" ~ "crt:76:K",
    crt_K76T == "K,T" ~ "crt:76:K|T",
    crt_K76T == "K,T" ~ NA,
    # Just marking as phased here. But if we come back to do extended haplotypes then change this
    crt_K76T == "T,Q*" ~ "crt:76:T|Q", 
    crt_K76T == "K,P" ~ "crt:76:K|P"
  )) %>% 
  mutate(mdr1_N86Y = case_when(
    mdr1_86_y == "N" ~ "mdr1:86:N",
    mdr1_86_y == "Y" ~ "mdr1:86:Y",
    # The mdr1 SNPs I had to grab from read counts off vcfs, but can't tell if they are phased.
    # Have left as unphased while we just consider the single haplotype at the moment 
    mdr1_86_y == "N,Y" ~ "mdr1:86:N|Y"
  )) %>% 
  mutate(mdr1_Y184F = case_when(
    mdr1_184_f == "F" ~ "mdr1:184:F",
    mdr1_184_f == "Y" ~ "mdr1:184:Y",
    # The mdr1 SNPs I had to grab from read counts off vcfs, but can't tell if they are phased.
    # Have left as unphased while we just consider the single haplotype at the moment 
    mdr1_184_f == "Y,F" ~ "mdr1:184:Y|F"
  )) %>% 
  mutate(k13_markers = convert_k13(k13_markers))
  




# First clean the misings
mdf <- mdf %>%
  # 1. All "-" in crt are missing (N = 26)
  mutate(crt_K76T = replace(crt_K76T, crt_K76T == "-", NA)) %>%
  # 2. convert our crt markers into numeric for crt76T
  mutate(crt_76T = bires_conv(crt_K76T, "T", "K")) %>%
  # 3. convert our mdr markers into numeric for mdr1_86Y
  mutate(mdr1_86Y = bires_conv(mdr1_N86Y, "Y", "N")) %>%
  # 3. convert our mdr markers into numeric for mdr1_184F
  mutate(mdr1_184F = bires_conv(mdr1_Y184F, "F", "Y"))

# C. Duplications - convert the dup calls to resi status and NAs
mdf <- mdf %>%
  mutate(mdr1_CNV = mdr1_dup_call) %>%
  mutate(pfpm23_CNV = pm2_dup_call) %>%
  mutate(mdr1_CNV = replace(mdr1_CNV, mdr1_CNV == -1, NA)) %>%
  mutate(pfpm23_CNV = replace(pfpm23_CNV, pfpm23_CNV == -1, NA))

# FINALLY:
# and work out totals etc
pf7k_res_df <- mdf %>%
  select(-k13_marker_in, -k13_markers, -k13_valid, -crt_K76T, -mdr1_N86Y, -mdr1_Y184F, -mdr1_dup_call, -pm2_dup_call) %>% 
  pivot_longer(starts_with(c("k13","crt","mdr1","pfpm23"))) %>%
  mutate(gene = gsub("(.*)_(.*)", "\\1", name)) %>%
  rename(mut = name) %>%
  mutate(url = "https://pubmed.ncbi.nlm.nih.gov/36864926/",
         pmid = "36864926",
         database = "Pf7k",
         site = paste("Pf7k Study:", study),
         source = NA) %>%
  group_by(
    site, source, admin_0, admin_1, lat, long, year,
    gene, mut, pmid, url, database
  ) %>%
  summarise(n = sum(!is.na(value)),
            x = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  mutate(study_start_year = NA, study_end_year = NA) %>%
  mutate(prev = x/n) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>% 
  mutate(source = "MalariaGEN") %>% 
  mutate(gene_mut = vapply(mut, clean_mutations, character(1)))

saveRDS(pf7k_res_df, here::here("analysis/data-derived/pf7k_res_df.rds"))

