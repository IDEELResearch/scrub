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
  mutate(k13_markers = replace(k13_markers, k13_markers == "", "WT")) %>%
  # 2. All "-" in k13 are missing (N = 843)
  mutate(k13_markers = replace(k13_markers, k13_markers == "-", NA)) %>%
  # 3. All "!" in k13 are frame-shift in the haplotype, consider missing (N = 4)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!", NA)) %>%
  # 4. All "!*" in k13 are frameshift and unphased het followed by het. consider missing (N = 1)
  mutate(k13_markers = replace(k13_markers, k13_markers == "!*", NA)) %>%
  # 5. All "*" in k13 are unphased het followed by het. consider missing (N = 8)
  mutate(k13_markers = replace(k13_markers, k13_markers == "*", NA))

# send example to Bob to comment on
# mdf %>% select(crt_76_k, mdr1_86_y, mdr1_184f, k13_markers) %>% 
#   unique() %>% 
#   write.csv("analysis/data-derived/pf7k_phase_examples.csv")


# FOR UNDERSTANDING
# Upper case = homozygous mutations
# lower case = heterozygous
# lower case without second mutation noted is heterozygous at SNP level but same AA
# , dictates separate and distinct haplotypes (different clones)
# / dictates additional NS changes in the same clone
# * in the context of mutations, these indicate that the sample could not be phased

# function to convert various pf7k k13 SNP formats for stave
pf7k_format_k13_for_stave <- function(marker) {
  
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
    return(clean_mutations(paste0("k13-", gsub("*", "", toupper(marker), fixed = TRUE))))
  }
  
  # Case 4:  Only ,
  if (grepl("\\,", marker) && !grepl("\\/", marker)) {
    
    # split these out and clean
    marker <- gsub("*", "", marker, fixed = TRUE)
    spl <- toupper(strsplit(marker, "\\,|\\/")[[1]])
    res <- vapply(paste0("k13-", spl), clean_mutations, character(1))
    
    # these will become two rows of data as these are different clones
    # so collapse like this to be split apart later
    res <- paste0(as.character(res), collapse = "&&")
    return(res)
  }
  
  # Case 5:  Only /
  if (grepl("\\/", marker) && !grepl("\\,", marker)) {
    
    # First work out all the loci and snps for this sample
    marker <- gsub("*", "", marker, fixed = TRUE)
    spl <- toupper(strsplit(marker, "\\,|\\/")[[1]])
    all_loci <- unique(gsub("(^\\w)(\\d*)(\\w)", "\\2", spl))
    all_snps <- unique(gsub("(^\\w)(\\d*)(\\w)", "\\3", spl))
    all_wts <- unique(gsub("(^\\w)(\\d*)(\\w)", "\\1", spl))
    
    # And concatenate and order it
    res <- paste0("k13:",
                           paste0(sort(as.numeric(all_loci)), collapse = "_"),
                           paste0(":", paste0(all_snps[order(as.numeric(all_loci))], collapse = "_"))
                    )
    
    # these will become two rows of data as these are different clones
    # so collapse like this to be split apart later
    res <- paste0(as.character(unlist(res)), collapse = "&&")
    return(res)
  }
  
  # Case 6 , and /
  if (grepl("\\,", marker) && grepl("\\/", marker)) {
  
  # First work out all the loci and snps for this sample
  marker <- gsub("*", "", marker, fixed = TRUE)
  spl <- toupper(strsplit(marker, "\\,")[[1]])
  spl2 <- strsplit(spl, "\\/")
  all_loci <- lapply(spl2, function(x) {unique(gsub("(^\\w)(\\d*)(\\w)", "\\2", x))})
  all_snps <- lapply(spl2, function(x) {unique(gsub("(^\\w)(\\d*)(\\w)", "\\3", x))})
  all_wts <- lapply(spl2, function(x) {unique(gsub("(^\\w)(\\d*)(\\w)", "\\1", x))})
  
  # And concatenate and order it
  res <- lapply(seq_along(all_loci), 
                function(x) {
                  paste0("k13:",
                paste0(sort(as.numeric(all_loci[[x]])), collapse = "_"),
                paste0(":", paste0(all_snps[[x]][order(as.numeric(all_loci[[x]]))], collapse = "_"))
                  )
                })
  
  # these will become two rows of data as these are different clones
  # so collapse like this to be split apart later
  res <- paste0(as.character(unlist(res)), collapse = "&&")
  return(res)
  
  }
  

}

# function to convert biallelic SNP formats for stave
pf7k_format_bia_for_stave <- function(mut, gls = deparse(substitute(mut))) {
  
  # what is the gene loci snp
  gl <- gsub("(^[[:alnum:]]*_\\d*)(_[[:alnum:]]$)", "\\1", gls)
  gl <- gsub("_", ":", gl)
  
  # convert to stave phased format
  res <- gsub("\\,", "\\|", mut)
  
  # if any unphased (*) handle this
  unphased <- grep("\\*", res)
  if (length(unphased) > 0) {
    res[unphased] <- gsub("\\|", "\\/", res[unphased])
    res[unphased] <- gsub("\\*", "", res[unphased])
  }
  
  # remove any missing (-)
  missing <- grep("\\-", res)
  if (length(missing) > 0) {
    res[missing] <- NA
  }  
  
  # remove if same SNP (G|G)
  rept <- strsplit(res, "|")
  if (any(lengths(rept) == 3)) {
    pos <- which(lengths(rept) == 3)
    fix <- unlist(lapply(rept[pos], function(x){x[1] == x[3]}))
    if (any(fix)){
    correction <- unlist(lapply(rept[pos[which(fix)]], function(x){x[1]}))  
    res[pos[which(fix)]] <- correction
    }
  }  
  
  # and paste together with gene 
  res[!is.na(res)] <- paste0(gl, ":", res[!is.na(res)])
  
  return(res)
}

# function to split duplicated entries
split_rows_by_and <- function(data, column) {
  # Use tidyr's separate_rows to split the specified column
  data %>%
    separate_rows(!!sym(column), sep = "\\s*&&\\s*")
}
# function to duplicate samples for 

## NB FOR STAVE: "|" reflects phased 
## NB FOR STAVE: "/" reflects unphased 

# make this k13 marker and biallelic data
mdf2 <- mdf %>% 
  mutate(across(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$"), pf7k_format_bia_for_stave)) %>% 
  rowwise() %>% 
  mutate(k13_markers = pf7k_format_k13_for_stave(k13_markers))
  
# now do the clean on the double ands
res <- rbind(
mdf2 %>% pivot_longer(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$|k13_markers")) %>% 
  select(sample, study, admin_0, admin_1, lat, long, year, name, value) %>% filter(grepl("&&", value)) %>% 
  split_rows_by_and("value"),
mdf2 %>% pivot_longer(matches("^[[:alnum:]]*_\\d*_[[:alnum:]]$|k13_markers")) %>% 
  select(sample, study, admin_0, admin_1, lat, long, year, name, value) %>% filter(!grepl("&&", value))
)

# FINALLY:
# and work out totals etc
pf7k_res_df <- res %>% 
  na.omit() %>% 
  group_by(across(study:name)) %>%
  mutate(total_num = sum(!is.na(value))) %>%
  ungroup() %>%
  group_by(across(study:total_num)) %>%
  summarise(variant_string = unique(value), 
            variant_num = sum(!is.na(value))) %>%
  rename(study_ID = study) %>%
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
  ) %>% 
  ungroup()

# and lastly sort out the WT calls
mutation_key_path <- here("analysis", "data-raw", "k13_ref_protein_codon_dictionary.csv")
mutation_key <- read.csv(mutation_key_path)
indices_to_transform <- which(pf7k_res_df$variant_string == "k13:WT")
# range <- "349-726" - this is the range noted in the original pf7k file
pf7k_res_df$variant_string[indices_to_transform] <- gsub("K13", "k13", collapse_k13_range("k13:349:726"))

saveRDS(pf7k_res_df, here::here("analysis/data-derived/pf7k_res.rds"))

