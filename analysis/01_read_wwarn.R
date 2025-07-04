library(tidyverse)
devtools::load_all()

# A. ART
# ---------------------------------------------------- o
# 1. WWARN database compile ----
# ---------------------------------------------------- o

# as of 2 Jan 2024 ART-R markers (valid and candidate)
## when we have new candidates of validated markers, simply update the csv
validated <- read_csv("analysis/data-raw/mutation_dictionary.csv")
validated <- validated %>%
  dplyr::filter(mut != "CNV") %>% 
  dplyr::filter(gene == "k13") %>% 
  dplyr::mutate(gene_mut = paste0(gene, "-", substring(mut,2)))
validated_mut <- paste0(validated$mut, collapse = "|")

# read in WWARN study info
k13ww <- readxl::read_xls("analysis/data-raw/WWARN_K13_database_04-12-2033.xls", sheet = 1)
pdww <- readxl::read_xls("analysis/data-raw/WWARN_partnerdrug_database_04-12-2023.xls", sheet = 1)

# add iso
k13ww$iso <- countrycode::countrycode(k13ww$country, "country.name.en", "iso3c")
pdww$iso <- countrycode::countrycode(pdww$country, "country.name.en", "iso3c")

## 1.1. WWARN K13 admin mapping ----

# We need to find admin 1 globally based on lat long.
# file sourced from Malaria Atlas Project
# map_1 <- malariaAtlas::getShp("ALL",admin_level = c("admin1")) %>% sf::st_as_sf()
# saveRDS(map_1, file = "analysis/data-raw/alt_map.rds")
map_1 <- readRDS("analysis/data-raw/alt_map.rds")
k13_map <- map_1[map_1$iso %in% unique(k13ww$iso), ]

# create coords
k13coords <- sf::st_as_sf(k13ww %>% select(lat, lon), coords = c("lon", "lat"), crs = sf::st_crs(goodmap))

# identify k13 matches
k13ins <- as.integer(sf::st_within(k13coords, k13_map, prepared = TRUE))

# There are 2 samples in French Guiana that the approach below for nearest admin makes no sense for so correct after
# Having gone to manuscript these represent travellers from FG that were sequenced in France
ggplot() +
  geom_sf(data = k13_map[,"name_0"], fill = "gray90", colour = "white") +
  geom_point(data = k13ww[which(is.na(k13ins)), ], size = 2, aes(color = country, x = as.numeric(lon), y = as.numeric(lat))) +
  ggrepel::geom_text_repel(data = k13ww[which(is.na(k13ins)), ], size = 2, aes(label = country, x = as.numeric(lon), y = as.numeric(lat)), max.overlaps = 30) +
  theme_minimal() +
  labs(title = "Distribution of samples with missing k13ins")

k13nears <- as.integer(sf::st_nearest_feature(k13coords, k13_map))
k13ins[which(is.na(k13ins))] <- k13nears[which(is.na(k13ins))]
k13ww$admin_1 <- k13_map$name_1[k13ins]

# These are the 2 samples that this approach makes no sense for
k13ww$admin_1[which(k13ww$site == "Saint-Georges")] <- "Cayenne"

## 1.2. WWARN PD admin mapping ----

# start with the map and the coords
pd_map <- map_1[map_1$iso %in% unique(pdww$iso), ]
pdcoords <- sf::st_as_sf(pdww %>% select(lat, lon), coords = c("lon", "lat"), crs = sf::st_crs(goodmap))

# identify pd matches
pdins <- as.integer(sf::st_within(pdcoords, pd_map, prepared = TRUE))

# There is a clump of samples in India that the approach below for nearest admin makes no sense for so correct after
ggplot() +
  geom_sf(data = pd_map[,"name_0"], fill = "gray90", colour = "white") +
  geom_point(data = pdww[which(is.na(pdins)), ], size = 2, aes(color = country, x = as.numeric(lon), y = as.numeric(lat))) +
  ggrepel::geom_text_repel(data = pdww[which(is.na(pdins)), ], size = 2, aes(label = country, x = as.numeric(lon), y = as.numeric(lat)), max.overlaps = 30) +
  theme_minimal() +
  labs(title = "Distribution of samples with missing pdins")

pdnears <- as.integer(sf::st_nearest_feature(pdcoords, pd_map))
pdins[which(is.na(pdins))] <- pdnears[which(is.na(pdins))]
pdww$admin_1 <- pd_map$name_1[pdins]

pdww$admin_1[which(pdww$site == "West Bengal")] <- "West Bengal"

# Lastly save to make downstream analysis quicker as the above is a little slow...
saveRDS(k13ww, "analysis/data-derived/k13ww.rds")
saveRDS(pdww, "analysis/data-derived/pdww.rds")

# ---------------------------------------------------- o
# 2. Sort k13 ----
# ---------------------------------------------------- o

# GINA: Just start here as the above is long to run
k13ww <- readRDS("analysis/data-derived/k13ww.rds")

# sort names as wanted
k13wwdf <- k13ww %>%
  rename(admin_0 = country) %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  rename(long = lon) %>%
  mutate(study_start_year = NA) %>% # start and end year do not exist in wwarn k13 - only "year"
  mutate(study_end_year = NA) %>%
  mutate(n = as.integer(tested)) %>%
  mutate(x = as.integer(present)) %>%
  mutate(prev = x/n) %>%
  mutate(gene = "k13") %>%
  rename(mut = mutation) %>%
  mutate(pmid = pubMedId) %>%
  mutate(url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid)) %>%
  mutate(url = replace(url, grepl("NA",url), NA)) %>%
  mutate(database = "WWARN") %>%
  rename(source = authors) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

# sort out val grouping issue
k13wwdf <- bind_rows(
  k13wwdf %>% group_by(across(c(-x, -n, -prev))) %>% 
    mutate(ns = n()) %>% 
    filter(ns==1), # because we removed val, we have extra entries , 
  k13wwdf %>% group_by(across(c(-x, -n, -prev))) %>% 
    mutate(ns = n()) %>% 
    filter(ns>1) %>% # because we removed val, we have extra entries 
    summarise(x = sum(x), n = sum(x), prev = x/n)
) %>% 
  select(-ns) %>% 
  ungroup()

# ---------------------------------------------------- o
## 2.1 Correct data extraction from source for known pubs  ----
# ---------------------------------------------------- o

# before doing manual cleaning - one study corresponds to many of the issues
# fails for Gina so commenting out for now
# TODO: make sure this is in the PR
# tf <- tempfile()
# download.file("https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/jid/223/6/10.1093_infdis_jiaa687/1/jiaa687_suppl_supplementary_table_s3.xlsx?Expires=1748452824&Signature=PjWn-A7qw1OTfQo4hrvOee1909kD~hyO-lT9jA2dQIhZmS8Jmcp8-HjqydW9iI3GuNx3sLPFz2GmT-Hm7i-58NoOSff~g84yASXtA3HcAJOy6maFSWXMNu3PtZibzklibTLZdXFqzh1n0aQCFsrzs8YEN70NE6hQZE9skxgQnrG6ipfqQWGICNkkzdCea-AcoP4Clgcma159BE3Mlac4R1RGtuzbTCYgEk-End4tck-Zk0ApCnjhnLcM7zA4DSLEyJ98rka29DEwKG8raNDj~MXR5yh1i6~6MScgAQpM6JIVB7UM~~euH4Qsc3XJlrR63kuVUUwzk12aLXrcuO6S7A__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA", tf)
# 
# # clean up and create our correct table
# df <- readxl::read_excel(tf) %>% janitor::clean_names()
df <- df %>% filter(gene_name == "k13") %>% filter(!(grepl("fs", mutation_name))) # remove frame shift
df_mut <- df %>% filter(gene_name=="k13") %>% filter(!is.na(genotype)) %>% pull(mutation_name) %>% unique()
aa_lookup <- setNames(
  variantstring::allowed_amino_acids()$IUPAC_amino_acid_code,
  variantstring::allowed_amino_acids()$three_letter_code
)
df <- df %>%
  mutate(mutation_name = str_remove(mutation_name, "^k13-")) %>%
  mutate(mutation_name = str_replace_all(mutation_name, aa_lookup))

# function to split duplicated entries
split_rows_by_and <- function(data, column) {
  # Use tidyr's separate_rows to split the specified column
  data %>%
    separate_rows(!!sym(column), sep = "\\s*&&\\s*")
}

# combine per sample, split any multiple loci samples as assume these are unphased
correct_st <- df %>% filter(!is.na(genotype)) %>% 
  group_by(district, year, sample_id) %>% 
  summarise(mut = paste(mutation_name[genotype==1], collapse = "&&")) %>% 
  split_rows_by_and("mut") %>% 
  group_by(district, year) %>% 
  mutate(n = length(unique(sample_id))) %>% 
  group_by(district, year, mut, n) %>% 
  summarise(x = n()) %>% 
  group_by(district, year) %>% 
  mutate(mut = replace(mut, mut == "", "wildtype")) %>% 
  mutate(prev = x/n) %>% 
  rename(site = district) %>% 
  mutate(pmid = "33146722") %>% 
  mutate(year = as.character(year))

# left join with the data in k13wwdf (no year as we have found new years for some sites)
correct_k13st <- k13wwdf %>% filter(pmid=="33146722") %>% 
  select(-x,-n,-prev, -mut, -year) %>% unique() %>% left_join(correct_st, .)

k13wwdf <- bind_rows(k13wwdf %>% filter(pmid!="33146722" | is.na(pmid)), 
                     correct_k13st)

# ---------------------------------------------------- o
## 2.2 Identify possible issues  ----
# ---------------------------------------------------- o

# assign to uuids 
assign_ids <- function(x) {
  x %>%
    group_by(across(c(-x, -n, -prev, -mut))) %>%
    mutate(uuid = cur_group_id()) %>%
    ungroup %>% 
    group_by(across(c(-x, -prev, -mut))) %>%
    mutate(nid = cur_group_id()) %>%
    ungroup() %>% 
    mutate(iid = seq_len(n()))
}

# ---------------------------------------------------- o
### 2.2.1 Creating list of possible concerns based on sum(x) != n  ----
# ---------------------------------------------------- o

# assign to uuids 
k13wwdf <- assign_ids(k13wwdf)

# Check 1 to scan over possible issues where sum of x equals n
iss1 <- k13wwdf %>%
  group_by(nid) %>%
  mutate(xs = sum(x)) %>% 
  filter(xs>n) %>% 
  arrange(uuid) 

# There are 62 nids where this fails... FML
nid_concerns <- unique(iss1$nid)

# Firstly, we can ignore any that have pmid 33146722 as we dealt with these before
nid_concerns <- k13wwdf %>% filter(nid %in% nid_concerns) %>% filter(pmid!="33146722" | is.na(pmid)) %>% 
  pull(nid) %>% unique()

# Great now only 40 nids

# Secondly, if there are only 3 entries, WT and two mutants, and the sum of x for the
# highest two observations equals n, then there is a mixed infection. 
# The only way to know if this is a double mutant or mixed infection is to look at the paper

# Most will likely not have been phased as PCR but to be correct should look into these
nids_maybe_okay <- k13wwdf %>% filter(nid %in% nid_concerns) %>% split(.$nid) %>% 
  lapply(function(x){x %>% mutate(n2 = n())}) %>% do.call(rbind, .) %>% filter(n2 == 3) %>% 
  group_by(nid) %>% summarise(g = sum(sort(x, decreasing = TRUE)[1:2]) == unique(n)) %>% 
  filter(g) %>% pull(nid)

# Okay nids based on paper's methods, i.e. they could not have detected double mutant
# read through the methodology for these papers and figured out which are mixed mutants and which are explicitly double
# nids_not_okay includes those with doubles to be fixed
# more details here: https://docs.google.com/document/d/1SwLdHDbkhi-eLvf9BKe2fWlk0DfDdOG4tNCppdzTNf0/edit?usp=sharing
nids_okay <- c(26, 227, 
               # Gina added -- assuming this means that these are definitely mixed
               939, 
               1118, 1121) # both studies reported no double k13

nids_not_okay <- c(192, 352, 353, 354, 738, 1189)
# c("192"="R645T_E668K",
# "352" = "M579T_N657H", # 12 double, 1 with M579T only
# "353" = "M579T_N657H", # 16 double, 1 with M579T only - same pmid as above different year,
# "354" = "M579T_N657H", # all double mutant - again, same pmid
# "738" = "F446I_P574L", # 1 double, 2 single
# "1189" = "S477F_T677I") # 1 double mutant

# fix these ones manually
k13wwdffix <- k13wwdf %>%
  filter(nid %in% nids_not_okay) %>%
  split(.$nid)

k13wwdffix$`192` <- k13wwdffix$`192` |> 
  dplyr::mutate(mut = c("wildtype", "R645T_E668K", "R645T_E668K")) |> 
  distinct(x, n, gene, mut, .keep_all = TRUE) 

k13wwdffix$`352` <- k13wwdffix$`352` |> 
  dplyr::mutate(mut = c("wildtype", "M579T", "M579T_N657H"),
                x = c(26, 1, 12)) |> 
  distinct(x, n, gene, mut, .keep_all = TRUE) 

k13wwdffix$`353` <- k13wwdffix$`353` |> 
  dplyr::mutate(mut = c("wildtype", "M579T", "M579T_N657H"),
                x = c(28, 1, 16)) |> 
  distinct(x, n, gene, mut, .keep_all = TRUE) 

# all double mutants
k13wwdffix$`354` <- k13wwdffix$`354` |> 
  dplyr::mutate(mut = c("wildtype", "M579T_N657H", "M579T_N657H"),
                x = c(29, 22, 22)) |> 
  distinct(x, n, gene, mut, .keep_all = TRUE) 

k13wwdffix$`738` <- k13wwdffix$`738` |> 
  dplyr::mutate(mut = c("wildtype", "F446I", "F446I_P574L"),
                x = c(2, 2, 1)) |># 2 WT, 2 single, 1 double
  distinct(x, n, gene, mut, .keep_all = TRUE) 

k13wwdffix$`1189` <- k13wwdffix$`1189` |> 
  dplyr::mutate(mut = c("wildtype", "S477F_T677I", "S477F_T677I")) |>
  distinct(x, n, gene, mut, .keep_all = TRUE) 

# combine back in with the rest of k13wwdf
k13wwdffix <- do.call(rbind, k13wwdffix) 
k13wwdf <- rbind(filter(k13wwdf, (nid %in% nids_not_okay) == FALSE), k13wwdffix) # checked we have same # of nids

# Last remove these from the list of concerns 
nid_concerns <- nid_concerns[!(nid_concerns %in% nids_maybe_okay)]

# ---------------------------------------------------- o
### 2.2.2 3 rows where sum(x) != n and look weird ----
# ---------------------------------------------------- o

# There are 29 uuids where this fails... FML

# Lets look at the remaning 3 ones
nids_d_not_okay <- k13wwdf %>% filter(nid %in% nid_concerns) %>% split(.$nid)  %>% 
  lapply(function(x){x %>% mutate(n2 = n())}) %>% do.call(rbind, .) %>% filter(n2 == 3) %>% group_by(nid) %>% 
  summarise(g = sum(sort(x, decreasing = TRUE)[1:2]) != unique(n)) %>% filter(g) %>% pull(nid)

# these actually look like data entry issues

# https://pubmed.ncbi.nlm.nih.gov/31132213/
new_nid <- rbind(k13wwdf %>% filter(site == "Epe" & pmid == "31132213"),
                 k13wwdf %>% filter(site == "Epe" & pmid == "31132213" & mut == "Q613H"))
new_nid$mut[4] <- "D464N"
new_nid$x[1] <- 190
new_nid$x[new_nid$mut == "Q613H"] <- 3
new_nid$n <- 195
new_nid$prev <- new_nid$x/new_nid$n

# update k13wwdf
k13wwdf <- bind_rows(k13wwdf[-which(k13wwdf$site == "Epe" & k13wwdf$pmid == "31132213"),], 
                     new_nid)
nid_concerns <- nid_concerns[which(nid_concerns != 833)]

# https://pubmed.ncbi.nlm.nih.gov/31132213/
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-018-2625-6/tables/1
# This study did actually work out whether infections were mixed or double mutations
# FML - reading through, Thateng had 1 C580Y_Y493H and Khong had 1 C580Y_R539T

# But below is correct correction for this site
new_nid2 <- rbind(k13wwdf %>% filter(site == "Pathoumphone DH, Pathoumphone, Champasak" & pmid == "30587196"))
new_nid2$x[new_nid2$mut == "wildtype"] <- 90
new_nid2$x[new_nid2$mut == "R539T"] <- 33
new_nid2$x[new_nid2$mut == "C580Y"] <- 219
new_nid2$prev <- new_nid2$x/new_nid2$n

# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-018-2625-6/tables/1
# see also additional file 8 with msp1 and msp2
# tot = 112; C580Y = 16; Y493H = 26; C580Y_R538T = 1 (double not mixed - double and mixed are conflated in table 1)
new_nid2b <- rbind(k13wwdf %>% filter(site == "Thateng DH, Thateng, Sekong" & pmid == "30587196"))
new_nid2b$mut[which(new_nid2b$mut == "C580C/Y")] <- "C580Y_Y493H"

new_nid2b <- rbind(k13wwdf %>% filter(site == "Thateng DH, Thateng, Sekong" & pmid == "30587196"))
# I think this was misextracted as mixed vs double. cross-referenced tab 1 and AF 8
new_nid2b$mut[which(new_nid2b$mut == "C580C/Y")] <- "C580Y_Y493H"
new_nid2b$prev <- new_nid2b$x/new_nid2b$n

new_nid2c <- rbind(k13wwdf %>% filter(site == "Khong DH, Khong, Champasak" & pmid == "30587196"))
new_nid2c$x[new_nid2$mut == "wildtype"] <- 20
# duplicate a row to add a double mutant 
new_nid2c <- new_nid2c %>%
  bind_rows(slice(., 2))
new_nid2c$mut[4] <- "C580Y_R539T"
new_nid2c$x[4] <- 1
new_nid2c$prev <- new_nid2c$x/new_nid2c$n

# update k13wwdf
k13wwdf <- bind_rows(k13wwdf[-which(k13wwdf$site %in% c("Pathoumphone DH, Pathoumphone, Champasak",
                                                        "Thateng DH, Thateng, Sekong",
                                                        "Khong DH, Khong, Champasak") & k13wwdf$pmid == "30587196"),], 
                     rbind(new_nid2, new_nid2b, new_nid2c))
nid_concerns <- nid_concerns[which(nid_concerns != 580)] # I think this was a typo - the nid of above is 580

# https://pubmed.ncbi.nlm.nih.gov/32795367/
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-020-03358-7/tables/1
# was sanger so no way to know how so leave as okay
nid_concerns <- nid_concerns[which(nid_concerns != 481)]

# ---------------------------------------------------- o
### 2.2.3 remaining sum(x) != n and look weird ----
# ---------------------------------------------------- o

# # pmid 99999999 doesn't seem to exist - I cannot find this paper. I'm going to exclude this paper. this removes a lot of the issues
nid_removed <- k13wwdf |>
  dplyr::filter(pmid == 99999999) |>
  dplyr::filter(nid %in% nid_concerns) |>
  dplyr::pull(nid) |> unique()
# remove data that corresponds to this study and has an nid causing issues
# TODO: check with OJ if this was actually the correct thing to do 
k13wwdf <- k13wwdf |>
  dplyr::filter((nid %in% nid_removed) == FALSE)

nid_concerns <- nid_concerns[which((nid_concerns %in% nid_removed) == FALSE)]

k13wwdffix <- k13wwdf %>% filter(nid %in% nid_concerns) %>% split(.$nid)

# prioritise African studies for now
# nids in Africa: 188, 186, 640, 644, 234, 274, 276, 452, 891, 895

# pmid 26667053; nid 188 186; incorrect extraction fixing
k13wwdffix$`186`$x[which(k13wwdffix$`186`$mut == "S477Y")] <- 0
k13wwdffix$`186` <- k13wwdffix$`186` |>
  dplyr::filter(x > 0) |>
  dplyr::mutate(prev = x/n)

# fixing incorrect extractions and double mutant
k13wwdffix$`188`$x[which(k13wwdffix$`188`$mut == "A578S")] <- 1
k13wwdffix$`188`$mut[which(k13wwdffix$`188`$mut == "I526M")] <- "I526M_A578S"
k13wwdffix$`188`$prev <- k13wwdffix$`188`$x/k13wwdffix$`188`$n

# pmid 29436339 nid 234 - two of the mutants were actually mixed infection, no doubles

# nid = 640 or 644 - prev sums > 1 so I think some of these are mixed
# impossible to disentangle what's happening

# 28797235; nid 274 and 276 = mixed markers 
# 452 is mixed also

# nid 891 https://academic.oup.com/jid/article/225/8/1411/6314293#supplementary-data
# incorrect extraction
k13wwdffix$`891`$x[which(k13wwdffix$`891`$mut == "Q661E")] <- 1
k13wwdffix$`891`$x[which(k13wwdffix$`891`$mut == "P667S")] <- 1
k13wwdffix$`891`$prev <- k13wwdffix$`891`$x/k13wwdffix$`891`$n

# nid 895 implies mixed infections vs double mutants - unclear from methods

# studies outside of Africa
# 665 - 668, 737 and 79 all have the same pmid and all mixed so leave 25927592

# 493 27001814
# https://journals.asm.org/doi/10.1128/aac.02370-15
# states one of the parasites had two mutants but doesn't specify haplotype


# 1180 27332904 https://pmc.ncbi.nlm.nih.gov/articles/instance/4955562/bin/NIHMS802595-supplement-app.pdf
# V520I & T474I in the same sample
k13wwdffix$`1180`$mut[which(k13wwdffix$`1180`$mut == "T474I")] <- "T474I_V520I"
k13wwdffix$`1180`$x[which(k13wwdffix$`1180`$mut == "V520I")] <- 0
k13wwdffix$`1180` <- k13wwdffix$`1180` |>
  dplyr::filter(x > 0) |>
  dplyr::mutate(prev = x/n)

# 592 30587196
k13wwdffix$`592`$x[which(k13wwdffix$`592`$mut == "wildtype")] <- 29
k13wwdffix$`592`$prev <- k13wwdffix$`592`$x/k13wwdffix$`592`$n

# 319 & 320  32393498 - talks about mixed infections and no mention of doubles 

nid_fixed <- c(186, 188, 891, 1180, 592)
nid_checked <- c(234, 640, 644,274, 276, 452, 895, 665:668, 737, 79, 493, 319, 320)

nid_concerns <- nid_concerns[intersect(which(!nid_concerns %in% nid_fixed), which(! nid_concerns %in% nid_checked))]
if(length(nid_concerns) == 0) {
  print("nid_concerns all resolved in this section")
}

# bind everything back together
k13wwdffix <- do.call(rbind, k13wwdffix) 
k13wwdf <- rbind(filter(k13wwdf, (nid %in% c(nid_fixed, nid_checked) == FALSE)), k13wwdffix) # checked we have same # of nids

# ---------------------------------------------------- o
### 2.2.4 TODO: Decide how much we care about looking for where sum(x) == n, but they have not reported double mutants----
# ---------------------------------------------------- o


# ---------------------------------------------------- o
### 2.3.1 Creating list of possible concerns based on sum(x) < n  ----
# ---------------------------------------------------- o

iss2 <- k13wwdf %>%
  group_by(nid) %>%
  mutate(xs = sum(x)) %>%
  filter(xs < n) %>%
  arrange(uuid) 
# TODO: decide what to do about the N-terminal SNPs and stop codons in WWARN

# write a function to add a row
add_a_row <- function(df, x_new, mut_new) {
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

# There are 23 nids where this fails... FML
nid_lconcerns <- unique(iss2$nid)

k13wwdf_lfix <- k13wwdf %>%
  filter(nid %in% nid_lconcerns) %>%
  split(.$nid)

# TODO: Go through these and work out what the denominator should actually be
# going through each pmid at a time because many from same study

# however https://journals.asm.org/doi/full/10.1128/aac.00802-19?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org
# nid 250; pmid 31427297 - specific issue with "forest region" in paper
# 10 mutations plus WT -- 11 rows

k13wwdf_lfix$`250` <- k13wwdf_lfix$`250` %>%
  bind_rows(slice(., 2:3))

fixed_cols <- data.frame(
  x = c(199,2,2,2,13,5,2,1,2,1),
  mut = c("wildtype", "N408S", "L422F",
          "S466N", "G496G", # removed the stop codon
          "R515I", "R539I", "C542W",
          "E606D", "K607E")
)

k13wwdf_lfix$`250` <- k13wwdf_lfix$`250` %>% 
  dplyr::mutate(x = fixed_cols$x,
                mut = fixed_cols$mut,
                prev = x/n)

# nid 414 pmid 25691626 - typo 
k13wwdf_lfix$`414` <- k13wwdf_lfix$`414` %>% 
  dplyr::mutate(x = c(15,1),
                prev = x/n)

# nid 648 pmid 28249583 - add wildtype
k13wwdf_lfix$`648` <- k13wwdf_lfix$`648` %>%
  bind_rows(slice(.,1)) # add a new top row
k13wwdf_lfix$`648`$x[1] <-  k13wwdf_lfix$`648`$n[1] - sum(k13wwdf_lfix$`648`$x[2:5])
k13wwdf_lfix$`648`$mut[1] <- "wildtype"
k13wwdf_lfix$`648`$prev <- k13wwdf_lfix$`648`$x/k13wwdf_lfix$`648`$n

# pmid 29582728 nid 458 and 466
k13wwdf_lfix$`458` <- add_a_row(k13wwdf_lfix$`458`, 1,
                                "A676V")

# nid 466 n is wrong, missing mutations and ones extracted are also wrong...
kisumu_cols <- data.frame(
  x = c(32,1,1,2,4,1),
  mut = c("wildtype", "E509G", "V568G",
          "D584Y", "N585K", "T677R")
)
k13wwdf_lfix$`466` <- k13wwdf_lfix$`466` %>%
  bind_rows(slice(.,5)) %>%
  dplyr::mutate(x = kisumu_cols$x,
                mut = kisumu_cols$mut, 
                n = 41, 
                prev = x/n)

# nid 705 pmid 25537878 
# n and tot mutants correct
# based on the text (plus a typo FML) there is one D516Y mut
k13wwdf_lfix$`705` <- add_a_row(k13wwdf_lfix$`705`, 1,
                                "D516Y")

# 743 pmid 31833468 - missing one mutation
k13wwdf_lfix$`743` <- add_a_row(k13wwdf_lfix$`743`, 1,
                                "G453D")

# 877 34270452 -- I can't even find this site in study
# remove this entry
k13wwdf_lfix$`877` <- NULL

# 1136 34551228 Table S7 https://www.nejm.org/doi/suppl/10.1056/NEJMoa2101746/suppl_file/nejmoa2101746_appendix.pdf
k13wwdf_lfix$`1136` <- add_a_row(k13wwdf_lfix$`1136`,
                                 x_new = 1, mut_new = "V661I")

# 682 - pmid 99999999 cannot find. 1 short of equal
# 683 is missing loads - leave 682 and omit 683

k13wwdf_lfix$`683` <- NULL

# 710 pmid 99999999 - no WT so assume missing val is WT
k13wwdf_lfix$`710` <- add_a_row(k13wwdf_lfix$`710`, 1,
                                "wildtype")

# 381 99% WT but cannot find out the missing mutant

# nids resolved so far
nid_lresolved <- c(250, 414, 648, 466, 458, 705, 743, 
                   877, 1136, 683, 710) |> as.character()

# add the resolved nids back into WWARN
k13wwdf_lfixed <- do.call(rbind, k13wwdf_lfix[nid_lresolved])

# pmid 31591113 nid 442:450
# all data in frequencies
# for now, ignore N-terminal SNPs 
# due to the different # sequenced, different ns.
# TODO: decide how to handle this

# it looks like 957 and all with source 
# "Suttipat Srisutham, Nguyen Than Thuy-Nhien, Ranitha Vongpromek, Teeradet Khomvarn, Mayfong Mayxay, Olivo Miotto, Francois Nosten, Frank Smithuis, Rob van der Pluijm, Lorenz von Seidlein, Carol H Sibley, Philippe J. Guérin, Nicholas P. J. Day, Arjen Dondorp, Mehul Dhorda, Mallika Imwong"
# may actually come from Das study: https://link.springer.com/article/10.1186/s12936-022-04095-9
# this is based on shared authors, study year and location
# looking at the Das study, almost all of the data is missing... ugh
# fixing an additional issue 
# Das study (PMID 35279140) is almost entirely missing (lots of datapoints in Chad)
# source "Suttipat Srisutham, Nguyen Than Thuy-Nhien, Ranitha Vongpromek, Teeradet Khomvarn, Mayfong Mayxay, Olivo Miotto, Francois Nosten, Frank Smithuis, Rob van der Pluijm, Lorenz von Seidlein, Carol H Sibley, Philippe J. Guérin, Nicholas P. J. Day, Arjen Dondorp, Mehul Dhorda, Mallika Imwong"
# looks like the same markers but also with missing info
# removed Srisutham et al. and full extracted Das et al.
# need to figure out a way to do this that won't break everything...

# https://link.springer.com/article/10.1186/s12936-022-04095-9/tables/2
fix_das <- data.frame(
  year = c(rep("2016", 10), rep("2017", 14)),
  x = c(6, 76, 8, 1, 1, 1, 1, 1, 6, 56,
        1, 1, 2, 83, 1, 1, 1, 17,1, 1, 1, 5, 1, 76),
  n = c(rep(157, 10), rep(192, 14)),
  mut = c("K189N", "K189T", "K189T/K",
          "K189T_I354V", # haplotypes
          "K189T_N197D",
          "L258M", "N197D", "Q633R", 
          "R255K", "wildtype", # 2016 mutations
          "A578S", "K189N", "K189N/K", "K189T",
          "K189T_N197D/N",
          "K189T_V636A/V",
          "K189T_W660C",
          "K189T/K",
          "K189T/K_N197D/N",
          "N195D/N",
          "N195K/N",
          "R255K", "S213G/S", "wildtype") # 2017 muts
) |> dplyr::mutate(pmid = "35279140")


k13wwdf_das <- k13wwdf |>
  dplyr::filter(pmid == 35279140) |>
  dplyr::select(-c(year, x, n, mut)) 

# between studies had nids 
k13wwdf_das <- rbind(k13wwdf_das,
                     k13wwdf_das[rep(1, 22),]) |>
  dplyr::mutate(x = fix_das$x,
                year = fix_das$year,
                n = fix_das$n,
                mut = fix_das$mut) |>
  dplyr::relocate(names(k13wwdf)) |>
  dplyr::mutate(x = as.integer(x),
                n = as.integer(n)) |>
  dplyr::select(-iid) |> # need to add this back in after its added to rest of df
  dplyr::mutate(nid = if_else(year == "2016", 952, 953)) |> # 952:959 is original nids
  dplyr::mutate(pmid = "35279140")
  
# add Das back in 
k13wwdf <- k13wwdf |>
  dplyr::select(-iid) |>
  dplyr::filter(pmid != "35279140") |>
  dplyr::filter(source != "Suttipat Srisutham, Nguyen Than Thuy-Nhien, Ranitha Vongpromek, Teeradet Khomvarn, Mayfong Mayxay, Olivo Miotto, Francois Nosten, Frank Smithuis, Rob van der Pluijm, Lorenz von Seidlein, Carol H Sibley, Philippe J. Guérin, Nicholas P. J. Day, Arjen Dondorp, Mehul Dhorda, Mallika Imwong") |>
  bind_rows(k13wwdf_das) |>
  dplyr::mutate(iid = seq_len(n()))

# TODO: still need to resolve pmid 31591113 - data still included but very wrong 
k13wwdf <- k13wwdf |>
  dplyr::filter(pmid != "31591113")


# ---------------------------------------------------- o
## 2.4 Now format after having assumed we have corrected all possible issues  ----
# ---------------------------------------------------- o

# this should now be TRUE
nrow(k13wwdf) == (k13wwdf %>% select(-uuid, -nid, -iid) %>%
                    group_by(across(c(-x, -n, -prev))) %>% 
                    group_indices() %>% unique() %>% length)

# and sort out the columns etc
k13ww_final_res_df <- k13wwdf %>% 
  select(-uuid, -nid, -iid) %>%
  mutate(mut = if_else(mut == "wildtype", "WT", mut)) %>% 
  ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gene_mut = wwarn_format_k13_for_stave(mut)) %>% 
  ungroup() %>% 
  mutate(mut = if_else(mut == "WT", "WT", substring(mut, 2))) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, gene_mut, database, pmid, url, source)

# TODO: update this so that it also includes the partner drug data
saveRDS(k13ww_final_res_df, here::here("analysis/data-derived/wwarn_res.rds"))

# ---------------------------------------------------- o
# 3. Sort PD Names ----
# ---------------------------------------------------- o

pdww <- readRDS("analysis/data-derived/pdww.rds") 

# sort names as wanted
pdwwdf <- pdww %>%
  rename(admin_0 = country) %>%
  mutate(iso3c = countrycode::countrycode(admin_0, "country.name.en", "iso3c")) %>%
  rename(long = lon) %>%
  mutate(study_start_year = as.integer(`study start`)) %>%
  mutate(study_end_year = as.integer(`study end`)) %>%
  mutate(year = as.integer(study_start_year + ((study_end_year - study_start_year)/2))) %>%
  mutate(n = as.integer(tested)) %>%
  mutate(x = as.integer(present)) %>%
  mutate(mix = as.integer(`Mixed present`)) %>%
  mutate(prev = as.integer(percentage)/100) %>%
  mutate(gene = "k13") %>%
  rename(mut = `marker group`) %>%
  mutate(pmid = notes) %>%
  mutate(url = `publication URL`) %>%
  mutate(database = "WWARN") %>%
  rename(source = author) %>%
  mutate(gene = NA) %>%
  mutate(gene = replace(gene, grepl("crt", mut), "crt")) %>%
  mutate(gene = replace(gene, grepl("mdr1", mut), "mdr1")) %>%
  mutate(gene = replace(gene, grepl("pm2", mut), "pfpm23")) %>%
  mutate(stuid = paste0(`study Id`, `site number`)) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source,
         prev, mix, stuid) %>%
  unique() %>%
  mutate(rowid = seq_len(n()))

# Growing List of typos in their data entry
# these have come from troubleshooting the cleans below and going back to papers

# All the mdr1 typos... ----------------
# A lot of from rounding errors from prev * n calculations
# Others are from looking at barplots / pie charts without any numbers available...
pdwwdf$x[which(pdwwdf$rowid == 7573)] <- pdwwdf$n[which(pdwwdf$rowid == 7573)]
pdwwdf$x[which(pdwwdf$rowid == 6995)] <- pdwwdf$n[which(pdwwdf$rowid == 6995)]
pdwwdf$x[which(pdwwdf$rowid == 1630)] <- 27
pdwwdf$x[which(pdwwdf$rowid == 64)] <- 57
pdwwdf$x[which(pdwwdf$rowid == 1372)] <- 0
pdwwdf$x[which(pdwwdf$rowid == 1377)] <- 58
pdwwdf$x[which(pdwwdf$rowid == 1380)] <- 2
pdwwdf$x[which(pdwwdf$rowid == 1653)] <- 34
pdwwdf$n[match(5131:5133, pdwwdf$rowid)] <- 74
pdwwdf$n[match(6573:6575, pdwwdf$rowid)] <- 22
pdwwdf$x[match(6573:6575, pdwwdf$rowid)] <- c(0,20,2)
pdwwdf$x[match(2806:2808, pdwwdf$rowid)] <- c(39,23,46)
pdwwdf$x[which(pdwwdf$rowid == 2742)] <- 8
pdwwdf$x[match(2914:2916, pdwwdf$rowid)] <- c(548,112,85)
pdwwdf$n[match(2983:2985, pdwwdf$rowid)] <- 19 # typo in the actual paper most likely
pdwwdf$x[which(pdwwdf$rowid == 3007)] <- 57
pdwwdf$x[which(pdwwdf$rowid == 3022)] <- 45
pdwwdf$x[which(pdwwdf$rowid == 3019)] <- 84
pdwwdf$x[which(pdwwdf$rowid == 1136)] <- 5
pdwwdf$x[which(pdwwdf$rowid == 5045)] <- 83
pdwwdf$x[which(pdwwdf$rowid == 4676)] <- 41
pdwwdf$x[which(pdwwdf$rowid == 5914)] <- 16
pdwwdf$n[match(4497:4499, pdwwdf$rowid)] <- 106
pdwwdf$n[match(4503:4505, pdwwdf$rowid)] <- 107
pdwwdf$x[which(pdwwdf$rowid == 6397)] <- 10
pdwwdf$x[which(pdwwdf$rowid == 6147)] <- 140
pdwwdf$x[which(pdwwdf$rowid == 5484)] <- 176
pdwwdf$x[which(pdwwdf$rowid == 5485)] <- 17
pdwwdf$x[which(pdwwdf$rowid == 5552)] <- 23
pdwwdf$n[match(c(5551,5552,5554), pdwwdf$rowid)] <- 145
pdwwdf$n[match(c(5458,5459,5461), pdwwdf$rowid)] <- 110
pdwwdf$x[which(pdwwdf$rowid == 5642)] <- 14
pdwwdf$n[which(pdwwdf$rowid == 5689)] <- 35

# the next block of typos are where prevalence is reported using mixed infections
# so we need to correct these
pdwwdf$x[which(pdwwdf$rowid == 13)] <- 42 # F is not Y so assign to N
pdwwdf$x[which(pdwwdf$rowid == 318)] <- 44
pdwwdf$x[which(pdwwdf$rowid == 718)] <- 114 # typo in their paper
pdwwdf$x[which(pdwwdf$rowid == 720)] <- 103

# typos again
pdwwdf$x[which(pdwwdf$rowid == 2475)] <- 5
pdwwdf$x[match(c(7102,7103), pdwwdf$rowid)] <- c(66, 94)
pdwwdf$n[match(c(7102,7103), pdwwdf$rowid)] <- 160
pdwwdf$x[which(pdwwdf$rowid == 7714)] <- 42
pdwwdf$x[which(pdwwdf$rowid == 7737)] <- 47
pdwwdf$n[match(c(8283,8286), pdwwdf$rowid)] <- 119
pdwwdf$x[match(c(8283,8286), pdwwdf$rowid)] <- c(7,112)
pdwwdf$x[which(pdwwdf$rowid == 4451)] <- 14 # F is not Y so assign to N

# duplicate studies
pdwwdf <- pdwwdf %>% filter(pmid != 80002356)
pdwwdf <- pdwwdf %>% filter(stuid != "265972541")
pdwwdf <- pdwwdf %>% filter(rowid != 1526) # duplicate error

# All the crt typos... ----------------
pdwwdf$n[which(pdwwdf$rowid == 8822)] <- 40
pdwwdf$n[which(pdwwdf$rowid == 8754)] <- 68
pdwwdf <- pdwwdf %>% filter(rowid != 308) # duplicate error
pdwwdf$x[which(pdwwdf$rowid == 7015)] <- 0.5 # mixed sample and swapping to mut type here
pdwwdf$x[which(pdwwdf$rowid == 5763)] <- 0 # mixed sample double accounted for
pdwwdf$mut[which(pdwwdf$rowid == 7015)] <- "pfcrt 76T" # mixed sample and swapping to mut type here

# one study's fixes
pdwwdf$x[match(c(1695,1696,1697,1698,1699,1700,2622,2623,2624,2625,2626,5650,5651,5652,5653,5654), pdwwdf$rowid)] <-
  c(35, 15, 35, 15, 0, 0, 1.5, 48.5, 1.5, 0, 48.5, 0.5, 37.5, 0.5, 0, 37.5)

# one study's fixes
pdwwdf$x[match(c(1608,1609,1610,1611,1612,1616,1843,1844,1845,1846,3105,3106,3107), pdwwdf$rowid)] <-
  c(27.5, 0, 21.5, 2, 29.5, 21.5, 3,0,4.5,25.5, 6,0,23)

# one study's fixes
pdwwdf$x[match(c(2792, 2816, 3016), pdwwdf$rowid)] <- c(33, 56, 73)

# one study's fixes
pdwwdf$x[match(c(1009,1010,1011,1015,1016,1017,8457,8458,8459), pdwwdf$rowid)] <-
  c(30.5, 1.25, 1.25, 3.5, 0, 9.5, 12, 5, 3)

# one study's fixes
# WWARN missed some rows here so easiest way is to just adding the counts
# frommissing haplotypes to the same 76 types
pdwwdf$x[match(c(807, 808), pdwwdf$rowid)] <- c(67, 7)

# one study's fixes
# WWARN missed some rows here so easiest way is to just adding the counts
# from missing haplotypes to the same 76 types
pdwwdf$x[match(c(3965), pdwwdf$rowid)] <- 3

# duplicate studies
pdwwdf <- pdwwdf %>% filter(pmid != 80002356)
pdwwdf <- pdwwdf %>% filter(stuid != "265972541")
pdwwdf <- pdwwdf %>% filter(rowid != 3906) # duplicate error
pdwwdf <- pdwwdf %>% filter(rowid != 3907) # duplicate error
pdwwdf <- pdwwdf %>% filter(rowid != 3957) # duplicate error
pdwwdf <- pdwwdf %>% filter(rowid != 3958) # duplicate error

# OTHER TYPOS
pdwwdf$study_start_year[match(c(3485,3486,9148,9161,9171,9205), pdwwdf$rowid)] <- 2003
pdwwdf$study_end_year[match(c(3485,3486,9148,9161,9171,9205), pdwwdf$rowid)] <- 2006
pdwwdf$study_start_year[match(c(6894,6895,6896,6897,6898,6899,6900,6901,6902,6903), pdwwdf$rowid)] <- 2011
pdwwdf$year[match(c(6894,6895,6896,6897,6898,6899,6900,6901,6902,6903), pdwwdf$rowid)] <- 2012

# split these out into each marker type
pdwwspl <- pdwwdf %>% split(pdwwdf$gene)

# ---------------------------------------------------- o
# 4. Sort CRT ----
# ---------------------------------------------------- o

# Grab 76 etc
pdcrt <- pdwwspl$crt %>%
  filter(grepl("76", mut))

# first put into study time lat points for looking at
pdcrt <- pdcrt %>%
  group_by(across(c(-x,-prev, -mut, -rowid, -mix))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  dplyr::distinct() # a lot of the issues with sum(x) =/= n are due to duplicate rows
# need to figure out how they have entered the EH 72

### TYPE 1 ------------------------------
# there is a group of entries where all x sum to n
# so for these convert the EHs down by grouping by mutation
# marker after converting the 72 muts to common type

# type 1 when there are multiple markers reported so we sum mutant and mixed
pdcrtspl1 <- pdcrt %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  mutate(mut = replace(mut, grepl("CxxxK", mut), "pfcrt K76")) %>% # convert haplotypes
  mutate(mut = replace(mut, grepl("CxxxT", mut), "pfcrt 76T")) %>% # convert haplotypes
  mutate(mut = replace(mut, grepl("SxxxT", mut), "pfcrt 76T")) %>% # convert haplotypes
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  ungroup() %>%
  group_by(across(c(-x, -n, -mut))) %>%
  mutate(l = n()) %>%
  filter(l > 1 | (l==1 & all(mut!="pfcrt K76"))) %>% # catch for when only one marker is reported
  summarise(x = ifelse(any(mut == "pfcrt 76T"), x[mut == "pfcrt 76T"], 0) +
              ifelse(any(mut == "pfcrt 76K/T"), x[mut %in% "pfcrt 76K/T"], 0), # add the mixed and mutant to get prev
            n = unique(n),
            prev = x/n) %>%
  mutate(mut = "crt_76T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

### TYPE 2 ------------------------------
# and type 2 is the catch for when it is only WT
pdcrtspl2 <- pdcrt %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  mutate(mut = replace(mut, grepl("CxxxK", mut), "pfcrt K76")) %>%
  mutate(mut = replace(mut, grepl("CxxxT", mut), "pfcrt 76T")) %>%
  mutate(mut = replace(mut, grepl("SxxxT", mut), "pfcrt 76T")) %>%
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  ungroup() %>%
  group_by(across(c(-x, -n, -mut))) %>%
  mutate(l = n()) %>%
  filter(l==1 & all(mut=="pfcrt K76")) %>% # catch for when only one marker and its WT is reported
  summarise(x = n - x,
            n = unique(n),
            prev = x/n) %>%
  mutate(mut = "crt_76T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

# now to focus on when sum of x does not equal n
# checked this covers all values of pdcrt$mut

mut_loc <- "pfcrt 76T"
wt_loc <- "pfcrt K76"
mix_loc <- "pfcrt 76K/T"

other_loc <- c("pfcrt 72-76 CxxxK", "pfcrt 72-76 CxxxT", "pfcrt 72-76 SxxxT")

### TYPE 3 ------------------------------
# The EH mutations for a number of samples are just extra information
# determined by filtering these out and rechecking if sum of x equals n
# These for these groups it is the same as above
pdcrtspl3 <- pdcrt %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(!xn) %>%
  filter((mut %in% c(mut_loc,wt_loc,mix_loc))) %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  ungroup() %>%
  group_by(across(c(-x, -n, -mut))) %>%
  mutate(l = n()) %>%
  filter(l > 1 | (l==1 & all(mut!="pfcrt K76"))) %>% # catch for when only one marker is reported
  summarise(x = ifelse(any(mut == "pfcrt 76T"), x[mut == "pfcrt 76T"], 0) +
              ifelse(any(mut == "pfcrt 76K/T"), x[mut %in% "pfcrt 76K/T"], 0), # add the mixed and mutant to get prev
            n = unique(n),
            prev = x/n) %>%
  mutate(mut = "crt_76T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

# There are no examples here where l == 1 and they are pfcrt k76 so no need to switch the
# prev around as for type 2


### TYPE 4 ------------------------------
# And conversely the same as above but filtering by the EH types and checking
# for sum of x equal to n
# However, these are all the same uuids as TYPE 4 so ignore
pdcrt %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(!xn) %>%
  filter(!(mut %in% c(mut_loc,wt_loc,mix_loc))) %>%
  group_by(uuid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  filter(!(uuid %in% pdcrtspl3$uuid))

# We have the vast majority now based on uuid capture
rbind(pdcrtspl1,pdcrtspl2,pdcrtspl3) %>%
  group_by(uuid)

# the remaining ones do not have x sum equal to n at all, so
# likely are mixed infections reported differently, they only report
# one locus or they are typos...


### TYPE 5 ------------------------------
# these are single record uuids so WWARN only captured one marker
pdcrtspl5 <- pdcrt %>%
  filter(!(uuid %in% c(pdcrtspl1$uuid,pdcrtspl2$uuid,pdcrtspl3$uuid))) %>%
  group_by(uuid) %>%
  filter(n()==1) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxT", "pfcrt 76T")) %>%
  group_by(uuid) %>%
  mutate(prev = x/n) %>%
  select(names(pdcrtspl3)) %>%
  mutate(mut = "crt_76T") %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)


### TYPE 6 ------------------------------
# the others...

# then the remaining weird ones
# simple cases of just two mutations being off.
# These are mixed infections being reallocated
# TODO: Convert this code to be prevalence instead of frequency and check the 
# specific studies to ensure this is correct
# TODO: I think I actually need to go back to the studies themselves and fix this one
pdcrtspl6 <- rbind(
  pdcrt %>%
    filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid, pdcrtspl5$uuid))) %>%
    group_by(uuid) %>%
    filter(n() == 2) %>%
    filter(sum(x) < n[1]) %>%
    mutate(x = x + (n[1] - sum(x))/2), # TODO: convert to prev
  pdcrt %>%
    filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid, pdcrtspl5$uuid))) %>%
    group_by(uuid) %>%
    filter(n() == 2) %>%
    filter(sum(x) > n[1]) %>%
    mutate(x = x - (sum(x) - n[1])/2) # TODO: convert to prev
) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 SxxxT", "pfcrt 76T")) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxT", "pfcrt 76T")) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxK", "pfcrt K76")) %>%
  filter(mut == "pfcrt 76T") %>%
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  mutate(prev = x/n) %>%
  group_by(uuid) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

### TYPE 7 ------------------------------

# as before but by getting rid of the mixed to then give 2 and adjust
# These are mixed infections being reallocated
# TODO: Check these manually later that this is the case for these
pdcrtspl7 <-
  rbind(
    pdcrt %>%
      filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid, pdcrtspl5$uuid, pdcrtspl6$uuid))) %>%
      group_by(uuid) %>%
      filter(n() == 3) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 SxxxT", "pfcrt 76T")) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxT", "pfcrt 76T")) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxK", "pfcrt K76")) %>%
      filter(mut %in% c(mut_loc, wt_loc)) %>%
      filter(sum(x) < n[1]) %>%
      mutate(x = x + (n[1] - sum(x))/2), # TODO: convert to prev
    pdcrt %>%
      filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid, pdcrtspl5$uuid, pdcrtspl6$uuid))) %>%
      group_by(uuid) %>%
      filter(n() == 3) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 SxxxT", "pfcrt 76T")) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxT", "pfcrt 76T")) %>%
      mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxK", "pfcrt K76")) %>%
      filter(mut %in% c(mut_loc, wt_loc)) %>%
      filter(sum(x) > n[1]) %>%
      mutate(x = x - (sum(x) - n[1])/2) # TODO: convert to prev
  ) %>%
  filter(mut == "pfcrt 76T") %>%
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  mutate(prev = x/n) %>%
  group_by(uuid) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)


### TYPE 8 ------------------------------
# the others...

# then the remaining weird ones
# remove the EH haplotypes and then these collapse into
# 2 as above
# TODO: Check these manually later that this is the case for these
pdcrtspl8 <- rbind(
  pdcrt %>%
    filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid,
                         pdcrtspl5$uuid, pdcrtspl6$uuid, pdcrtspl7$uuid))) %>%
    group_by(uuid) %>%
    filter(mut %in% c(mut_loc, wt_loc)) %>%
    filter(n() == 2) %>%
    filter(sum(x) < n[1]) %>%
    mutate(x = x + (n[1] - sum(x))/2), # TODO: convert to prev
  pdcrt %>%
    filter(!(uuid %in% c(pdcrtspl1$uuid, pdcrtspl2$uuid, pdcrtspl3$uuid,
                         pdcrtspl5$uuid, pdcrtspl6$uuid, pdcrtspl7$uuid))) %>%
    group_by(uuid) %>%
    filter(mut %in% c(mut_loc, wt_loc)) %>%
    filter(n() == 2) %>%
    filter(sum(x) > n[1]) %>%
    mutate(x = x - (sum(x) - n[1])/2) # TODO: convert to prev
) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 SxxxT", "pfcrt 76T")) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxT", "pfcrt 76T")) %>%
  mutate(mut = replace(mut, mut == "pfcrt 72-76 CxxxK", "pfcrt K76")) %>%
  filter(mut == "pfcrt 76T") %>%
  group_by(across(c(-x, -n, -prev, -mix, -rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  mutate(prev = x/n) %>%
  group_by(uuid) %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source)

## Complete. Bring all together again-------------------

crtww_final_res_df <-
  list(pdcrtspl1,pdcrtspl2,pdcrtspl3,pdcrtspl5, pdcrtspl6, pdcrtspl7, pdcrtspl8) %>%
  lapply(function(x){
    x %>% ungroup %>% select(iso3c, admin_0, admin_1, site, lat, long,
                             year, study_start_year, study_end_year,
                             x, n, prev, gene, mut, database, pmid, url, source)
  }) %>% do.call(rbind,.) %>%
  ungroup %>%
  mutate(mut = replace(mut, mut == "pfcrt 76T", "crt_76T")) %>%
  dplyr::mutate(gene_mut = str_replace(mut, "_","-")) %>%
  dplyr::left_join(select(validated, c("gene_mut", "annotation")), by = "gene_mut") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gene_mut = clean_mutations(gene_mut)) %>%
  dplyr::relocate(names(k13ww_final_res_df))

# and the sanity check
(list(pdcrtspl1,pdcrtspl2,pdcrtspl3,pdcrtspl5, pdcrtspl6, pdcrtspl7, pdcrtspl8) %>%
    lapply(function(x){
      x %>% ungroup %>% select(iso3c, admin_0, admin_1, site, lat, long,uuid,
                               year, study_start_year, study_end_year,
                               x, n, prev, gene, mut, database, pmid, url, source)
    }) %>% do.call(rbind,.)  %>%
    pull(uuid) %>% length()) ==
  (pdcrt$uuid %>% unique()%>% length())


# ---------------------------------------------------- o
# 5. Sort MDR1 ----
# ---------------------------------------------------- o

# grab and fitler to just important SNPs
pdmdr1 <- pdwwspl$mdr1 %>%
  filter(grepl("184|86|copy|NFD|NxxxD|YYXXY|YYY", mut))

# first put into study time lat points for looking at
pdmdr1 <- pdmdr1 %>%
  group_by(across(c(-x,-prev, -mut, -rowid, -mix))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup()

## STEP 1: Work out how to sort out the markers ------------------------------

# Viewed through example splits like this to determine that we can filter out some more markers
# pdmdr1 %>% ungroup %>%
#   filter(stuid %in%
#            (pdmdr1 %>% filter(mut == "pfmdr1 YYY") %>% pull(stuid) %>% unique())) %>%
#   split(.$stuid)

#  [1,] "pfmdr1 86Y"
#  [2,] "pfmdr1 184F"
#  [3,] "pfmdr1 copy number >1"
#  [4,] "pfmdr1 NxxxD" : the same as"pfmdr1 N86" : Just EH info so remove
#  [5,] "pfmdr1 N86"
#  [6,] "pfmdr1 86N/Y"
#  [7,] "pfmdr1 Y184"
#  [8,] "pfmdr1 184Y/F"
#  [9,] "pfmdr1 YYXXY" : These are just giving EH info but in all studies prevalence can be identified from the other markers : remove
#  [10,] "pfmdr1 YYY" : Likewise with YYY and NFD
#  [11,] "pfmdr1 NFD"

# filter out unneeded EH info and now can easily assign locus groups
pdmdr1 <- pdmdr1 %>%
  filter(grepl("184|86|copy", mut)) %>%
  mutate(locus = NA) %>%
  mutate(locus = replace(locus, grepl("86", mut), "86")) %>%
  mutate(locus = replace(locus, grepl("184", mut), "184")) %>%
  mutate(locus = replace(locus, grepl("copy", mut), "CNV"))

res_loc <- c("pfmdr1 184F", "pfmdr1 86Y", "pfmdr1 copy number >1")
mix_loc <- c("pfmdr1 184Y/F", "pfmdr1 86N/Y")
wt_loc <- c("pfmdr1 N86", "pfmdr1 Y184")

# have they all been grouped - Yes
pdmdr1$locus %>% table(useNA = "a")

## STEP 2: Figure out how mixed infections work ------------------------------

# Some studies the mixed infections (184Y/F and 86N/Y) are reported separately
# in prevalence, i.e. prev of N86, 86Y and 86N/Y > 1

# put into study time lat points and work out which have different mixed numbers
pdmdr1 <- pdmdr1 %>%
  group_by(across(c(-x, -n, -prev, -mut, -mix, -rowid))) %>%
  mutate(newid = cur_group_id()) %>%
  mutate(non_mix = all(mix == x))

# non mix, i.e. where the mix counts equal to the x, all have the same n
pdmdr1 %>% filter(!non_mix) %>%
  group_by(newid) %>% summarise(n = length(unique(n))) %>% pull(n) %>% all


### TYPE 1 ------------------------------
# this means that these ones we can use x for x but must account for mix
# i.e. grab marker values from x for 86Y 184F and CNV
# and we then will need to add half from the mixed genotype

# for viewing each group
# pdmdr1 %>% filter(!non_mix) %>%
#   group_by(newid) %>%
#   mutate(xn = all(sum(x) == n[1])) %>%
#   filter(xn) %>%
#   split(.$newid)

mdrsplit1 <- pdmdr1 %>% filter(!non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  group_by(across(c(-x, -n, -prev, -mut, -mix, -rowid))) %>%
  summarise(x = x[mut %in% res_loc] + ifelse(any(mut %in% mix_loc), 0.5*x[mut %in% mix_loc], 0),
            n = unique(n),
            prev = x/n) %>%
  mutate(mut = "mdr1_86Y") %>%
  mutate(mut = replace(mut, locus == "184", "mdr1_184F"))

# Note: There used to be a type 2 here where the sum of x didn't equal n
# here but these are all due to typos at WWARN or in studies themselves
# All the above corrections earlier were by iterating through this type
### TYPE 2 ------------------------------
# pdmdr1 %>% filter(!non_mix) %>%
#   group_by(newid) %>%
#   mutate(xn = all(sum(x) == n[1])) %>%
#   filter(!xn) %>%
#   group_by(across(c(-x, -n, -prev, -mut, -mix, -rowid)))


# and for those where mix is the same and all x add up to n
# i.e. grab marker values from x for 86Y 184F and CNV
# checked these work with:
### TYPE 3 ------------------------------
pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  group_by(newid) %>% mutate(prev = x/n) %>% filter(grepl("184",mut)) %>% mutate(p = sum(prev)) %>% filter(p!=1)

pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  group_by(newid) %>% mutate(prev = x/n) %>% filter(grepl("86",mut)) %>% mutate(p = sum(prev)) %>% filter(p!=1)

# CNV the WT is not reported so this would not be expected to add to 1

# so just grab 86Y 184F and CNV rows from this filter
mdrsplit2 <- pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(xn) %>%
  group_by(newid) %>%
  mutate(prev = x/n) %>%
  filter(mut %in% res_loc)

### TYPE 4 ------------------------------

# This last group is now composed of samples where the
# x does not sum to n but the mix equals the n,
# Many of these groups only report the resistance prevalence
# so where there is only 1 record per newid and it is a resistance
# locus we can use those directly
mdrsplit3 <- pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(!xn) %>%
  filter(n()==1) %>%
  group_by(newid) %>%
  mutate(prev = x/n) %>%
  filter(mut %in% res_loc)

### TYPE 5 ------------------------------
# Some entries are for wildtype rather than the mutant
# so we need to change these round to reflect the mut prevalence
mdrsplit4 <- pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(!xn) %>%
  filter(n()==1) %>%
  group_by(newid) %>%
  filter(mut %in% wt_loc) %>%
  mutate(mut = replace(mut, mut == "pfmdr1 N86", "pfmdr1 86Y")) %>%
  mutate(mut = replace(mut, mut == "pfmdr1 Y184", "pfmdr1 184F")) %>%
  mutate(x = n-x) %>%
  mutate(prev = x/n)

### TYPE 6 ------------------------------
# The last type are multiple samples of CNV in the same
# locations/study etc but they are split out
# so group and sum these
mdrsplit5 <- pdmdr1 %>% filter(non_mix) %>%
  group_by(newid) %>%
  mutate(xn = all(sum(x) == n[1])) %>%
  filter(!xn) %>%
  filter(n()>1) %>%
  group_by(across(c(-x, -n, -prev, -mut,-mix,-rowid))) %>%
  summarise(x = sum(x), n = unique(n)) %>%
  mutate(mut = "mdr1_CNV") %>%
  mutate(prev = x/n)

## STEP 3. Bring all together againb-------------------

# and group by to record prevalence of each mdr1 marker type
mdr1ww_final_res_df <- rbind(mdrsplit1, mdrsplit2, mdrsplit3, mdrsplit4, mdrsplit5) %>%
  ungroup %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup %>%
  mutate(mut = replace(mut, mut == "pfmdr1 86Y", "mdr1_86Y")) %>%
  mutate(mut = replace(mut, mut == "pfmdr1 184F", "mdr1_184F")) %>%
  mutate(mut = replace(mut, mut == "pfmdr1 copy number >1", "mdr1_CNV")) %>%
  dplyr::mutate(gene_mut = str_replace(mut, "_","-")) %>%
  dplyr::left_join(select(validated, c("gene_mut", "annotation")), by = "gene_mut") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(gene_mut = clean_mutations(gene_mut)) %>%
  dplyr::relocate(names(k13ww_final_res_df))

# and the sanity check
(rbind(mdrsplit1, mdrsplit2, mdrsplit3, mdrsplit4, mdrsplit5) %>%
    pull(uuid) %>% length()) ==
  (pdmdr1$newid %>% unique() %>% length())

# ---------------------------------------------------- o
# 6. Sort PFPM23 ----
# ---------------------------------------------------- o

# because this is true we can ignore mix
all(pdwwspl$pfpm23$mix == pdwwspl$pfpm23$x)

pfpm23res <- pdwwspl$pfpm23 %>%
  select(-mix) %>%
  mutate(prev = x/n) %>%
  mutate(mut = replace(mut, mut == "pfpm2 copy number=1", "WT")) %>%
  mutate(mut = replace(mut, mut == "pfpm2 copy number >1", "pm23_CNV"))

# however because these aren't equal we now that the database has not been
# reporting prevalence
pfpm23res$mut %>% table

# add in the prev catch
pfpm23res <- pfpm23res %>%
  group_by(across(c(-x,-prev, -mut, -rowid))) %>%
  mutate(uuid = cur_group_id()) %>%
  ungroup() %>%
  group_by(uuid) %>%
  mutate(p = sum(prev))

# okay so all the non p == 1 appear to just have one row
pfpm23res %>% filter(p != 1) %>%
  split(.$uuid) %>% lapply(nrow) %>% unlist %>% table

# As a result it's actually easy as we can just filter to the CNV
# entries. WWARN has for some records reported just the CNV and for others
# has repoted both the CNV and the WT. But because the WT and CNV entries
# all sum to 1 then we can just filter to CNV
pfpm23ww_final_res_df <- pfpm23res %>% filter(mut == "pm23_CNV") %>%
  ungroup %>%
  select(iso3c, admin_0, admin_1, site, lat, long,
         year, study_start_year, study_end_year,
         x, n, prev, gene, mut, database, pmid, url, source) %>%
  ungroup

# ---------------------------------- o
# LAST -----------------------------
# ---------------------------------- o

# bring it all back together
wwarn_res_df <- rbind(crtww_final_res_df, mdr1ww_final_res_df, k13ww_final_res_df)
saveRDS(wwarn_res_df, here::here("analysis/data-derived/wwarn_res.rds"))
