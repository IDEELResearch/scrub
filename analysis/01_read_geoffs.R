# now convert all of the geoff objects into a list
geoff_path <- "analysis/data-geoff/" # filepath to all of the geoff tsvs which are in separate folders
folders <- list.dirs(path = geoff_path)

geoff <- data.frame()

for(i in 2:length(folders)) { # the first path is the root which we want to omit
  files <- list.files(path = folders[i],
                      all.files = TRUE, full.names = TRUE) %>%
    str_subset(".tsv")
  # relies on a standardised naming convention of the .tsv files
  # TODO: update this to reflect the real unique string
  study_path <- str_subset(files, "study")
  site_path <- str_subset(files, "site")
  prev_path <- str_subset(files, "prevalence")
  
  geoff_obj <- read_geoff(study_path = study_path,
                          site_path = site_path,
                          prev_path = prev_path)
  
  # this is now a list of data frames ready for STAVE
  geoff <- rbind(geoff, geoff_obj)
  
}
