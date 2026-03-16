# clean RIS files

# Load library
library(revtools)

# function to clean RIS files
clean_ris <- function(file) {
  # 1. Read raw bytes
  raw_bytes <- readBin(file, what = "raw", n = file.info(file)$size)
  
  # 2. Remove all NUL bytes at the byte level
  raw_bytes <- raw_bytes[raw_bytes != as.raw(0)]
  
  # 3. Convert all remaining bytes into one UTF-8 string
  # use iconv to safely handle any invalid bytes
  full_text <- rawToChar(raw_bytes, multiple = FALSE)
  full_text <- iconv(full_text, from = "", to = "UTF-8", sub = "")
  
  # 4. Split into lines manually
  lines <- unlist(strsplit(full_text, "\r?\n"))
  
  # 5. Replace problematic characters
  lines <- gsub("\u2043", "-", lines)                    # hyphen bullet
  lines <- gsub("[\u2010\u2011\u2012\u2013\u2014]", "-", lines)  # other dash variants
  lines <- gsub("<[0-9a-fA-F]+>", "", lines)             # RIS control tags like <8f>
  
  # 6. Fix mojibake
  lines <- gsub("Ã©", "é", lines)
  lines <- gsub("Ã¨", "è", lines)
  lines <- gsub("Ã¢", "â", lines)
  lines <- gsub("Ãª", "ê", lines)
  lines <- gsub("Ã°", "ð", lines)
  lines <- gsub("Ã§", "ç", lines)
  lines <- gsub("Ã±", "ñ", lines)
  lines <- gsub("Ã´", "ô", lines)
  lines <- gsub("Ã¹", "ù", lines)
  
  # 7. Write cleaned RIS
  out <- sub("\\.ris$", "_clean.ris", file)
  writeLines(lines, out, useBytes = TRUE)
  
  return(out)
}




