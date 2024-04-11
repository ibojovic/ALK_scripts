# install.packages(c("xml2", "dplyr", "purrr"))

library(xml2)
library(dplyr)
library(purrr)

# Function to parse a single XML file and return a data frame
parse_xml_to_df <- function(file_path) {
  xml_data <- read_xml(file_path)
  all_elements <- xml_find_all(xml_data, "//*")
  
  element_names <- list()
  element_texts <- list()
  
  for (i in seq_along(all_elements)) {
    element_names[[i]] <- xml_name(all_elements[[i]])
    element_texts[[i]] <- xml_text(all_elements[[i]])
  }
  
  df <- data.frame(
    FileName = basename(file_path),
    ElementName = unlist(element_names),
    ElementText = unlist(element_texts),
    stringsAsFactors = FALSE
  )
  
  return(df)
}

# List all XML files recursively in the target directory
xml_files <- list.files("~/downloads/TCGA/LUAD/clinical_new/", pattern = "\\.xml$", recursive = TRUE, full.names = TRUE)

# Parse each XML file and combine into one data frame
combined_df <- map_dfr(xml_files, parse_xml_to_df)

# Show the combined data frame
print(combined_df)
