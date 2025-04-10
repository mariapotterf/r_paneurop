# Get climate data from the DWD: Deutsche weather institute:


library(rvest)
library(stringr)
library(fs)  # For directory creation

# Base URLs for different climate variables
base_urls <- list(
  temp = "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/",
  precip = "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation/"
)

# List of month folder names
month_names <- c("01_Jan", "02_Feb", "03_Mar", "04_Apr", "05_May", "06_Jun",
                 "07_Jul", "08_Aug", "09_Sep", "10_Oct", "11_Nov", "12_Dec")

# Function to download data for a specific variable and month
download_data <- function(variable, month) {
  print(paste("Processing:", variable, "-", month))
  
  # Construct page URL
  page_url <- paste0(base_urls[[variable]], month, "/")
  
  # Read HTML page
  page <- tryCatch(read_html(page_url), error = function(e) return(NULL))
  if (is.null(page)) {
    message(paste("Failed to read:", page_url))
    return(NULL)
  }
  
  # Extract file links and names
  link_text <- page %>% html_elements("a") %>% html_text()
  links <- paste0(page_url, link_text)[-1]
  
  # Extract filenames and filter for years 1950-2024
  filenames <- stringr::str_extract(link_text[-1], "\\d+")
  year_filtered_files <- grep("^(195[0-9]|196[0-9]|197[0-9]|198[0-9]|199[0-9]|200[0-9]|201[0-9]|202[0-4])", filenames, value = TRUE)
  year_filtered_links <- links[filenames %in% year_filtered_files]
  
  if (length(year_filtered_links) == 0) {
    message(paste("No files found for", variable, month, "in the years 1950-2024"))
    return(NULL)
  }
  
  # Create directories: 'rawData/variable/month'
  save_dir <- paste0("rawData/", variable, "/", month)
  dir_create(save_dir, recurse = TRUE)
  
  # Set download timeout
  options(timeout = max(600, getOption("timeout")))
  
  # Download files
  for (i in seq_along(year_filtered_links)) {
    file_path <- paste0(save_dir, "/", year_filtered_files[i], "asc.gz")
    if (!file_exists(file_path)) {
      download.file(year_filtered_links[i], file_path, mode = "wb")
      message(paste("Downloaded:", file_path))
    } else {
      message(paste("File already exists:", file_path))
    }
  }
}

# Run the function for each variable and month
for (variable in names(base_urls)) {
  lapply(month_names, function(month) download_data(variable, month))
}
