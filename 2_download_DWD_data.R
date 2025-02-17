# Get climate data from the DWD: Deutsche weather institute:


# STRAT --------------
library(rvest)
library(stringr)
library(fs)  # For directory creation

# Base URLs for different climate variables
base_urls <- list(
  temperature = "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/",
  precipitation = "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation/"
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
  
  # Extract filenames and filter for year 2000
  filenames <- stringr::str_extract(link_text[-1], "\\d+")
  year_2000_files <- grep("^2000", filenames, value = TRUE)
  year_2000_links <- links[filenames %in% year_2000_files]
  
  if (length(year_2000_links) == 0) {
    message(paste("No files found for", variable, month, "in the year 2000"))
    return(NULL)
  }
  
  # Create directories: 'rawData/variable/month'
  save_dir <- paste0("rawData/", variable, "/", month)
  dir_create(save_dir, recurse = TRUE)
  
  # Set download timeout
  options(timeout = max(600, getOption("timeout")))
  
  # Download files
  for (i in seq_along(year_2000_links)) {
    file_path <- paste0(save_dir, "/", year_2000_files[i], "asc.gz")
    if (!file_exists(file_path)) {
      download.file(year_2000_links[i], file_path, mode = "wb")
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


# --- END 

































library(rvest)
library(stringr)

page_link <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/"
page_link_precip <- 'https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation/'

month_name <- c("01_Jan",
                "02_Feb",
                "03_Mar",
                "04_Apr",
                "05_May",
                "06_Jun",
                "07_Jul",
                "08_Aug",
                "09_Sep",
                "10_Oct",
                "11_Nov",
                "12_Dec")

# you can set month_name as "05_May" to get the data from 05_May
download_data <- function(x, ...) {
  print(x)
  
  month_name = x
  print(month_name)
  
  # getting the html page for 01_Jan folder
  page <- read_html(paste0(page_link_precip, month_name, "/"))
  
  # getting the link text
  link_text <- page %>%
    html_elements("a") %>%
    html_text()
  
  # creating links
  links <- paste0(page_link_precip, month_name, "/", link_text)[-1]
  
  # extracting the numbers for filename
  filenames <- stringr::str_extract(pattern = "\\d+", string = link_text[-1])
  
  # creating a directory
  dir.create(month_name)
  
  # setting the option for maximizing time limits for downloading
  options(timeout = max(600, getOption("timeout")))
  
  # downloading the file
  for (i in seq_along(links)) {
    download.file(links[i], paste0(month_name, "/", filenames[i], "asc.gz"))
  }
  
}


lapply(month_name, download_data)
