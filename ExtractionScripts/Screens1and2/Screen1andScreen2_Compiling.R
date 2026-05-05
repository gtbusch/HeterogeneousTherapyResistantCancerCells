library(tidyverse)
library(dplyr)
library(here)

##

##Task a - Load the Data

#get the path to the parent folder
current_path <- here()

parent_path <- dirname(dirname(current_path))

new_path <- file.path(parent_path, "RawData", "HighThroughputScreens")

data_path <- file.path(parent_path, "ExtractedData", "HighThroughputScreens")

#list all subfolders in the parent folder
subfolders <- list.dirs(new_path, full.names = TRUE, recursive = TRUE)

#initialize a vector to store the paths of library RDS files
library_rds_files <- c()

#loop through each subfolder and find the library plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "^screen[12]_library.*\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  library_rds_files <- c(library_rds_files, files_in_subfolder)
}

for (file_path in library_rds_files) {
  # Extract the file name without the extension to use as the variable name
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read the RDS file and assign it to a dynamically named variable
  assign(file_name, readRDS(file_path))
}

#initialize a vector to store the paths of control RDS files
control_rds_files <- c()

#loop through each subfolder and find the control plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "^screen[12]_control.*\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  control_rds_files <- c(control_rds_files, files_in_subfolder)
}

for (file_path in control_rds_files) {
  # Extract the file name without the extension to use as the variable name
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Read the RDS file and assign it to a dynamically named variable
  assign(file_name, readRDS(file_path))
}

##

## Task 2 - Process the data to combine 

screen1_librarydata_long$Screen <- rep(c(1), times = length(screen1_librarydata_long$Name))
screen2_librarydata_long$Screen <- rep(c(2), times = length(screen2_librarydata_long$Name))

screen1_librarydata_long = screen1_librarydata_long %>% filter(Name %in% screen2_librarydata_long$Name)
panscreen_long = rbind(screen1_librarydata_long, screen2_librarydata_long)

screen1_librarydata$Screen <- rep(c(1), times = length(screen1_librarydata$Name))
screen2_librarydata$Screen <- rep(c(2), times = length(screen2_librarydata$Name))

screen1_librarydata = screen1_librarydata %>% filter(Name %in% screen2_librarydata$Name)
panscreen = rbind(screen1_librarydata, screen2_librarydata)

saveRDS(panscreen,file.path(data_path, "panscreen.rds"))
saveRDS(panscreen_long, file.path(data_path, "panscreen_long.rds"))

##

## Task 3 - Process the control data to combine

screen1_controldata$Screen <- rep(c(1), times = length(screen1_controldata$Drug))
screen2_controldata$Screen <- rep(c(2), times = length(screen2_controldata$Drug))

screen2_controldata <- screen2_controldata %>%
  mutate(Drug = case_when(
    Drug == "Clofarabine" ~ "Dasatinib",  # Change Clofarabine to Dasatinib
    Drug == "Dasatinib" ~ "Clofarabine",  # Change Dasatinib to Clofarabine
    TRUE ~ Drug                           # Leave all other drugs unchanged
  ))

panscreen_controls_long = rbind(screen1_controldata, screen2_controldata)
saveRDS(panscreen_controls_long, file.path(data_path, "panscreen_controls_long.rds"))

##
