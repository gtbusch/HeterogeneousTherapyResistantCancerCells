library(tidyverse)
library(here)
library(plotly)
library(ggrepel)
library(circlize)
library(matrixStats)
library(data.table)
library(openxlsx)
library(gridExtra)
library(readxl)
library(purrr)

##

# This is a script to compile and further process the extracted data generated from the raw data excel files for Follow Up Panel A from Wistar HTS. 

# Created by GTB on 1/12/2024

# Last edited 6/27/2024

##

## Task 0: Load the library plate RDS data and control RDS data from each subfolder

#get the path to the parent folder
current_path <- here()

parent_path <- dirname(dirname(current_path))

new_path <- file.path(parent_path, "RawData", "FollowUpPanelA")

data_path <- file.path(parent_path, "ExtractedData", "FollowUpPanelA")

#list all subfolders in the parent folder
subfolders <- list.dirs(new_path, full.names = TRUE, recursive = TRUE)

#initialize a vector to store the paths of library RDS files
library_rds_files <- c()

#loop through each subfolder and find the library plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "^followup.*\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  library_rds_files <- c(library_rds_files, files_in_subfolder)
}

#check if any files were found
print(library_rds_files)

#load and combine the RDS files
followup_panela <- map(library_rds_files, readRDS) %>%
  bind_rows()

#initialize a vector to store the paths of control RDS files
control_rds_files <- c()

#loop through each subfolder and find the control plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "^controls_followup.*\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  control_rds_files <- c(control_rds_files, files_in_subfolder)
}

#check if any files were found
print(control_rds_files)

#load and combine the RDS files
controls_followup_panela <- map(control_rds_files, readRDS) %>%
  bind_rows()

##

## Task 1: Adjust the format of the library data to match analysis code and clean up doses

#mutate the data to long form by turning the 16 toxicity columns into Toxicity, Dose, and Replicate Columns
followup_panela_long <- followup_panela %>%
  pivot_longer(
    cols = starts_with("Dose"),
    names_to = c("Dose", "Replicate"),
    names_pattern = "Dose_([A-H])_(\\d)",
    values_to = "Toxicity"
  ) %>%
  mutate(
    Dose = factor(Dose, levels = c("A", "B", "C", "D", "E", "F", "G", "H")),
    DoseIndex = as.numeric(Dose)  # Convert the Dose letter to a numeric index
  )

#use the dose letter index and the starting dose to calculate the actual dose for each row
followup_panela_long <- followup_panela_long %>%
  rowwise() %>%  # Ensure the operation is performed row by row
  mutate(
    StartingDose = as.numeric(StartingDose),
    # Calculate the actual dose based on the dilution series
    ActualDose = StartingDose / (3.16 ^ (DoseIndex - 1))
  ) %>%
  select(-DoseIndex) %>%  # Remove the DoseIndex column
  ungroup()

#adjust the format of the numbers in the Actual Dose Column and remove dose information columns
followup_panela_long = followup_panela_long %>%
  mutate(ActualDose = formatC(ActualDose, format = "f", digits = 6, flag = "-"),
         Dose = NULL) %>%  # Remove the original Dose column)
         rename(Dose = ActualDose) %>% #rename the actual dose column
         select(-StartingDose) #remove the starting dose column

##

## Task 2 - Save the new follow up panel long variables

saveRDS(followup_panela_long, file.path(data_path, "followup_panela_long.rds"))

saveRDS(controls_followup_panela, file.path(data_path, "controls_followup_panela_long.rds"))

##






