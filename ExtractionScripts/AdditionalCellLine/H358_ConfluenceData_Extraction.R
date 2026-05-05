library(tidyverse)
library(dplyr)
library(here)
library(readxl)

##

## Task 0 - Load the Data

current_path <- here()

parent_path <- dirname(dirname(current_path))

new_path <- file.path(parent_path, "RawData", "AdditionalCellLine")

data_path <- file.path(parent_path, "ExtractedData", "AdditionalCellLine")

h358_data = read_excel(file.path(new_path,"PR_20240807_H358D2_sotorasib_adaptation_confluence.xlsx"))

#remove the rows with metadata
h358_data_filtered = h358_data[-(1:6), ]

#set the column names as the values of the first row
colnames(h358_data_filtered) = h358_data_filtered[1, ]

#remove the first row, since it is no longer needed
h358_data_filtered = h358_data_filtered[-1, ]

#remove the date/time column, since we only need elapsed time
h358_data_filtered = h358_data_filtered[ , -1]

#rename the elapsed column as time
colnames(h358_data_filtered)[1] <- "Time"

#pivot the data frame to a longer format
h358_data_long = h358_data_filtered %>% pivot_longer(!Time, names_to = "PlateIndex", values_to = "Confluency")

#save the data frame as an rds file for easier loading later
saveRDS(h358_data_long, file.path(data_path, "h358_sotorasib_adaptation_confluence.rds")) 
