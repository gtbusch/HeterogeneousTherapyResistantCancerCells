## This is a script to take the generated TPM and metadata files and combine them

library(ggplot2)
library(tidyverse)
library(here)

## 

# Task 0: Load the appropriate data files

#get the path to the parent folder
current_path <- here()

parent_path <- dirname(dirname(current_path))

new_path <- file.path(parent_path, "RawData", "BulkRNASeq_Run1")

data_path <- file.path(parent_path, "ExtractedData", "BulkRNASeq_Run1")

#list all subfolders in the parent folder
subfolders <- list.dirs(new_path, full.names = TRUE, recursive = TRUE)

#initialize a vector to store the paths of library RDS files
tpm_rds_files <- c()

#loop through each subfolder and find the library plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "^GTB_BulkRNASeq_Run1_Step1_TxImport_(metadata|myTPM|myCounts)\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  tpm_rds_files <- c(tpm_rds_files, files_in_subfolder)
}

for (file_path in tpm_rds_files) {
  # Determine which file we're dealing with based on the filename
  if (grepl("metadata", basename(file_path), ignore.case = TRUE)) {
    var_name <- "metadata"
  } else if (grepl("myTPM", basename(file_path), ignore.case = TRUE)) {
    var_name <- "myTPM"
  } else if (grepl("myCounts", basename(file_path), ignore.case = TRUE)) {
    var_name <- "myCounts"
  } else {
    warning(paste("Unexpected file:", basename(file_path)))
    next  # Skip this file and continue with the next iteration
  }
  
  # Read the RDS file and assign it to the determined variable name
  assign(var_name, readRDS(file_path))
}

## 

# Task 1: Manipulate the data frames to reflect the metadata

myTPM = myTPM %>% select(-gene_name)

myCounts = myCounts %>% select(-gene_name)

task1_TPM <- myTPM %>%
  # First, gather all columns except target_id
  pivot_longer(cols = -target_id, names_to = "sample", values_to = "expression") %>%
  # Then spread the data so genes become columns
  pivot_wider(names_from = target_id, values_from = expression) %>%
  # Finally, convert the sample column to row names
  column_to_rownames("sample")

task1_Counts <- myCounts %>%
  # First, gather all columns except target_id
  pivot_longer(cols = -target_id, names_to = "sample", values_to = "expression") %>%
  # Then spread the data so genes become columns
  pivot_wider(names_from = target_id, values_from = expression) %>%
  # Finally, convert the sample column to row names
  column_to_rownames("sample")

task1_annotatedTPM <- cbind(metadata, task1_TPM)

task1_annotatedCounts <- cbind(metadata, task1_Counts)

saveRDS(task1_annotatedTPM, file = file.path(data_path, "bulkrnaseq_run1_tpm_plusmetadata.rds"))

saveRDS(task1_annotatedCounts, file = file.path(data_path, "bulkrnaseq_run1_counts_plusmetadata.rds"))

write_tsv(task1_annotatedTPM, file.path(data_path, "bulkrnaseq_run1_tpm_plusmetadata.tsv"))

write_tsv(task1_annotatedCounts, file.path(data_path, "bulkrnaseq_run1_counts_plusmetadata.tsv"))

##