library(tidyverse)
library(dplyr)
library(here)
library(readxl)
library(ggplot2)
library(ggpubr)

##

# This is a script to analyze Annexin V and Confluency data obtained via IncuCyte software

# Created by GTB on 10/9/2025

# Last edited by GTB on 10/9/2025

##

## Task 0 - Load the Data ----

# Initialize an empty master dataframe to store all combined data
combined_experiment_data <- data.frame()

current_path <- here()
parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "RCGrowthCurve")

figure_path <- file.path(parent_path, "Plots")

# Get a list of all subfolders in the current directory
subfolders <- list.dirs(path = new_path, full.names = TRUE, recursive = FALSE)

# Loop through each subfolder
for (folder in subfolders) {
  # Find the confluence, annexin, and metadata files in the current subfolder
  confluence_file <- list.files(path = folder, pattern = "Confluence\\.xlsx$", full.names = TRUE)
  metadata_file <- list.files(path = folder, pattern = "Metadata\\.csv$", full.names = TRUE)
  
  if (length(confluence_file) > 0 && length(metadata_file) > 0) {
    # Read the confluence data
    confluence_data <- read_excel(confluence_file)
    confluence_data <- confluence_data[-(1:6), -(1:2)]
    confluence_data <- t(confluence_data)
    colnames(confluence_data) <- c("XY", "Confluence_Day0", "ConfluenceDay7")
    rownames(confluence_data) <- 1:length(confluence_data[,1])
    confluence_data  <- as.data.frame(confluence_data)
    
    # Read the metadata
    metadata <- read.csv(metadata_file)
    
    # Add metadata columns
    experiment_data <- left_join(confluence_data, metadata, by = "XY")
    
    # Append to master dataframe
    combined_experiment_data <- rbind(combined_experiment_data, experiment_data)
    
    cat("Processed data from folder:", folder, "\n")
  } else {
    cat("Missing required files in folder:", folder, "\n")
  }
}

custom_figure_theme <- function() {
  theme(
    text = element_text(size = 25.6, family = "Helvetica", color = "black"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "white", linewidth = 1.2),
    panel.grid.minor = element_line(color = "white", linewidth = 0.6),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 19.2),
    axis.title = element_text(size = 25.6, face = "bold"),
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 25.6, hjust = 0.5),
    plot.caption = element_text(size = 19.2, face = "bold", hjust = 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 19.2),
    legend.title = element_text(size = 25.6, face = "bold")
  )
}

custom_names <- c("Parental" = "Drug-Naive", "E2" = "Resistant Clone 1", "D13" = "Resistant Clone 2", "D4" = "Resistant Clone 3", "C12" = "Resistant Clone 4")

custom_colors <- c("Parental" = "#5E5E5E", "D13" = "#FF7E79", "E2" = "#4294F8", "D4" = "#9E6D15", "C12" = "#662D91")

##

## Task 1 - Plot RC Growth Curve for Negative Controls ----

growthcurve_data = combined_experiment_data %>% filter(Condition == "DMSO") %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

growthcurve_data_foldchange = growthcurve_data %>% mutate(ConfluenceDay7 = as.numeric(ConfluenceDay7),
                                                       Confluence_Day0 = as.numeric(Confluence_Day0),
                                                       FoldChange = ConfluenceDay7/Confluence_Day0) %>% mutate(CellLine = factor(CellLine, levels = c("Parental", "E2", "D13", "D4", "C12")))

ggplot((growthcurve_data_foldchange %>%
          group_by(CellLine) %>%
          summarise(
            mean_FoldChange = mean(FoldChange, na.rm = TRUE),
            se_FoldChange = sd(FoldChange, na.rm = TRUE) / sqrt(n()),
            n = n()
          )), aes(x = CellLine, y = mean_FoldChange)) +
  geom_col(aes(fill = CellLine)) +
  scale_fill_manual(values = custom_colors) +
  geom_errorbar(aes(ymin = mean_FoldChange - se_FoldChange, 
                    ymax = mean_FoldChange + se_FoldChange),
                width = 0.2) +
  labs(
    title = "Average Confluence Change",
    x = "Cell Line",
    y = "Confluence Change (Day 0 to Day 7)"
  ) +
  scale_x_discrete(labels = custom_names) +
  custom_figure_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(figure_path,"rc_growthcurve.pdf"), width = 9.6, height = 9.6, units = "in")

