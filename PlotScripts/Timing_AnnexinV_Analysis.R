library(tidyverse)
library(dplyr)
library(here)
library(readxl)
library(ggplot2)
library(ggpubr)

##

# This is a script to analyze Annexin V and Confluency data obtained via IncuCyte software

# Created by GTB on 10/9/2025

# Last edited by GTB on 11/3/2025

##

## Task 0 - Load the Data ----

# Initialize an empty master dataframe to store all combined data
timing_experiment_data <- data.frame()

current_path <- here()
parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "TimingValidation")

figure_path <- file.path(parent_path, "Plots")

# Get a list of all subfolders in the current directory
subfolders <- list.dirs(path = new_path, full.names = TRUE, recursive = FALSE)

# Loop through each subfolder
for (folder in subfolders) {
  # Find the confluence, annexin, and metadata files in the current subfolder
  confluence_file <- list.files(path = folder, pattern = "Confluence\\.xlsx$", full.names = TRUE)
  annexin_file <- list.files(path = folder, pattern = "Annexin\\.xlsx$", full.names = TRUE)
  metadata_file <- list.files(path = folder, pattern = "Metadata\\.csv$", full.names = TRUE)
  
  if (length(confluence_file) > 0 && length(metadata_file) > 0) {
    # Read the confluence data
    confluence_data <- read_excel(confluence_file)
    confluence_data <- confluence_data[-(1:6), -(1:2)]
    confluence_data <- t(confluence_data)
    colnames(confluence_data) <- c("XY", "PhaseConfluence")
    rownames(confluence_data) <- 1:length(confluence_data[,1])
    confluence_data  <- as.data.frame(confluence_data)
    
    #Read the annexin data
    annexin_data <- read_excel(annexin_file)
    annexin_data <- annexin_data[-(1:6), -(1:2)]
    annexin_data <- t(annexin_data)
    colnames(annexin_data) <- c("XY", "AnnexinConfluence")
    rownames(annexin_data) <- 1:length(annexin_data[,1])
    annexin_data <- as.data.frame(annexin_data)
    
    # Read the metadata
    metadata <- read.csv(metadata_file)
    
    # Add metadata columns
    experiment_data <- left_join(confluence_data, metadata, by = "XY")
    experiment_data <- left_join(experiment_data, annexin_data, by = "XY")
    
    # Append to master dataframe
    timing_experiment_data <- rbind(timing_experiment_data, experiment_data)
    
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

condition_labels <- c(
  "Vemurafenib Only" = "Vemurafenib for 5 Weeks",       
  "21 Days Vemurafenib" = "Vem. to Combination at Week 3", 
  "Combination Treatment" = "Combination for 5 Weeks",  
  "7 Days Vemurafenib" = "Vem. to Combination at Week 1",       
  "14 Days Vemurafenib" = "Vem. to Combination at Week 2",
  "DMSO Only" = "Negative Control",
  "Saracatinib Only" = "Saracatinib for 5 Weeks",
  "SR18662" = "SR18662 for 5 Weeks"
)

##

## Task 1 - Plot Annexin V area vs confluency for Saracatinib ----

saracatinib_data = timing_experiment_data %>% filter((Plate == 1) | (Plate == 2))

# Reshape the data from wide to long format
singleagent_long <- saracatinib_data %>%
  pivot_longer(
    cols = c(PhaseConfluence, AnnexinConfluence),
    names_to = "ConfluenceType",
    values_to = "Confluence"
  ) %>%
  mutate(
    Confluence = as.numeric(Confluence),
    ConfluenceType = gsub("Confluence", "", ConfluenceType)
  )

# Calculate means for each group
singleagent_summary <- singleagent_long %>%
  group_by(Condition, ConfluenceType) %>%
  summarise(
    mean_confluence = mean(Confluence, na.rm = TRUE),
    sd_confluence = sd(Confluence, na.rm = TRUE),
    n = n(),
    se_confluence = sd_confluence / sqrt(n),
    .groups = "drop"
  )

# Create the grouped bar plot
ggplot((singleagent_summary %>% mutate(Condition = factor(Condition, levels = c("Combination", 
                                                                                "7 Days Vemurafenib", 
                                                                                "14 Days Vemurafenib",
                                                                                "21 Days Vemurafenib",
                                                                                "Vemurafenib Only",
                                                                                "Saracatinib Only",
                                                                                "DMSO Only")))), 
       aes(x = Condition, y = mean_confluence, 
                            fill = ConfluenceType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_confluence - se_confluence, 
                    ymax = mean_confluence + se_confluence),
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.5) +
  labs(
    x = "Condition",
    y = "Confluence (%)",
    fill = "Confluence Type",
    title = "Timing of 1 µM Saracatinib and Vemurafenib: Phase and Annexin V Confluence by Condition"
  ) +
  scale_fill_manual(values = c("Phase" = "#7f8c8d", "Annexin" = "#c0392b")) +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.position = c(0.27, 0.85),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    strip.text = element_text(size = 25.6, face = "bold")
  )

ggsave(file.path(figure_path,"saracatinibtiming_annexinandconfluence.pdf"), width = 24, height = 7, units = "in")

## Task 2 - Plot Annexin V area vs confluency for SR18662 ----

sr18662_data = timing_experiment_data %>% filter((Plate == 3) | (Plate == 4))

# Reshape the data from wide to long format
singleagent_long <- sr18662_data %>%
  pivot_longer(
    cols = c(PhaseConfluence, AnnexinConfluence),
    names_to = "ConfluenceType",
    values_to = "Confluence"
  ) %>%
  mutate(
    Confluence = as.numeric(Confluence),
    ConfluenceType = gsub("Confluence", "", ConfluenceType)
  )

# Calculate means for each group
singleagent_summary <- singleagent_long %>%
  group_by(Condition, ConfluenceType) %>%
  summarise(
    mean_confluence = mean(Confluence, na.rm = TRUE),
    sd_confluence = sd(Confluence, na.rm = TRUE),
    n = n(),
    se_confluence = sd_confluence / sqrt(n),
    .groups = "drop"
  )

# Create the grouped bar plot
ggplot((singleagent_summary %>% mutate(Condition = factor(Condition, levels = c("Combination", 
                                                                                "7 Days Vemurafenib", 
                                                                                "14 Days Vemurafenib",
                                                                                "21 Days Vemurafenib",
                                                                                "Vemurafenib Only",
                                                                                "SR18662 Only",
                                                                                "DMSO Only")))), 
       aes(x = Condition, y = mean_confluence, 
           fill = ConfluenceType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_confluence - se_confluence, 
                    ymax = mean_confluence + se_confluence),
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.5) +
  labs(
    x = "Condition",
    y = "Confluence (%)",
    fill = "Confluence Type",
    title = "Timing of 1 µM SR18662 and Vemurafenib: Phase and Annexin V Confluence by Condition"
  ) +
  scale_fill_manual(values = c("Phase" = "#7f8c8d", "Annexin" = "#c0392b")) +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.position = c(0.27, 0.85),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    strip.text = element_text(size = 25.6, face = "bold")
  )

ggsave(file.path(figure_path,"sr18662timing_annexinandconfluence.pdf"), width = 24, height = 7, units = "in")
