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

new_path <- file.path(parent_path, "ExtractedData", "HTSValidation", "AnnexinVData")

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

custom_names <- c("Parental" = "Drug-Naive", "E2" = "Resistant Clone 1", "D13" = "Resistant Clone 2", "D4" = "Resistant Clone 3", "C12" = "Resistant Clone 4", "C14" = "Resistant Clone 5",
                  "D5" = "Resistant Clone 6", "F9" = "Resistant Clone 7", "G12" = "Resistant Clone 8", "I11" = "Resistant Clone 9", "A2" = "Resistant Clone 10", "A11" = "Resistant Clone 11",
                  "C3" = "Resistant Clone 12", "E8" = "Resistant Clone 13")

##

## Task 1 - Plot Annexin V area vs confluency for SR18662 and Saracatinib ----

singleagent_data = combined_experiment_data %>% filter((Round == 1 & Plate == 1) | (Round == 1 & Plate == 2)) %>% filter(Condition %in% c("Bortezomib", "SR18662", "Saracatinib", "DMSO")) %>%
  filter(!(CellLine == "Parental" & Treatment == "Vemurafenib"))

# Reshape the data from wide to long format
singleagent_long <- singleagent_data %>%
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
  group_by(CellLine, Condition, ConfluenceType) %>%
  summarise(
    mean_confluence = mean(Confluence, na.rm = TRUE),
    sd_confluence = sd(Confluence, na.rm = TRUE),
    n = n(),
    se_confluence = sd_confluence / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Recode the condition labels
    Condition = case_when(
      Condition == "DMSO" ~ "Negative Control",
      Condition == "Bortezomib" ~ "Positive Control",
      TRUE ~ Condition
    ),
    # Set the factor order for conditions
    Condition = factor(Condition, levels = c("Positive Control", "SR18662", "Saracatinib", "Negative Control")),
    CellLine = factor(CellLine, levels = c("Parental", "E2", "D13", "D4", "C12"))
  )

# Create the grouped bar plot
ggplot(singleagent_summary, aes(x = Condition, y = mean_confluence, 
                            fill = ConfluenceType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_confluence - se_confluence, 
                    ymax = mean_confluence + se_confluence),
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.5) +
  facet_wrap(~ CellLine, ncol = 5, labeller = labeller(CellLine = custom_names)) +
  labs(
    x = "Condition",
    y = "Confluence (%)",
    fill = "Confluence Type",
    title = "Phase and Annexin V Confluence by Condition and Cell Line"
  ) +
  scale_fill_manual(values = c("Phase" = "#7f8c8d", "Annexin" = "#c0392b")) +
  custom_figure_theme() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.position = c(0.27, 0.85),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    strip.text = element_text(size = 25.6, face = "bold")
  )

ggsave(file.path(figure_path,"singleagent_annexinandconfluence.pdf"), width = 24, height = 7, units = "in")

## Task 2 - Plot Annexin V area vs confluency for Deguelin, SB273005, and the combo ----

combo1_data = combined_experiment_data %>% filter((Round == 2 & Plate == 1) | (Round == 2 & Plate == 2)) %>% filter(Condition %in% c("Bortezomib", "Deguelin", "SB273005", "Deguelin + SB273005", "DMSO")) %>%
  filter(!(CellLine == "Parental" & Treatment == "Vemurafenib"))

# Reshape the data from wide to long format
combo1_long <- combo1_data %>%
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
combo1_summary <- combo1_long %>%
  group_by(CellLine, Condition, ConfluenceType) %>%
  summarise(
    mean_confluence = mean(Confluence, na.rm = TRUE),
    sd_confluence = sd(Confluence, na.rm = TRUE),
    n = n(),
    se_confluence = sd_confluence / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Recode the condition labels
    Condition = case_when(
      Condition == "DMSO" ~ "Negative Control",
      Condition == "Bortezomib" ~ "Positive Control",
      TRUE ~ Condition
    ),
    # Set the factor order for conditions
    Condition = factor(Condition, levels = c("Positive Control", "Deguelin", "SB273005", "Deguelin + SB273005", "Negative Control")),
    CellLine = factor(CellLine, levels = c("Parental", "E2", "D13", "D4", "C12"))
  )

# Create the grouped bar plot
ggplot(combo1_summary, aes(x = Condition, y = mean_confluence, 
                                fill = ConfluenceType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_confluence - se_confluence, 
                    ymax = mean_confluence + se_confluence),
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.5) +
  facet_wrap(~ CellLine, ncol = 5, labeller = labeller(CellLine = custom_names)) +
  labs(
    x = "Condition",
    y = "Confluence (%)",
    fill = "Confluence Type",
    title = "Phase and Annexin V Confluence by Condition and Cell Line"
  ) +
  scale_fill_manual(values = c("Phase" = "#7f8c8d", "Annexin" = "#c0392b")) +
  custom_figure_theme() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.position = c(0.07, 0.70),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(size = 25.6, face = "bold")
  )

ggsave(file.path(figure_path,"combo1_annexinandconfluence.pdf"), width = 24, height = 5.5, units = "in")




## Task 3 - Plot Annexin V area vs confluency for IWP-O1, SB273005, and the combo ----

combo2_data = combined_experiment_data %>% filter((Round == 2 & Plate == 1) | (Round == 2 & Plate == 2)) %>% filter(Condition %in% c("Bortezomib", "IWP-O1", "SB273005", "IWP-O1 + SB273005", "DMSO")) %>%
  filter(!(CellLine == "Parental" & Treatment == "Vemurafenib"))

# Reshape the data from wide to long format
combo2_long <- combo2_data %>%
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
combo2_summary <- combo2_long %>%
  group_by(CellLine, Condition, ConfluenceType) %>%
  summarise(
    mean_confluence = mean(Confluence, na.rm = TRUE),
    sd_confluence = sd(Confluence, na.rm = TRUE),
    n = n(),
    se_confluence = sd_confluence / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Recode the condition labels
    Condition = case_when(
      Condition == "DMSO" ~ "Negative Control",
      Condition == "Bortezomib" ~ "Positive Control",
      TRUE ~ Condition
    ),
    # Set the factor order for conditions
    Condition = factor(Condition, levels = c("Positive Control", "IWP-O1", "SB273005", "IWP-O1 + SB273005", "Negative Control")),
    CellLine = factor(CellLine, levels = c("Parental", "E2", "D13", "D4", "C12"))
  )

# Create the grouped bar plot
ggplot(combo2_summary, aes(x = Condition, y = mean_confluence, 
                           fill = ConfluenceType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), 
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_confluence - se_confluence, 
                    ymax = mean_confluence + se_confluence),
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 0.5) +
  facet_wrap(~ CellLine, ncol = 5, labeller = labeller(CellLine = custom_names)) +
  labs(
    x = "Condition",
    y = "Confluence (%)",
    fill = "Confluence Type",
    title = "Phase and Annexin V Confluence by Condition and Cell Line"
  ) +
  scale_fill_manual(values = c("Phase" = "#7f8c8d", "Annexin" = "#c0392b")) +
  custom_figure_theme() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    legend.background = element_blank(),
    legend.position = c(0.07, 0.70),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(size = 25.6, face = "bold")
  )

ggsave(file.path(figure_path,"combo2_annexinandconfluence.pdf"), width = 24, height = 5.5, units = "in")
