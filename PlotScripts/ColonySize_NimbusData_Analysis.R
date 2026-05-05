library(tidyverse)
library(dplyr)
library(here)
library(readxl)
library(ggplot2)
library(ggpubr)

##

# This is a script to analyze colony count data obtained via NimbusImage. This script is designed to analyze all data taken for combination experiments

# Created by GTB on 4/4/2024

# Last edited by GTB on 10/20/2025

##

## Task 0 - Load the Data

# Initialize an empty master dataframe to store all combined data
combination_experiment_data <- data.frame()

all_metadata <- data.frame()

current_path <- here()
parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "CombinationExperiments")

figure_path <- file.path(parent_path, "Plots")

# Get a list of all subfolders in the current directory
subfolders <- list.dirs(path = new_path, full.names = TRUE, recursive = FALSE)

# Loop through each subfolder
for (folder in subfolders) {
  # Find the count and metadata files in the current subfolder
  area_file <- list.files(path = folder, pattern = "objectarea\\.csv$", full.names = TRUE)
  metadata_file <- list.files(path = folder, pattern = "metadata\\.csv$", full.names = TRUE)
  
  if (length(area_file) > 0 && length(metadata_file) > 0) {
    # Read the count data
    area_data <- read.csv(area_file)
    
    # Read the metadata
    metadata <- read.csv(metadata_file)
    
    # Add metadata columns to count data
    updated_area_data <- left_join(area_data, metadata, by = "XY")
    
    # Append to master dataframe
    combination_experiment_data <- rbind(combination_experiment_data, updated_area_data)
    all_metadata <- rbind(all_metadata, metadata)
    
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

# Define the custom labels for the x-axis
time_labels <- c("5" = "4 weeks", "6" = "5 weeks", "7" = "6 weeks")

#remove columns that don't contain relevant data
combination_experiment_data = combination_experiment_data %>% select(-c(Name,Channel,Z,Shape,XY))

# Replace Default blob or blank values in the Tags column with "Live Colony"
combination_experiment_data$Tags[combination_experiment_data$Tags == "Default blob"] <- "Live Colony"
combination_experiment_data$Tags[combination_experiment_data$Tags == ""] <- "Live Colony"
combination_experiment_data <- combination_experiment_data %>% mutate(Tags = replace_na(Tags, "Live Colony"))

##

## Task 1 - Make the data frame wider to get a count of the number of annotations per condition per replicate, add metadata

combination_experiment_data = combination_experiment_data %>% select(-c(All.Blob.metrics...Solidity,All.Blob.metrics...Rectangularity, All.Blob.metrics...Perimeter, All.Blob.metrics...Elongation,
                                                                        All.Blob.metrics...Eccentricity, All.Blob.metrics...Convexity, All.Blob.metrics...Circularity, All.Blob.metrics...Centroid...y,
                                                                        All.Blob.metrics...Centroid...x))

combination_experiment_data = combination_experiment_data %>% rename(Area = All.Blob.metrics...Area)

#filter out the colonies that are dead
dead_combination_experiment_data <- combination_experiment_data %>% filter(str_detect(Tags, "Dead Colony"))
live_combination_experiment_data <- combination_experiment_data %>% filter(str_detect(Tags, "Live Colony"))

# Group by XY and Time, then count unique IDs
areadata <- live_combination_experiment_data %>%
  group_by(Time, Replicate, Condition, Drug1Name, Drug1Dose, Drug2Name, Drug2Dose, Set, Plate) %>%
  summarise(AverageArea = mean(Area)) %>%
  ungroup()

normalized_areadata <- areadata %>%
  group_by(Set, Plate, Replicate, Condition) %>%
  mutate(NormalizedArea = AverageArea / AverageArea[Time == 5]) %>%
  ungroup()

# Plot combinations 

ggplot((normalized_areadata %>% filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% group_by(Time, Condition) %>%
          summarize(mean_area = mean(NormalizedArea, na.rm = TRUE), se = sd(NormalizedArea, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_area, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - se, ymax = mean_area + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.) \n+ 0.1 µM SB273005 (Integrin inhib.)"),
              y = mean_area*1.03), # Position above the line
            hjust = 1.1, vjust = 0.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.)"),
              y = mean_area*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
           aes(label = case_when(
             Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"),
           y = mean_area*1.03), # Position above the line
           hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
              y = mean_area*1.03), # Position above the line
            hjust = -0.15, vjust = -1.15, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Size Over Time by Condition",
       x = "Time",
       y = "Average Live Colony Size, Normalized to Size at 4 Weeks") +
  scale_x_continuous(breaks = 5:7, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = normalized_areadata %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")),
             aes(x = Time, y = NormalizedArea, color = Condition), 
             size = 3, alpha = 0.5) 

ggsave(file.path(figure_path,"deguelinsb273005_combination_normalizedarea.pdf"), width = 12.48, height = 8, units = "in")

ggplot((areadata %>% filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% group_by(Time, Condition) %>%
          summarize(mean_area = mean(AverageArea, na.rm = TRUE), se = sd(AverageArea, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_area, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - se, ymax = mean_area + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.) \n+ 0.1 µM SB273005 (Integrin inhib.)"),
              y = mean_area*1.03), # Position above the line
            hjust = 1.1, vjust = 0.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.)"),
              y = mean_area*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"),
              y = mean_area*1.03), # Position above the line
            hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
              y = mean_area*1.03), # Position above the line
            hjust = -0.15, vjust = -1.15, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Size Over Time by Condition",
       x = "Time",
       y = "Average Live Colony Size") +
  scale_x_continuous(breaks = 5:7, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = areadata %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")),
             aes(x = Time, y = AverageArea, color = Condition), 
             size = 3, alpha = 0.5) 

ggsave(file.path(figure_path,"deguelinsb273005_combination_area.pdf"), width = 12.48, height = 8, units = "in")

ggplot((areadata %>% filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% group_by(Time, Condition) %>%
          summarize(mean_area = mean(AverageArea, na.rm = TRUE), se = sd(AverageArea, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_area, color = Condition)) +
  geom_hline(data = . %>% filter(Time == 5),
             aes(yintercept = mean_area),
             color = "red", linetype = "dashed", linewidth = 1) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - se, ymax = mean_area + se), 
                width = 0.1, linewidth = 0.8) +
  facet_wrap(~ Condition, 
             labeller = labeller(Condition = c(
               "1 \xb5M Deguelin" = "1 µM Deguelin\n(PI3K/Akt inhib.)",
               "0.1 \xb5M SB273005" = "0.1 µM SB273005\n(Integrin inhib.)",
               "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005\n(PI3K/Akt + Integrin inhib.)",
               "Vemurafenib Only" = "Vemurafenib Only"
             ))) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Size Over Time by Condition",
       x = "Time",
       y = "Average Live Colony Size") +
  scale_x_continuous(breaks = 5:7, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = areadata %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")),
             aes(x = Time, y = AverageArea, color = Condition), 
             size = 3, alpha = 0.5) 

ggsave(file.path(figure_path,"deguelinsb273005_combination_facetedarea.pdf"), width = 16.2, height = 8, units = "in")

