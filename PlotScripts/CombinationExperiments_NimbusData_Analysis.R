library(tidyverse)
library(dplyr)
library(here)
library(readxl)
library(ggplot2)
library(ggpubr)

##

# This is a script to analyze colony count data obtained via NimbusImage. This script is designed to analyze all data taken for combination experiments

# Created by GTB on 4/4/2024

# Last edited by GTB on 4/17/2025

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
  count_file <- list.files(path = folder, pattern = "counts\\.csv$", full.names = TRUE)
  metadata_file <- list.files(path = folder, pattern = "metadata\\.csv$", full.names = TRUE)
  
  if (length(count_file) > 0 && length(metadata_file) > 0) {
    # Read the count data
    count_data <- read.csv(count_file)
    
    # Read the metadata
    metadata <- read.csv(metadata_file)
    
    # Add metadata columns to count data
    updated_count_data <- left_join(count_data, metadata, by = "XY")
    
    # Append to master dataframe
    combination_experiment_data <- rbind(combination_experiment_data, updated_count_data)
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
time_labels <- c("5" = "4 weeks", "6" = "5 weeks", "7" = "6 weeks", "8" = "7 weeks", "9" = "8 weeks")

#remove columns that don't contain relevant data
combination_experiment_data = combination_experiment_data %>% select(-c(Name,Channel,Z,Shape,XY))

# Replace Default blob or blank values in the Tags column with "Live Colony"
combination_experiment_data$Tags[combination_experiment_data$Tags == "Default blob"] <- "Live Colony"
combination_experiment_data$Tags[combination_experiment_data$Tags == ""] <- "Live Colony"
combination_experiment_data <- combination_experiment_data %>% mutate(Tags = replace_na(Tags, "Live Colony"))

##

## Task 1 - Make the data frame wider to get a count of the number of annotations per condition per replicate, add metadata

#filter out the colonies that are dead
dead_combination_experiment_data <- combination_experiment_data %>% filter(str_detect(Tags, "Dead Colony"))
live_combination_experiment_data <- combination_experiment_data %>% filter(str_detect(Tags, "Live Colony"))

# Group by XY and Time, then count unique IDs
countdata <- live_combination_experiment_data %>%
  group_by(Time, Replicate, Condition, Drug1Name, Drug1Dose, Drug2Name, Drug2Dose, Set, Plate) %>%
  summarise(Count = n_distinct(Id)) %>%
  ungroup()

#Add Counts of 0 for any xy/time combo that does not exist in the count data table
expanded_metadata <- all_metadata %>%
  crossing(Time = 5:9) %>%
  mutate(Count = 0)

complete_countdata <- full_join(
  countdata,
  expanded_metadata,
  by = c("Time", "Replicate", "Condition", "Drug1Name", "Drug1Dose", 
         "Drug2Name", "Drug2Dose", "Set", "Plate")
)

complete_countdata <- complete_countdata %>% filter(Time %in% 5:9)

complete_countdata <- complete_countdata %>% mutate(Count.x = ifelse(is.na(Count.x), 0, Count.x)) %>% mutate(Count = Count.x + Count.y) %>% select(-c(XY, Count.x, Count.y))

# Plot combinations 

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.) \n+ 0.1 µM SB273005 (Integrin inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = 1.1, vjust = 0.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
           aes(label = case_when(
             Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"),
           y = mean_count*1.03), # Position above the line
           hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -1.15, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
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
  geom_point(data = complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% 
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
             
             # t-tests
             combo_vs_deguelin_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
             )$p.value
             
             combo_vs_sb_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
             )$p.value
             
             paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                    "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
             
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
             
             # t-tests
             combo_vs_deguelin_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
             )$p.value
             
             combo_vs_sb_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
             )$p.value
             
             paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                    "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
             
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
             
             # t-tests
             combo_vs_deguelin_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
             )$p.value
             
             combo_vs_sb_p <- t.test(
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin"],
               time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
             )$p.value
             
             paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                    "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
             
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"deguelinsb273005_combination.pdf"), width = 12.48, height = 8, units = "in")

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%
          filter(Condition == "Vemurafenib Only")), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.) \n+ 0.1 µM SB273005 (Integrin inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = 1.1, vjust = 0.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -1.15, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_manual(values = c(
    "Vemurafenib Only" = "#C77CFF",         
    "0.1 \xb5M SB273005" = "#F8766D",       
    "1 \xb5M Deguelin" = "#00BFC4",         
    "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "#7CAE00"  
  ),
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none")

ggsave(file.path(figure_path,"deguelinsb273005_vemonly_combination.pdf"), width = 12.48, height = 8, units = "in")

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4) | (Set == 3 & Plate == 1)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%
          filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin"))), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.) \n+ 0.1 µM SB273005 (Integrin inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = 1.1, vjust = 0.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin (PI3K/Akt inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -1.15, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_manual(values = c(
    "Vemurafenib Only" = "#C77CFF",         
    "0.1 \xb5M SB273005" = "#F8766D",       
    "1 \xb5M Deguelin" = "#00BFC4",         
    "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "#7CAE00"  
  ),
  labels = c(
    "1 \xb5M Deguelin" = "1 µM Deguelin",
    "0.1 \xb5M SB273005" = "0.1 µM SB273005",
    "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
    "Vemurafenib Only" = "Vemurafenib Only"
  )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none")

ggsave(file.path(figure_path,"deguelinsb273005_vemandsingledrugs_combination.pdf"), width = 12.48, height = 8, units = "in")

ggplot((complete_countdata %>% filter((Set == 7 & Plate == 1) | (Set == 7 & Plate == 2)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005")) %>% filter(Time %in% 5:7) %>%
          mutate(Condition = factor(Condition, levels = c("0.1 \xb5M SB273005", 
                                                          "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005", 
                                                          "3.16 \xb5M IWP-O1",
                                                          "Vemurafenib Only"))) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005 \n(Integrin inhib.)"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005" ~ "3.16 µM IWP-O1 (Porcn inhib.) + \n0.1 µM SB273005 (Integrin inhib.)"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0, vjust = -0.4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1" ~ "3.16 µM IWP-O1 (Porcn inhib.)",
              ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
            y = mean_count*1.03), # Position above the line
            hjust = -0.15, vjust = -1.8, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "3.16 \xb5M IWP-O1" = "3.16 µM IWP-O1",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005" = "3.16 µM IWP-O1 + 0.1 µM SB273005",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() + 
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 7 & Plate == 1) | (Set == 7 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005")) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 7 & Plate == 1) | (Set == 7 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               combo_vs_sb_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
               )$p.value
               
               paste0("Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)), "\n",
                      "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 7 & Plate == 1) | (Set == 7 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               combo_vs_sb_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
               )$p.value
               
               paste0("Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)), "\n",
                      "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 7 & Plate == 1) | (Set == 7 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "0.1 \xb5M SB273005", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               combo_vs_sb_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 0.1 \xb5M SB273005"],
                 time_data$Count[time_data$Condition == "0.1 \xb5M SB273005"]
               )$p.value
               
               paste0("Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)), "\n",
                      "Combo vs SB273005: p=", ifelse(combo_vs_sb_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sb_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"iwpo1sb273005_combination.pdf"), width = 12.48, height = 8, units = "in")

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "5 nM YM155", "5 nM GMX1778 + 5 nM YM155")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only",
              Condition == "5 nM YM155" ~ "5 nM YM155"),
            y = mean_count*1.03), # Position above the line
            hjust = 0.15, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "5 nM GMX1778 + 5 nM YM155" ~ "5 nM GMX1778 + 5 nM YM155"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 1, vjust = 2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "5 nM GMX1778" ~ "5 nM GMX1778"
              ),
            y = mean_count*1.03), # Position above the line
            hjust = -.5, vjust = 2, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "5 nM YM155" = "5 nM YM155",
      "5 nM GMX1778" = "5 nM GMX1778",
      "5 nM YM155 + 5 nM GMX1778" = "5 nM YM155 + 5 nM GMX1778",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "5 nM YM155", "5 nM GMX1778 + 5 nM YM155")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "5 nM YM155", "5 nM GMX1778 + 5 nM YM155")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "5 nM YM155", "5 nM GMX1778 + 5 nM YM155")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "5 nM YM155", "5 nM GMX1778 + 5 nM YM155")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "5 nM GMX1778 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"ym155gmx1778_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM GMX1778")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only",
              Condition == "5 nM GMX1778" ~ "5 nM GMX1778",
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.3, vjust = 1.25, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin + 5 nM GMX1778" ~ "5 nM GMX1778 + 1 µM Deguelin"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 6, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "5 nM GMX1778" = "5 nM GMX1778",
      "1 \xb5M Deguelin + 5 nM GMX1778" = "1 µM Deguelin + 5 nM GMX1778",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM GMX1778")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM GMX1778")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM GMX1778")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM GMX1778")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"deguelingmx1778_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM YM155")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin + 5 nM YM155" ~ "1 µM Deguelin + 5 nM YM155"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.05, vjust = -0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin"),
            y = mean_count*1.03), # Position above the line
            hjust = -0.25, vjust = -0.75, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "5 nM YM155" ~ "5 nM YM155"),
            y = mean_count*1.03), # Position above the line
            hjust = 0.1, vjust = 3.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = -2, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "5 nM YM155" = "5 nM YM155",
      "1 \xb5M Deguelin + 5 nM YM155" = "1 µM Deguelin + 5 nM YM155",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "1 \xb5M Deguelin", "1 \xb5M Deguelin + 5 nM YM155")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "5 nM YM155", "1 \xb5M Deguelin + 5 nM YM155")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "5 nM YM155", "1 \xb5M Deguelin + 5 nM YM155")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "5 nM YM155", "1 \xb5M Deguelin + 5 nM YM155")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deguelin_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deguelin_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deguelin_p)), "\n",
                      "Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"deguelinym155_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM GMX1778")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only",
              Condition == "5 nM GMX1778" ~ "5 nM GMX1778",
              Condition == "3.16 \xb5M IWP-O1" ~ "3.16 µM IWP-O1"),
            y = mean_count*1.03), # Position above the line
            hjust = -0.3, vjust = 1.25, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778" ~ "5 nM GMX1778 + 3.16 µM IWP-O1"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 5, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "3.16 \xb5M IWP-O1" = "3.16 µM IWP-O1",
      "5 nM GMX1778" = "5 nM GMX1778",
      "3.16 \xb5M IWP-O1 + 5 nM GMX1778" = "3.16 µM IWP-O1 + 5 nM GMX1778",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM GMX1778")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM GMX1778")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM GMX1778")) %>%
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 1) | (Set == 2 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM GMX1778")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)
  

ggsave(file.path(figure_path,"iwpo1gmx1778_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M SR18662", "1 \xb5M SR18662 + 5 nM GMX1778")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only",
              Condition == "5 nM GMX1778" ~ "5 nM GMX1778",
              Condition == "1 \xb5M SR18662" ~ "1 µM SR18662"),
            y = mean_count*1.03), # Position above the line
            hjust = -0.3, vjust = 1.25, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M SR18662 + 5 nM GMX1778" ~ "5 nM GMX1778 + 1 µM SR18662"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 2, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M SR18662" = "1 µM SR18662",
      "5 nM GMX1778" = "5 nM GMX1778",
      "1 \xb5M SR18662 + 5 nM GMX1778" = "1 µM SR18662 + 5 nM GMX1778",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M SR18662", "1 \xb5M SR18662 + 5 nM GMX1778")) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M SR18662", "1 \xb5M SR18662 + 5 nM GMX1778")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M SR18662", "1 \xb5M SR18662 + 5 nM GMX1778")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 3) | (Set == 1 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM GMX1778", "1 \xb5M SR18662", "1 \xb5M SR18662 + 5 nM GMX1778")) %>%  
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_gmx_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "5 nM GMX1778"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 5 nM GMX1778"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs GMX1778: p=", ifelse(combo_vs_gmx_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_gmx_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"sr18662gmx1778_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM YM155")) %>% filter(Time %in% 5:7) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = -3, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "5 nM YM155" ~ "5 nM YM155"
             ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1" ~ "3.16 µM IWP-O1"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.5, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155" ~ "3.16 µM IWP-O1 + \n        5 nM YM155"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.25, vjust = 1.5, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "3.16 \xb5M IWP-O1" = "3.16 µM IWP-O1",
      "5 nM YM155" = "5 nM YM155",
      "3.16 \xb5M IWP-O1 + 5 nM YM155" = "3.16 µM IWP-O1 + 5 nM YM155",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM YM155")) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM YM155")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM YM155")) %>%
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 1 & Plate == 1) | (Set == 1 & Plate == 2)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "5 nM YM155", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 5 nM YM155")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_ym155_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "5 nM YM155"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 5 nM YM155"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs YM155: p=", ifelse(combo_vs_ym155_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_ym155_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"iwpo1ym155_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "1 \xb5M SR18662", "1 \xb5M SR18662 + 1 \xb5M Deguelin")) %>% filter(Time %in% 5:7) %>%
          mutate(Condition = factor(Condition, levels = c("1 \xb5M Deguelin", 
                                                          "1 \xb5M SR18662 + 1 \xb5M Deguelin", 
                                                          "1 \xb5M SR18662",
                                                          "Vemurafenib Only"))) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = -1.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin"),
            y = mean_count*1.03), # Position above the line
            hjust = 1.1, vjust = 3, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M SR18662" ~ "1 µM SR18662"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin" ~ "1 µM Deguelin + 1 µM SR18662"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 0.1, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M SR18662" = "1 µM SR18662",
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "1 \xb5M SR18662 + 1 \xb5M Deguelin" = "1 µM SR18662 + 1 µM Deguelin",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "1 \xb5M SR18662", "1 \xb5M SR18662 + 1 \xb5M Deguelin")) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "1 \xb5M SR18662", "1 \xb5M SR18662 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 5)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "1 \xb5M SR18662", "1 \xb5M SR18662 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 6)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 2 & Plate == 3) | (Set == 2 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M Deguelin", "1 \xb5M SR18662", "1 \xb5M SR18662 + 1 \xb5M Deguelin")) %>%   
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_sr_p <- t.test(
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M SR18662"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs SR18662: p=", ifelse(combo_vs_sr_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_sr_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"sr18662deguelin_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 9 & Plate == 3) | (Set == 9 & Plate == 4)) %>% 
          filter(Condition %in% c("Dabrafenib/Trametinib Only", "1 \xb5M Deguelin", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin")) %>% filter(Time %in% 7:9) %>%
          mutate(Condition = factor(Condition, levels = c("1 \xb5M Deguelin", 
                                                          "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin", 
                                                          "3.16 \xb5M IWP-O1",
                                                          "Dabrafenib/Trametinib Only"))) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "Dabrafenib/Trametinib Only" ~ "Dabrafenib/Trametinib Only"
              ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.1, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin"),
            y = mean_count*1.03), # Position above the line
            hjust = 0.1, vjust = 4.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1" ~ "3.16 µM IWP-O1"),
            y = mean_count*1.03), # Position above the line
            hjust = -0.4, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin" ~ "1 µM Deguelin + \n3.16 µM IWP-O1"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 1.25, vjust = 1.1, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "3.16 \xb5M IWP-O1" = "3.16 µM IWP-O1",
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin" = "3.16 µM IWP-O1 + 1 µM Deguelin",
      "Dabrafenib/Trametinib Only" = "Dabrafenib/Trametinib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 9 & Plate == 3) | (Set == 9 & Plate == 4)) %>% 
               filter(Condition %in% c("Dabrafenib/Trametinib Only", "1 \xb5M Deguelin", "3.16 \xb5M IWP-O1", "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin")) %>%
               filter(Time %in% c(7:9)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5) +
  # Time 5 ANOVA + t-tests
  annotate("text", x = 4.5, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 9 & Plate == 3) | (Set == 9 & Plate == 4)) %>% 
               filter(Condition %in% c("Dabrafenib/Trametinib Only", "1 \xb5M Deguelin", "3.16 \xb5M IWP-O", "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 7)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 6 ANOVA + t-tests
  annotate("text", x = 5.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 9 & Plate == 3) | (Set == 9 & Plate == 4)) %>% 
               filter(Condition %in% c("Dabrafenib/Trametinib Only", "1 \xb5M Deguelin", "3.16 \xb5M IWP-O", "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin")) %>% 
               filter(Time == 8)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6) +
  
  # Time 7 ANOVA + t-tests
  annotate("text", x = 6.25, y = Inf, 
           label = {
             time_data <- complete_countdata %>% 
               filter((Set == 9 & Plate == 3) | (Set == 9 & Plate == 4)) %>% 
               filter(Condition %in% c("Dabrafenib/Trametinib Only", "1 \xb5M Deguelin", "3.16 \xb5M IWP-O", "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin")) %>%   
               filter(Time == 9)
             
             # ANOVA
             anova_result <- aov(Count ~ Condition, data = time_data)
             anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
             
             if(anova_p < 0.001) {
               
               # t-tests
               combo_vs_deg_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "1 \xb5M Deguelin"]
               )$p.value
               
               combo_vs_iwpo1_p <- t.test(
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1 + 1 \xb5M Deguelin"],
                 time_data$Count[time_data$Condition == "3.16 \xb5M IWP-O1"]
               )$p.value
               
               paste0("Combo vs Deguelin: p=", ifelse(combo_vs_deg_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_deg_p)), "\n",
                      "Combo vs IWP-O1: p=", ifelse(combo_vs_iwpo1_p < 0.001, "<0.001", sprintf("%.3f", combo_vs_iwpo1_p)))
               
             } else {}
           },
           vjust = 1, hjust = 0, size = 6)

ggsave(file.path(figure_path,"dt_iwpo1deguelin_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 9 & Plate == 1) | (Set == 9 & Plate == 2)) %>% 
          filter(Condition %in% c("Dabrafenib/Trametinib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>% filter(Time %in% 7:9) %>%
          mutate(Condition = factor(Condition, levels = c("1 \xb5M Deguelin", 
                                                          "0.1 \xb5M SB273005 + 1 \xb5M Deguelin", 
                                                          "0.1 \xb5M SB273005",
                                                          "Dabrafenib/Trametinib Only"))) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "Dabrafenib/Trametinib Only" ~ "Dabrafenib/Trametinib Only"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.1, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "1 \xb5M Deguelin" ~ "1 µM Deguelin"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.1, vjust = 4.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005" ~ "0.1 µM SB273005"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.4, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 8), 
            aes(label = case_when(
              Condition == "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" ~ "1 µM Deguelin + \n0.1 µM SB273005"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 1.25, vjust = 1.1, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M Deguelin" = "1 µM Deguelin",
      "0.1 \xb5M SB273005" = "0.1 µM SB273005",
      "0.1 \xb5M SB273005 + 1 \xb5M Deguelin" = "1 µM Deguelin + 0.1 µM SB273005",
      "Dabrafenib/Trametinib Only" = "Dabrafenib/Trametinib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 9 & Plate == 1) | (Set == 9 & Plate == 2)) %>% 
               filter(Condition %in% c("Dabrafenib/Trametinib Only", "0.1 \xb5M SB273005", "1 \xb5M Deguelin", "0.1 \xb5M SB273005 + 1 \xb5M Deguelin")) %>%
               filter(Time %in% c(7:9)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"dt_sb273005deguelin_combination.pdf"), width = 10, height = 9.6, units = "in")

ggplot((complete_countdata %>% filter((Set == 8 & Plate == 3) | (Set == 8 & Plate == 4)) %>% 
          filter(Condition %in% c("Vemurafenib Only", "1 \xb5M SR18662", "1 \xb5M S3I-201", "1 \xb5M SR18662 + 1 \xb5M S3I-201")) %>% filter(Time %in% 5:7) %>%
          mutate(Condition = factor(Condition, levels = c("1 \xb5M SR18662", 
                                                          "1 \xb5M SR18662 + 1 \xb5M S3I-201", 
                                                          "1 \xb5M S3I-201",
                                                          "Vemurafenib Only"))) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib Only"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = -0.1, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M SR18662" ~ "1 µM SR18662"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.1, vjust = 4.5, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M S3I-201" ~ "1 µM S3I-201"),
              y = mean_count*1.03), # Position above the line
            hjust = -0.4, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 6), 
            aes(label = case_when(
              Condition == "1 \xb5M SR18662 + 1 \xb5M S3I-201" ~ "1 µM SR18662 + \n1 µM S3I-201"
            ),
            y = mean_count*1.03), # Position above the line
            hjust = 1.25, vjust = 1.1, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Colony Count Over Time by Condition",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition") +
  scale_x_continuous(breaks = 5:9, labels = time_labels) +
  scale_color_discrete(
    labels = c(
      "1 \xb5M SR18662" = "1 µM SR18662",
      "1 \xb5M S3I-201" = "1 µM S3I-201",
      "1 \xb5M SR18662 + 1 \xb5M S3I-201" = "1 µM SR18662 + 1 µM S3I-201",
      "Vemurafenib Only" = "Vemurafenib Only"
    )
  ) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = complete_countdata %>% 
               filter((Set == 8 & Plate == 3) | (Set == 8 & Plate == 4)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "1 \xb5M SR18662", "1 \xb5M S3I-201", "1 \xb5M SR18662 + 1 \xb5M S3I-201")) %>%
               filter(Time %in% c(5:7)),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"sr18662s3i201_combination.pdf"), width = 10, height = 9.6, units = "in")

##



