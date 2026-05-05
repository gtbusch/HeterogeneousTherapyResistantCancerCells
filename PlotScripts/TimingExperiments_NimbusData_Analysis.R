library(tidyverse)
library(dplyr)
library(here)
library(readxl)
library(ggplot2)
library(ggpubr)

##

# This is a script to analyze colony count data obtained via NimbusImage. This script is designed to analyze all data taken for timing experiments

# Created by GTB on 4/4/2024

# Last edited by GTB on 4/8/2025

##

## Task 0 - Load the Data

# Initialize an empty master dataframe to store all combined data
timing_experiment_data <- data.frame()

current_path <- here()
parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "TimingExperiments")

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
    timing_experiment_data <- rbind(timing_experiment_data, updated_count_data)
    
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
time_labels <- c("3" = "2 weeks", "4" = "3 weeks", "5" = "4 weeks", "6" = "5 weeks")

timing_colors <- c(
    "Vemurafenib Only" = "#F8766D",       # Default ggplot red
    "3 Weeks of Vemurafenib" = "#7CAE00", # Default ggplot green
    "Combination Treatment" = "#00BFC4",  # Default ggplot cyan
    "1 Week of Vemurafenib" = "#C77CFF",       # Default ggplot purple
    "2 Weeks of Vemurafenib" = "#FF8C00"         # Default ggplot orange
  )

condition_labels <- c(
  "Vemurafenib Only" = "Vemurafenib for 5 Weeks",       
  "3 Weeks of Vemurafenib" = "Vem. to Combination at Week 3", 
  "Combination Treatment" = "Combination for 5 Weeks",  
  "1 Week of Vemurafenib" = "Vem. to Combination at Week 1",       
  "2 Weeks of Vemurafenib" = "Vem. to Combination at Week 2"         
)

#remove columns that don't contain relevant data
timing_experiment_data = timing_experiment_data %>% select(-c(Name,Channel,Z,Shape,XY))

#remove conditions that won't be included in the analysis
timing_experiment_data = timing_experiment_data %>% filter(Condition %in% c("Vemurafenib Only", "1 Week of Vemurafenib", "2 Weeks of Vemurafenib", "3 Weeks of Vemurafenib", "Combination Treatment"))

# Replace Default blob or blank values in the Tags column with "Live Colony"
timing_experiment_data$Tags[timing_experiment_data$Tags == "Default blob"] <- "Live Colony"
timing_experiment_data$Tags[timing_experiment_data$Tags == ""] <- "Live Colony"
timing_experiment_data <- timing_experiment_data %>% mutate(Tags = replace_na(Tags, "Live Colony"))

##

## Task 1 - Make the data frame wider to get a count of the number of annotations per condition per replicate, add metadata

#filter out the colonies that are dead
dead_timing_experiment_data <- timing_experiment_data %>% filter(Tags == "Dead Colony")
live_timing_experiment_data <- timing_experiment_data %>% filter(Tags != "Dead Colony")

# Group by XY and Time, then count unique IDs
countdata <- live_timing_experiment_data %>%
  group_by(Time, Replicate, Condition, DrugName, Dose, Set) %>%
  summarise(Count = n_distinct(Id)) %>%
  ungroup()

#adjust the time indeces for Saracatinib to match the rest of the datasets
countdata <- countdata %>% mutate(Time = ifelse((DrugName == "Saracatinib" & Set == 1), 
                       Time / 2 + 0.5, Time))

#Add Counts of 0 for any xy/time combo that does not exist in the count data table
countdata <- countdata %>%
  complete(Replicate = 1:2, Condition = unique(countdata$Condition), 
           Time = 3:6, fill = list(Count = 0), DrugName = c("Dasatinib", "(+)-JQ1", "S3I-201", "SR18662"))

countdata <- countdata %>%
  complete(Replicate = 1:2, Condition = c("1 Week of Vemurafenib","2 Weeks of Vemurafenib","3 Weeks of Vemurafenib", "Vemurafenib Only"), 
           Time = 3:6, fill = list(Count = 0), DrugName = c("Saracatinib"), Set = 1)

countdata <- countdata %>%
  complete(Replicate = 1:6, Condition = c("Combination Treatment"), 
           Time = 3:6, fill = list(Count = 0), DrugName = c("Saracatinib"), Set = 1)

countdata <- countdata %>%
  complete(Replicate = 1:2, Condition = c("Combination Treatment", "1 Week of Vemurafenib","2 Weeks of Vemurafenib","3 Weeks of Vemurafenib", "Vemurafenib Only"), 
           Time = 3:6, fill = list(Count = 0), DrugName = c("Saracatinib"), Set = 2)

countdata <- countdata %>%
  complete(Replicate = 1:2, Condition = c("Combination Treatment", "1 Week of Vemurafenib","2 Weeks of Vemurafenib","3 Weeks of Vemurafenib", "Vemurafenib Only"), 
           Time = 3:6, fill = list(Count = 0), DrugName = c("Saracatinib"), Set = 3)

#Make plots showing the following conditions: Vemurafenib Only and 3 Weeks of Vemurafenib ----

ggplot((countdata %>% filter(DrugName == "Saracatinib") %>% filter(Time %in% c(4:6)) %>% filter (Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "3 Weeks of Vemurafenib" ~ "Vem. to Combination at Week 3"),
            y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 3, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib for 5 Weeks"),
              y = mean_count*1.03), # Position above the line
            hjust = .9, vjust = -1, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Timing of 1 µM Saracatinib \nand Vemurafenib",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_continuous(breaks = 3:6, labels = time_labels) +
  scale_color_manual(values = timing_colors, labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = countdata %>% 
               filter(DrugName == "Saracatinib") %>% 
               filter(Time %in% c(4:6)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"saracatinib_3weekonly_timingexperiment.pdf"), width = 7.36, height = 6.4, units = "in")

ggplot((countdata %>% filter(DrugName == "S3I-201") %>% filter(Time %in% c(4:6)) %>% filter (Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "3 Weeks of Vemurafenib" ~ "Vem. to Combination at Week 3"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 0.1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib for 5 Weeks"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 5, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Timing of 1 µM S3I-201 \nand Vemurafenib",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_continuous(breaks = 3:6, labels = time_labels) +
  scale_color_manual(values = timing_colors, labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = countdata %>% 
               filter(DrugName == "S3I-201") %>% 
               filter(Time %in% c(4:6)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"s3i201_3weekonly_timingexperiment.pdf"), width = 7.36, height = 6.44, units = "in")

ggplot((countdata %>% filter(DrugName == "SR18662") %>% filter(Time %in% c(4:6)) %>% filter (Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "3 Weeks of Vemurafenib" ~ "Vem. to Combination at Week 3"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 4, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib for 5 Weeks"),
              y = mean_count*1.03), # Position above the line
            hjust = .75, vjust = -2, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Timing of 1 µM SR18662 \nand Vemurafenib",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_continuous(breaks = 3:6, labels = time_labels) +
  scale_color_manual(values = timing_colors, labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = countdata %>% 
               filter(DrugName == "SR18662") %>% 
               filter(Time %in% c(4:6)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"sr18662_3weekonly_timingexperiment.pdf"), width = 7.36, height = 6.4, units = "in")

ggplot((countdata %>% filter(DrugName == "Dasatinib") %>% filter(Time %in% c(4:6)) %>% filter (Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "3 Weeks of Vemurafenib" ~ "Vem. to Combination at Week 3"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = -2, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib for 5 Weeks"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = -2.75, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Timing of 0.1 µM Dasatinib \nand Vemurafenib",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_continuous(breaks = 3:6, labels = time_labels) +
  scale_color_manual(values = timing_colors, labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = countdata %>% 
               filter(DrugName == "Dasatinib") %>% 
               filter(Time %in% c(4:6)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"dasatinib_3weekonly_timingexperiment.pdf"), width = 7.36, height = 6.4, units = "in")

ggplot((countdata %>% filter(DrugName == "(+)-JQ1") %>% filter(Time %in% c(4:6)) %>% filter (Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")) %>% group_by(Time, Condition) %>%
          summarize(mean_count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Time, y = mean_count, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se), 
                width = 0.1, linewidth = 0.8) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "Vemurafenib Only" ~ "Vemurafenib for 5 Weeks"),
              y = mean_count*1.03), # Position above the line
            hjust = .75, vjust = -1, size = 8, show.legend = FALSE) +
  geom_text(data = . %>% group_by(Condition) %>% filter(Time == 5), 
            aes(label = case_when(
              Condition == "3 Weeks of Vemurafenib" ~ "Vem. to Combination at Week 3"),
              y = mean_count*1.03), # Position above the line
            hjust = 0.5, vjust = 2.5, size = 8, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Timing of 0.1 µM (+)-JQ1 \nand Vemurafenib",
       x = "Time",
       y = "Number of Live Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_continuous(breaks = 3:6, labels = time_labels) +
  scale_color_manual(values = timing_colors, labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none") +
  geom_point(data = countdata %>% 
               filter(DrugName == "(+)-JQ1") %>% 
               filter(Time %in% c(4:6)) %>% 
               filter(Condition %in% c("Vemurafenib Only", "3 Weeks of Vemurafenib")),
             aes(x = Time, y = Count, color = Condition), 
             size = 3, alpha = 0.5)

ggsave(file.path(figure_path,"jq1_3weekonly_timingexperiment.pdf"), width = 7.36, height = 6.4, units = "in")



#Make plots showing all conditions at the final timepoint ----

timing_comparisons = list(c("Combination Treatment", "1 Week of Vemurafenib"), c("Combination Treatment", "Vemurafenib Only"), c("Combination Treatment", "2 Weeks of Vemurafenib"), c("Combination Treatment", "3 Weeks of Vemurafenib"),
                          c("3 Weeks of Vemurafenib", "1 Week of Vemurafenib"), c("Vemurafenib Only", "1 Week of Vemurafenib"),
                          c("2 Weeks of Vemurafenib", "3 Weeks of Vemurafenib"),c("2 Weeks of Vemurafenib", "Vemurafenib Only"), c("3 Weeks of Vemurafenib", "Vemurafenib Only"))

ggplot((countdata %>% filter(DrugName == "Saracatinib") %>% filter(Time %in% c(6)) %>% group_by(Time, Condition) %>% 
          mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                          "1 Week of Vemurafenib", 
                                                          "2 Weeks of Vemurafenib", 
                                                          "3 Weeks of Vemurafenib", 
                                                          "Vemurafenib Only"))) %>%
          summarize(Count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Condition, y = Count)) +
  geom_bar(stat = "identity", width = 0.25, fill = "darkgreen") +
  geom_hline(yintercept = 26.166667, color = "red", linetype = "dashed") +   
  geom_errorbar(aes(ymin = Count - se, ymax = Count + se), 
                width = 0.1, linewidth = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "1 µM Saracatinib and Vemurafenib Timing", 
       x = "Time",
       y = "Final Number of Live \nResistant Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  geom_point(data = countdata %>% 
               filter(DrugName == "Saracatinib") %>% 
               filter(Time %in% c(6)),
             aes(x = Condition, y = Count), 
             size = 3, color = "black", alpha = 0.5) +
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "Saracatinib") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "anova", 
                     size = 6,
                     label = "p.format",
                     label.y = 55,
                     inherit.aes = FALSE) +
  # Pairwise t-tests with brackets
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "Saracatinib") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "t.test", 
                     comparisons = timing_comparisons,
                     size = 6,
                     label = "p.format",
                     step.increase = 0.125,
                     tip.length = 0.01)


ggsave(file.path(figure_path,"saracatinib_timingexperiment.pdf"), width = 24, height = 6.4, units = "in")

timing_comparisons = list(c("Combination Treatment", "1 Week of Vemurafenib"), c("Combination Treatment", "Vemurafenib Only"), c("Combination Treatment", "2 Weeks of Vemurafenib"), c("Combination Treatment", "3 Weeks of Vemurafenib"),
                          c("2 Weeks of Vemurafenib", "1 Week of Vemurafenib"), c("3 Weeks of Vemurafenib", "1 Week of Vemurafenib"), c("Vemurafenib Only", "1 Week of Vemurafenib"),
                          c("2 Weeks of Vemurafenib", "3 Weeks of Vemurafenib"),c("2 Weeks of Vemurafenib", "Vemurafenib Only"), c("3 Weeks of Vemurafenib", "Vemurafenib Only"))

ggplot((countdata %>% filter(DrugName == "S3I-201") %>%filter(Time %in% c(6)) %>% group_by(Time, Condition) %>% 
          mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                          "1 Week of Vemurafenib", 
                                                          "2 Weeks of Vemurafenib", 
                                                          "3 Weeks of Vemurafenib", 
                                                          "Vemurafenib Only"))) %>%
          summarize(Count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Condition, y = Count)) +
  geom_bar(stat = "identity", width = 0.25, fill = "darkgreen") +
  geom_hline(yintercept = 28.16667, color = "red", linetype = "dashed") +
  geom_errorbar(aes(ymin = Count - se, ymax = Count + se), 
                width = 0.1, linewidth = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
 labs(title = "1 µM S3I-201 and Vemurafenib Timing", 
       x = "Time",
       y = "Final Number of Live \nResistant Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none", panel.spacing = unit(.45, "in"), plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  geom_point(data = countdata %>% 
               filter(DrugName == "S3I-201") %>% 
               filter(Time %in% c(6)),
             aes(x = Condition, y = Count), 
             size = 3, color = "black", alpha = 0.5) +
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "S3I-201") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "anova", 
                     size = 6,
                     label = "p.format",
                     label.y = 45,
                     inherit.aes = FALSE)

ggsave(file.path(figure_path,"s3i201_timingexperiment.pdf"), width = 24, height = 5, units = "in")

timing_comparisons = list(c("Combination Treatment", "Vemurafenib Only"), c("Combination Treatment", "3 Weeks of Vemurafenib"),
                          c("3 Weeks of Vemurafenib", "1 Week of Vemurafenib"), c("Vemurafenib Only", "1 Week of Vemurafenib"),
                          c("2 Weeks of Vemurafenib", "3 Weeks of Vemurafenib"),c("2 Weeks of Vemurafenib", "Vemurafenib Only"), c("3 Weeks of Vemurafenib", "Vemurafenib Only"))

ggplot((countdata %>% filter(DrugName == "(+)-JQ1") %>%filter(Time %in% c(6)) %>% group_by(Time, Condition) %>% 
          mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                          "1 Week of Vemurafenib", 
                                                          "2 Weeks of Vemurafenib", 
                                                          "3 Weeks of Vemurafenib", 
                                                          "Vemurafenib Only"))) %>%
          summarize(Count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Condition, y = Count)) +
  geom_bar(stat = "identity", width = 0.25, fill = "darkgreen") +
  geom_hline(yintercept = 25.8333333, color = "red", linetype = "dashed") +   
  geom_errorbar(aes(ymin = Count - se, ymax = Count + se), 
                width = 0.1, linewidth = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "0.1 µM (+)-JQ1 and Vemurafenib Timing", 
       x = "Time",
       y = "Final Number of Live \nResistant Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  geom_point(data = countdata %>% 
               filter(DrugName == "(+)-JQ1") %>% 
               filter(Time %in% c(6)),
             aes(x = Condition, y = Count), 
             size = 3, color = "black", alpha = 0.5) +
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "(+)-JQ1") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "anova", 
                     size = 6,
                     label = "p.format",
                     label.y = 55,
                     inherit.aes = FALSE) +
  # Pairwise t-tests with brackets
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "(+)-JQ1") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "t.test", 
                     comparisons = timing_comparisons,
                     size = 6,
                     label = "p.format",
                     step.increase = 0.125,
                     tip.length = 0.01)

ggsave(file.path(figure_path,"jq1_timingexperiment.pdf"), width = 24, height = 6.4, units = "in")

timing_comparisons = list(c("Combination Treatment", "1 Week of Vemurafenib"), c("Combination Treatment", "Vemurafenib Only"), c("Combination Treatment", "2 Weeks of Vemurafenib"), c("Combination Treatment", "3 Weeks of Vemurafenib"),
                          c("Vemurafenib Only", "1 Week of Vemurafenib"), c("2 Weeks of Vemurafenib", "Vemurafenib Only"), c("3 Weeks of Vemurafenib", "Vemurafenib Only"))

ggplot((countdata %>% filter(DrugName == "SR18662") %>%filter(Time %in% c(6)) %>% group_by(Time, Condition) %>% 
          mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                          "1 Week of Vemurafenib", 
                                                          "2 Weeks of Vemurafenib", 
                                                          "3 Weeks of Vemurafenib", 
                                                          "Vemurafenib Only"))) %>%
          summarize(Count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Condition, y = Count)) +
  geom_bar(stat = "identity", width = 0.25, fill = "darkgreen") +
  geom_hline(yintercept = 22, color = "red", linetype = "dashed") + 
  geom_errorbar(aes(ymin = Count - se, ymax = Count + se), 
                width = 0.1, linewidth = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "1 µM SR18662 and Vemurafenib Timing", 
       x = "Time",
       y = "Final Number of Live \nResistant Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  geom_point(data = countdata %>% 
               filter(DrugName == "SR18662") %>% 
               filter(Time %in% c(6)),
             aes(x = Condition, y = Count), 
             size = 3, color = "black", alpha = 0.5) +
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "SR18662") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "anova", 
                     size = 6,
                     label = "p.format",
                     label.y = 40,
                     inherit.aes = FALSE) +
  # Pairwise t-tests with brackets
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "SR18662") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "t.test", 
                     comparisons = timing_comparisons,
                     size = 6,
                     label = "p.format",
                     step.increase = 0.125,
                     tip.length = 0.01)

ggsave(file.path(figure_path,"sr18662_timingexperiment.pdf"), width = 24, height = 6.4, units = "in")

timing_comparisons = list(c("Combination Treatment", "Vemurafenib Only"), c("Combination Treatment", "3 Weeks of Vemurafenib"),
                          c("3 Weeks of Vemurafenib", "1 Week of Vemurafenib"), c("Vemurafenib Only", "1 Week of Vemurafenib"),
                          c("2 Weeks of Vemurafenib", "Vemurafenib Only"), c("3 Weeks of Vemurafenib", "Vemurafenib Only"))


ggplot((countdata %>% filter(DrugName == "Dasatinib") %>%filter(Time %in% c(6)) %>% group_by(Time, Condition) %>% 
          mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                          "1 Week of Vemurafenib", 
                                                          "2 Weeks of Vemurafenib", 
                                                          "3 Weeks of Vemurafenib", 
                                                          "Vemurafenib Only"))) %>%
          summarize(Count = mean(Count, na.rm = TRUE), se = sd(Count, na.rm = TRUE) / sqrt(n()), .groups = 'drop')), 
       aes(x = Condition, y = Count)) +
  geom_bar(stat = "identity", width = 0.25, fill = "darkgreen") +
  geom_hline(yintercept = 32.75, color = "red", linetype = "dashed") +   
  geom_errorbar(aes(ymin = Count - se, ymax = Count + se), 
                width = 0.1, linewidth = 0.8) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "0.1 µM Dasatinib and Vemurafenib Timing", 
       x = "Time",
       y = "Final Number of Live \nResistant Colonies",
       color = "Condition",
       linetype = "Replicate") +
  scale_x_discrete(labels = condition_labels) +
  custom_figure_theme() +
  theme(legend.position = "none", plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  geom_point(data = countdata %>% 
               filter(DrugName == "Dasatinib") %>% 
               filter(Time %in% c(6)),
             aes(x = Condition, y = Count), 
             size = 3, color = "black", alpha = 0.5) +
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "Dasatinib") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "anova", 
                     size = 6,
                     label = "p.format",
                     label.y = 55,
                     inherit.aes = FALSE) +
  # Pairwise t-tests with brackets
  stat_compare_means(data = countdata %>% 
                       filter(DrugName == "Dasatinib") %>% 
                       filter(Time %in% c(6)) %>%
                       mutate(Condition = factor(Condition, levels = c("Combination Treatment", 
                                                                       "1 Week of Vemurafenib", 
                                                                       "2 Weeks of Vemurafenib", 
                                                                       "3 Weeks of Vemurafenib", 
                                                                       "Vemurafenib Only"))),
                     aes(x = Condition, y = Count),
                     method = "t.test", 
                     comparisons = timing_comparisons,
                     size = 6,
                     label = "p.format",
                     step.increase = 0.125,
                     tip.length = 0.01)

ggsave(file.path(figure_path,"lowdas_timingexperiment.pdf"), width = 24, height = 6.4, units = "in")

##