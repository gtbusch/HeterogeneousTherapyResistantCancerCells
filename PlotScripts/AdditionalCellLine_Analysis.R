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

current_path <- here()
parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "AdditionalCellLine")

figure_path <- file.path(parent_path, "Plots")

h358_data_long <- readRDS(file.path(new_path, "h358_sotorasib_adaptation_confluence.rds"))

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

indexmapping = c("A1" = "Control", "A2" = "Low Dose", "A3" = "JNKi", "A4" = "Low Dose + JNKi",
                "B1" = "High Dose", "B2" = "Low to High Dose at 7 Days", "B3" = "High Dose + JNKi", "B4" = "Low to High Dose + JNKi at 7 Days",
                "C1" = "Low to High Dose at 14 Days", "C2" = "Low to High Dose at 21 Days", "C3" = "Low to High Dose + JNKi at 14 Days", "C4" = "Low to High Dose + JNKi at 21 Days")

h358_data_long = h358_data_long %>% mutate(Time = round(as.numeric(Time)))

h358_data_long = h358_data_long %>% mutate(Time = as.numeric(Time)/24)

h358_data_long = h358_data_long %>% mutate(Time = ceiling(as.numeric(Time)))

##

## Task 1 - Plot the data frame

h358_data_indexed <- h358_data_long %>%
  mutate(Condition = indexmapping[PlateIndex])

ggplot((h358_data_indexed %>% filter(Condition %in% c("Control", "Low Dose", "High Dose")) %>% filter(Time %in% c(0, 7, 14, 21, 28))), 
       aes(x = factor(Time), y = as.numeric(Confluency))) +
  stat_summary(fun = mean, geom = "bar", aes(fill = Condition), 
               position = position_dodge(width = 0.8), width = 0.7, alpha = 0.3) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "Confluency at Different Doses of Sotorasib", 
       x = "Time (Days)",
       y = "Confluency",
       color = "Condition") +
  custom_figure_theme() +
  theme(legend.background = element_blank())

ggsave(file.path(figure_path,"highlow_sotorasib_experiment.pdf"), width = 12.8, height = 9.6, units = "in")

##