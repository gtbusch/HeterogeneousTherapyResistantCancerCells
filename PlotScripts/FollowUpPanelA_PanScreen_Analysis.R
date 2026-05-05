library(tidyverse)
library(here)
library(scales)
library(RColorBrewer)
library(circlize)
library(viridis)
library(gridExtra)
library(magrittr)
library(Cairo)
library(VennDiagram)
library(grid)
library(gridExtra)
library(plotly)
library(cowplot)
library(purrr)
library(fs)
library(ComplexHeatmap)
library(nsprcomp)
library(ggrepel)
library(gtools)
library(ggplotify)
library(ggbreak)
library(vegan)
library(writexl)
library(gtable)

##

current_path <- here()

parent_path <- dirname(current_path)

figure_path <- file.path(parent_path, "Plots")

## Task 0a - Load the Screens 1 & 2 Data----

new_path <- file.path(parent_path, "ExtractedData", "HighThroughputScreens")

panscreen <- readRDS(file.path(new_path, "panscreen.rds"))
panscreen_long <- readRDS(file.path(new_path, "panscreen_long.rds"))
panscreen_controls_long <- readRDS(file.path(new_path, "panscreen_controls_long.rds"))

custom_colors <- c("Parental" = "#5E5E5E", "D13" = "#FF7E79", "E2" = "#4294F8", "C12" = "#662D91", "A15" = "#3E8D27", "E2 + D13" = "#8E1356", "C14" = "#953877", "D5" = "#6ADCC1", "D4" = "#9E6D15",
                   "A2" = "#EF8BF9", "A11" = "#CC5C76", "C3" = "#F57946", "E8" = "#F9AD2A", "F9" = "#1D457F", "I11" = "#625A94", "G12" = "#0013F0", "E2 + D13" = "#39FF14")

custom_names <- c("Parental" = "Drug-Naive", "E2" = "Resistant Clone 1", "D13" = "Resistant Clone 2", "D4" = "Resistant Clone 3", "C12" = "Resistant Clone 4", "C14" = "Resistant Clone 5",
                  "D5" = "Resistant Clone 6", "F9" = "Resistant Clone 7", "G12" = "Resistant Clone 8", "I11" = "Resistant Clone 9", "A2" = "Resistant Clone 10", "A11" = "Resistant Clone 11",
                  "C3" = "Resistant Clone 12", "E8" = "Resistant Clone 13", "E2 + D13" = "Resistant Clones 1 and 2")

custom_cellline_pairs <- c(
  "D13 - E2" = "Resistant Clones 2 and 1",
  "D13 - D4" = "Resistant Clones 2 and 3",
  "D13 - C12" = "Resistant Clones 2 and 4",
  "E2 - D4" = "Resistant Clones 1 and 3",
  "E2 - C12" = "Resistant Clones 1 and 4",
  "D4 - C12" = "Resistant Clones 3 and 4"
)

custom_cellline_pairs_pcadist <- c(
  "D13 - E2" = "Clones 2 and 1",
  "D13 - D4" = "Clones 2 and 3",
  "C12 - D13" = "Clones 2 and 4",
  "D4 - E2" = "Clones 1 and 3",
  "C12 - E2" = "Clones 1 and 4",
  "C12 - D4" = "Clones 3 and 4"
)

custom_doses <- c(
  "0.01" = "0.01 \u00B5M",
  "0.1" = "0.1 \u00B5M",
  "1" = "1 \u00B5M",
  "10" = "10 \u00B5M"
)

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

#Clean up the pathway annotations
panscreen_long <- panscreen_long %>%
  mutate(Pathway = case_when(
    is.na(Pathway) ~ "Other Category",
    Pathway == "Others" | Pathway == "others" | Pathway == "Unknown" ~ "Other Category",
    Pathway == "DNA Damage" | Pathway == "DNA damage" ~ "DNA Damage/DNA Repair",
    Pathway == "NF-κB" | Pathway == "NF-?B" ~ "NF-\u03BAB",
    Pathway == "GPCR" ~ "GPCR & G Protein",
    Pathway == "Stem Cells &  Wnt" ~ "Stem Cells & Wnt",
    Pathway == "Endocrinology" ~ "Endocrinology & Hormones",
    Pathway == "Immunology" ~ "Immunology & Inflammation",
    Pathway == "antiangiogenic" ~ "Antiangiogenic",
    Pathway == "Angiogenesis∩┐╜&∩┐╜Epigenetics" ~ "Angiogenesis & Epigenetics",
    Pathway == "Tyrosine-Kinase" | Pathway == "Protein Tyrosine Kinase" ~ "Protein Tyrosine Kinase",
    Pathway == "Inflammation/ImmunologyCancer" ~ "Inflammation & Immunology & Cancer",
    TRUE ~ Pathway  # Keep original if no match
  ))

pathway_palette <- c(
  brewer.pal(8, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Paired"),
  brewer.pal(8, "Accent"),
  # Add a few more colors to reach 44
  "#FF34B3", "#2A52BE", "#32CD32", "#8B4513"  # Custom colors
)

pathway_color_mapping <- setNames(pathway_palette, unique(panscreen_long$Pathway))

# Define the HEX code you want to use
new_hex_code <- "#CCCCCC"  # This is gray - replace with your desired color
# Add the new named entry to your color mapping
pathway_color_mapping["Non-Top Pathways"] <- new_hex_code

create_readable_label <- function(cell_line) {
  # Extract base name (everything before underscore)
  base_name <- sub("_.*$", "", cell_line)
  
  # Map the base name using custom_names
  base_label <- custom_names[base_name]
  if(is.na(base_label)) base_label <- base_name  # fallback if no mapping exists
  
  # Combine the two parts
  paste(base_label)
}

##

## Task 0b - Load the In-Lab Validation Data ----

new_path <- file.path(parent_path, "ExtractedData", "HTSValidation", "CellTiterGloData")

lab_validation_data <- readRDS(file.path(new_path, "validationdata_20240403.rds"))

revision_round1_validation_data <- readRDS(file.path(new_path, "revision_round1data.rds"))

revision_round2_validation_data <- readRDS(file.path(new_path, "revision_round2data.rds"))

revision_round3_validation_data <- readRDS(file.path(new_path, "revision_round3data.rds"))

jq1_coculture_data <- readRDS(file.path(new_path, "jq1_coculture_data.rds"))


##

## Task 0c - Load the Data from (follow up panel a) ----

new_path <- file.path(parent_path, "ExtractedData", "FollowUpPanelA")

followup_panela_long <- readRDS(file.path(new_path, "followup_panela_long.rds"))

controls_followup_panela_long <- readRDS(file.path(new_path, "controls_followup_panela_long.rds"))

##

## Task 0d - Load the bulk rna seq data ----

new_path <- file.path(parent_path, "ExtractedData", "BulkRNASeq_Run1")

#list all subfolders in the parent folder
subfolders <- list.dirs(new_path, full.names = TRUE, recursive = TRUE)

#initialize a vector to store the paths of library RDS files
tpm_rds_files <- c()

#loop through each subfolder and find the library plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "bulkrnaseq_run1_tpm_plusmetadata\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  tpm_rds_files <- c(tpm_rds_files, files_in_subfolder)
}

# Load the file if found
if (length(tpm_rds_files) > 0) {
  run1_tpmdata <- readRDS(tpm_rds_files[1])
  print("File loaded successfully")
} else {
  print("No files found")
}

new_path <- file.path(parent_path, "ExtractedData", "BulkRNASeq_Run2")

#list all subfolders in the parent folder
subfolders <- list.dirs(new_path, full.names = TRUE, recursive = TRUE)

#initialize a vector to store the paths of library RDS files
tpm_rds_files <- c()

#loop through each subfolder and find the library plate RDS files
for (subfolder in subfolders) {
  # List RDS files in the subfolder
  files_in_subfolder <- list.files(subfolder, pattern = "bulkrnaseq_run2_tpm_plusmetadata\\.rds$", full.names = TRUE)
  
  #combine the files found with the main list
  tpm_rds_files <- c(tpm_rds_files, files_in_subfolder)
}

# Load the file if found
if (length(tpm_rds_files) > 0) {
  run2_tpmdata <- readRDS(tpm_rds_files[1])
  print("File loaded successfully")
} else {
  print("No files found")
}

##

## Task 0e - Process data into forms that will be useful for multiple tasks/figures ----

dosenorm_followup_panela_long = followup_panela_long %>%
  # Convert Dose to numeric and round to 5 decimal places
  mutate(Dose = round(as.numeric(Dose), digits = 5))

# Apply the logarithmic transformation, rounding, and re-exponentiation directly
dosenorm_followup_panela_long <- dosenorm_followup_panela_long %>%
  mutate(Dose = 10^round(log10(Dose), 2))

# Optional: Check the unique binned doses to see the effect of binning
unique_doses <- unique(dosenorm_followup_panela_long$Dose)
print(unique_doses)

#Rename the Drug column in the control follow up panel a dataset to match the taxonomy used for everything else
controls_followup_panela_long = controls_followup_panela_long %>% rename(Name = Drug)

#Rename the Drug column in the lab validatation dataset to match the taxonomy used for everything else
lab_validation_data = lab_validation_data %>% rename(Name = Drug)

revision_round1_validation_data = revision_round1_validation_data %>% rename(Name = Drug)

revision_round2_validation_data = revision_round2_validation_data %>% rename(Name = Drug)

##

## Task 0f - Plot the data for a single drug with all conditions, and in combination with vemurafenib, from the pan screen data only ----

#Figure 1 plot

selected_drug <- c("Dasatinib")

selected_data <- panscreen_long %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Treatment == "Vemurafenib") %>% filter(CellLine != "Parental") %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
        aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(title = paste("Response to 1 µM", selected_drug, "\nin High-Throughput Screens"),
       x = "",
       y = "% Viability") +
 custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names)

ggsave(file.path(figure_path, "dasatinibviability_1uM_hts_resistantlinesonly.pdf"), width = 7, height = 7.68, units = "in")

# Additional plots

selected_drug <- c("Deguelin")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-25, 105)) +
  labs(title = paste("Response to 1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names)+
  scale_y_continuous(breaks = c(-25, 0, 25, 50, 75, 100))

ggsave(file.path(figure_path, "deguelin_sharedkiller_hts_allcelllines.pdf"), width = 9.6, height = 9.6, units = "in")

selected_drug <- c("IACS-010759 (IACS-10759)")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-25, 105)) +
  labs(title = paste("Response to 1 µM IACS-010759"),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names)+
  scale_y_continuous(breaks = c(-25, 0, 25, 50, 75, 100))

ggsave(file.path(figure_path, "iacs010759_sharedkiller_hts_allcelllines.pdf"), width = 9.6, height = 9.6, units = "in")

selected_drug <- c("SR18662")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-25, 105)) +
  labs(title = paste("Response to 1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names)+
  scale_y_continuous(breaks = c(-25, 0, 25, 50, 75, 100))

ggsave(file.path(figure_path, "sr18662_uniquekiller_hts_allcelllines.pdf"), width = 9.6, height = 9.6, units = "in")

selected_drug <- c("SB273005")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 0.1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-25, 105)) +
  labs(title = paste("Response to 0.1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names)+
  scale_y_continuous(breaks = c(-25, 0, 25, 50, 75, 100))

ggsave(file.path(figure_path, "sb273005_sharedkiller_hts_allcelllines.pdf"), width = 9.6, height = 9.6, units = "in")

selected_drug <- c("S3I-201")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = paste("Response to 1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "s3i201dosecurve_hts_allcelllines.pdf"), width = 7.36, height = 6.4, units = "in")

selected_drug <- c("Dasatinib")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == .1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = paste("Response to 0.1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "dasatinibdosecurve_hts_allcelllines.pdf"), width = 7.36, height = 6.4, units = "in")

#Supplemental Figure 2 Plots

selected_drug <- c("S3I-201")

selected_data = panscreen_controls_long %>% pivot_longer(
  cols = c(Toxicity_Rep1, Toxicity_Rep2),
  names_to = "Replicate",
  values_to = "Toxicity") %>%
  mutate(Replicate = str_replace(Replicate, "Toxicity_Rep", "")) %>% 
  filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib")) %>%
  rename(Name = Drug)

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = Dose, y = (100 - Toxicity), color = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add points at each mean value
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Keep the error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(0, 125)) +
  labs(title = paste("Response to", selected_drug),
       x = "Dose",
       y = "% Viability",
       color = "Cell Line") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),
    legend.position = c(0.7,0.3)) +
  scale_x_log10(breaks = as.numeric(names(custom_doses)), labels = custom_doses) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path,"s3i_alldoses_htscontrols_allcelllines.pdf"), width = 7.04, height = 6.144, units = "in")

selected_drug <- c("Epirubicin HCl")

selected_data = panscreen_controls_long %>% pivot_longer(
  cols = c(Toxicity_Rep1, Toxicity_Rep2),
  names_to = "Replicate",
  values_to = "Toxicity") %>%
  mutate(Replicate = str_replace(Replicate, "Toxicity_Rep", "")) %>% 
  filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib")) %>%
  rename(Name = Drug)

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = Dose, y = (100 - Toxicity), color = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add points at each mean value
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Keep the error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(0, 125)) +
  labs(title = paste("Response to", selected_drug),
       x = "Dose",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_log10(breaks = as.numeric(names(custom_doses)), labels = custom_doses) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "epirubicinhcl_alldoses_htscontrols_allcelllines.pdf"), width = 7.04, height = 6.144, units = "in")

selected_drug <- c("Clofarabine")

selected_data = panscreen_controls_long %>% pivot_longer(
  cols = c(Toxicity_Rep1, Toxicity_Rep2),
  names_to = "Replicate",
  values_to = "Toxicity") %>%
  mutate(Replicate = str_replace(Replicate, "Toxicity_Rep", "")) %>% 
  filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib")) %>%
  rename(Name = Drug)

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = Dose, y = (100 - Toxicity), color = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add points at each mean value
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Keep the error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(0, 125)) +
  labs(title = paste("Response to", selected_drug),
       x = "Dose",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_log10(breaks = as.numeric(names(custom_doses)), labels = custom_doses) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "clofarabine_alldoses_htscontrols_allcelllines.pdf"), width = 7.04, height = 6.144, units = "in")

selected_drug <- c("Dasatinib")

selected_data = panscreen_controls_long %>% pivot_longer(
  cols = c(Toxicity_Rep1, Toxicity_Rep2),
  names_to = "Replicate",
  values_to = "Toxicity") %>%
  mutate(Replicate = str_replace(Replicate, "Toxicity_Rep", "")) %>% 
  filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib")) %>%
  rename(Name = Drug)

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = Dose, y = (100 - Toxicity), color = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  # Add points at each mean value
  stat_summary(fun = mean, geom = "point", size = 3) +
  # Keep the error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(0, 125)) +
  labs(title = paste("Response to", selected_drug),
       x = "Dose",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_log10(breaks = as.numeric(names(custom_doses)), labels = custom_doses) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "dasatinib_alldoses_htscontrols_allcelllines.pdf"), width = 7.04, height = 6.144, units = "in")

selected_drug <- c("Saracatinib (AZD0530)")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 10) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = "Response to 10 µM Saracatinib",
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "saracatinib10uM_hts_allcelllines.pdf"), width = 7.36, height = 6.4, units = "in")

selected_drug <- c("Saracatinib (AZD0530)")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = "Response to 1 µM Saracatinib",
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "saracatinib1uM_hts_allcelllines.pdf"), width = 7.36, height = 6.4, units = "in")

selected_drug <- c("(+)-JQ1")

selected_data = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(Dose == 0.1) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = paste("Response to 0.1 µM", selected_drug),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "jq101uM_hts_allcelllines.pdf"), width = 7.36, height = 6.4, units = "in")


##

## Task 0g - plot the data for a single drug from the follow up panel A data  ----

selected_drugs <- c("Deguelin", "SB273005")

selected_data = followup_panela_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment != "DMSO"))

selected_data <- selected_data %>% filter(Name %in% selected_drugs)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% 
          filter((Name == "Deguelin" & round(Dose, digits = 2) == 1.0) | 
                   (Name == "SB273005" & round(Dose, digits = 2) == 0.100)) %>% 
          mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = Name)) +
  # Plot bars grouped by drug
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.8) +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.9), width = 0.4, linewidth = 1) +
  # Use manually defined colors for the drugs
  scale_fill_manual(values = c("Deguelin" = "#02BFC4", "SB273005" = "#F8766D"), 
                    labels = c("Deguelin" = "Deguelin (PI3K/Akt inhib.)", "SB273005" = "SB273005 (Integrin inhib.)"),
                    name = "") +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = "Response to 1 µM Deguelin or \n0.1 µM SB273005 in the Follow-Up Panel",
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Angled text for better readability
        legend.position = c(0.7, 1.0),
        legend.background = element_blank()) +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "deguelin_1uM_sb273005_0_1uM_followuppanel_allcelllines.pdf"), width = 10, height = 7.68, units = "in")

selected_drugs <- c("IWP-O1", "SB273005")

selected_data = followup_panela_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment != "DMSO"))

selected_data <- selected_data %>% filter(Name %in% selected_drugs)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% 
          filter((Name == "IWP-O1" & round(Dose, digits = 2) == 3.16) | 
                   (Name == "SB273005" & round(Dose, digits = 2) == 0.100)) %>% 
          mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = Name)) +
  # Plot bars grouped by drug
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.8) +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(width = 0.9), width = 0.4, linewidth = 1) +
  # Use manually defined colors for the drugs
  scale_fill_manual(values = c("IWP-O1" = "#02BFC4", "SB273005" = "#F8766D"), 
                    labels = c("IWP-O1" = "IWP-O1 (Porcn inhib.)", "SB273005" = "SB273005 (Integrin inhib.)"),
                    name = "") +
  coord_cartesian(ylim = c(-5, 125)) +
  labs(title = "Response to 3.16 µM IWP-O1 or \n0.1 µM SB273005 in the Follow-Up Panel",
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Angled text for better readability
        legend.position = c(0.7, 1.0),
        legend.background = element_blank()) +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "iwpo1_316uM_sb273005_1uM_followuppanel_allcelllines.pdf"), width = 10, height = 7.68, units = "in")

selected_drug <- c("YM155 (Sepantronium Bromide)")

selected_data = followup_panela_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment != "DMSO"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(round(Dose, digits = 3) == 0.01) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 130)) +
  labs(title = paste("Response to 0.01 µM YM155 \nin the Follow-Up Panel"),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "ym155_01uM_followuppanel_allcelllines.pdf"), width = 10, height = 6.144, units = "in")

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(round(Dose, digits = 3) == 0.003) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 130)) +
  labs(title = paste("Response to 0.003 µM YM155 \nin the Follow-Up Panel"),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "ym155_003uM_followuppanel_allcelllines.pdf"), width = 10, height = 4.45, units = "in")

selected_drug <- c("GMX1778 (CHS828)")

selected_data = followup_panela_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment != "DMSO"))

selected_data <- selected_data %>% filter(Name %in% selected_drug)

selected_data$Dose <- as.numeric(selected_data$Dose)
selected_data$Treatment <- as.factor(selected_data$Treatment)

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(round(Dose, digits = 3) == 0.01) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 130)) +
  labs(title = paste("Response to 0.01 µM GMX1778 \nin the Follow-Up Panel"),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "gmx1778_01uM_followuppanel_allcelllines.pdf"), width = 10, height = 6.144, units = "in")

#plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
ggplot((selected_data %>% filter(round(Dose, digits = 3) == 0.003) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
       aes(x = CellLine, y = (100 - Toxicity), fill = CellLine)) +
  # Plot lines for the average of the points
  stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  coord_cartesian(ylim = c(-5, 130)) +
  labs(title = paste("Response to 0.003 µM GMX1778 \nin the Follow-Up Panel"),
       x = "",
       y = "% Viability") +
  custom_figure_theme() + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") +
  scale_x_discrete(labels = custom_names) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))

ggsave(file.path(figure_path, "gmx1778_003uM_followuppanel_allcelllines.pdf"), width = 10, height = 4.45, units = "in")

##

## Task 0h - Process data for Stephanie Huang's group ----

controls_followup_panela_longer = controls_followup_panela_long %>% pivot_longer(
  cols = starts_with("Toxicity_Rep"),
  names_to = "Replicate",
  values_to = "Toxicity"
) %>%
  mutate(
    Replicate = as.numeric(gsub("Toxicity_Rep", "", Replicate))
  )

controls_followup_panela_longer = controls_followup_panela_longer %>% mutate(Screen = "A")

followup_panela_data = rbind(dosenorm_followup_panela_long, controls_followup_panela_longer)

drug_info <- panscreen_long %>%
  select(Name, Target, Pathway) %>%
  distinct()

# Merge into your long_data
followup_panela_data <- followup_panela_data %>%
  left_join(drug_info, by = "Name")

write_tsv(followup_panela_data, "followup_panela_data.tsv")
write_tsv(panscreen_long, "highthroughput_data.tsv")

AUC_results <- panscreen_long %>%
  mutate(LogDose = log10(Dose)) %>%
  group_by(Name, Target, Pathway, CellLine, Treatment) %>%
  arrange(desc(LogDose)) %>%
  summarise(AUC = sum((lag(Toxicity, default = first(Toxicity)) + Toxicity) * 
                        (lag(LogDose, default = first(LogDose)) - LogDose) / 2, na.rm = TRUE))

write_tsv(AUC_results, "AUC_data_panscreen.tsv")

AUC_panela_results <- followup_panela_data %>%
  mutate(LogDose = log10(Dose)) %>%
  group_by(Name, CellLine, Treatment) %>%
  arrange(desc(LogDose)) %>%
  summarise(AUC = sum((lag(Toxicity, default = first(Toxicity)) + Toxicity) * 
                        (lag(LogDose, default = first(LogDose)) - LogDose) / 2, na.rm = TRUE))

write_tsv(AUC_panela_results, "AUC_data_followuppanela.tsv")

##

## Task 2 - Plot the toxicity for each drug/dose pair in Panel A, Replicate 1 vs Replicate 2 ----

task2_dosenorm_followup_panela_long <- dosenorm_followup_panela_long %>% select(-Screen) %>%
  pivot_wider(
    names_from = Replicate,   # Create columns from Replicate #
    values_from = Toxicity)   # Fill cells with toxicity values

task2_dosenorm_followup_panela_long$Toxicity_Rep1 = task2_dosenorm_followup_panela_long$"1"
task2_dosenorm_followup_panela_long$Toxicity_Rep2 = task2_dosenorm_followup_panela_long$"2"

task2_dosenorm_followup_panela_long = task2_dosenorm_followup_panela_long %>% select(-c("1","2"))

task2_all_followup_panela_long = rbind(task2_dosenorm_followup_panela_long, controls_followup_panela_long)

#Plot the toxicity for each dose/pair in Rep 1 vs Rep 2 in Panel A, colored by drug
ggplot(task2_all_followup_panela_long, aes(x = (100 - Toxicity_Rep1), y = (100 - Toxicity_Rep2), color = as.factor(Name))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +  # Adds the y=x line 
  facet_wrap(~Name, nrow = 3, ncol = 8, labeller = labeller(Name = function(x) {
    # Create a named vector just for the drugs you want to rename
    drug_renames <- c(
      "YM155 (Sepantronium Bromide)" = "YM155",
      "Delanzomib (CEP-18770)" = "Delanzomib",
      "Saracatinib (AZD0530)" = "Saracatinib",
      "GMX1778 (CHS828)" = "GMX1778",
      "Carfilzomib (PR-171)" = "Carfilzomib"
      
      # Add only the drugs you want to rename
    )
    # Replace only specified drugs, leave others unchanged
    ifelse(x %in% names(drug_renames), drug_renames[x], x)
  })) +
  labs(title = "Comparison of Viability Across Replicates in Panel A",
       x = "% Viability of Replicate 1",
       y = "% Viability of Replicate 2",
       color = "Name") +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-5, 150), ylim = c(-5, 150)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

ggsave(file.path(figure_path, "scatter_followuppanela_replicates.pdf"), width = 20.2, height = 10, units = "in")

## Task 3 - Create a heatmap of the data from the initial screen showing the panel developed for Follow Up Panel A and the corresponding follow up panel A data, as well as a distance matrix ----

drug_dose_mapping <- data.frame(
  Name = sort(unique(followup_panela_long$Name)), 
  DesiredDose = c(1, 1, .01, 1, .01, 1, 1, .01, 1, 1, 1, 1, 1, 1, 10, .1, 1, 10, .1, .01)  
)

heatmap_colors = colorRamp2(seq(0, 125, length.out = 125), viridis(125))

# generate a heatmap of the hts responses to the drugs at the doses of interest that were used in follow-up panel a 

#filter the panscreen data to only include cell line/treatment pairs used in the validation screen A
task3_panscreen_long <- panscreen_long %>%
  filter((CellLine == "Parental" & Treatment == "DMSO") | 
           (CellLine != "Parental" & Treatment == "Vemurafenib")
  )

#filter the pan screen data to only include drugs used in the validation screen A and remove the Target and Pathway columns
task3_panscreen_long <- task3_panscreen_long %>%
  filter(Name %in% followup_panela_long$Name) %>% select(-Target) %>% select(-Pathway)

#group Pan Screen Data by Name, CellLine, Treatment, and Dose then calculating the average Toxicity for each group
task3_panscreen_long_avg <- task3_panscreen_long %>%
  group_by(Name, CellLine, Treatment, Dose) %>%
  summarise(PanScreen_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

task3_panscreen_long_avg = task3_panscreen_long_avg %>% mutate(PanScreen_Tox = (100 - PanScreen_Tox))

task3_panscreen_heatmap_data <- task3_panscreen_long_avg %>%
  # Join with the mapping dataframe
  inner_join(drug_dose_mapping, by = "Name") %>%
  # Filter rows where the Dose matches the DesiredDose from the mapping
  filter(Dose == DesiredDose) %>%
  # Optionally, remove the DesiredDose column if it's no longer needed
  select(-DesiredDose)

task3_panscreen_heatmap_data <- task3_panscreen_heatmap_data %>% select(-Treatment) %>%
  pivot_wider(
    names_from = CellLine,   # Create columns from CellLine
    values_from = PanScreen_Tox,  # Fill the cells with toxicity values
    values_fill = list(PanScreen_Tox = NA)  # Fill missing values with NA
  ) %>% select(-Dose)

# Step 2: Convert the resulting dataframe to a matrix
task3_panscreen_toxicity_matrix <- as.matrix(task3_panscreen_heatmap_data[, -1])  # Exclude the first column (Name) for the matrix
rownames(task3_panscreen_toxicity_matrix) <- task3_panscreen_heatmap_data$Name  # Set row names as drug names
colnames(task3_panscreen_toxicity_matrix) <- c("Resistant Clone 4", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 1", "Drug-Naive")

# Now you can use the toxicity_matrix with ComplexHeatmap
drug_order = c("VLX1570", "Delanzomib (CEP-18770)", "Saracatinib (AZD0530)","EC330","Autophinib","IWP-O1","SR18662","Deguelin","MCB-613","NSC228155","YM155 (Sepantronium Bromide)",
               "Carfilzomib (PR-171)", "IMD 0354", "JH-RE-06", "Fluvastatin Sodium", "TPX-0005", "GMX1778 (CHS828)", "SB273005", "(+)-JQ1", "OTX015")

cell_order = c("Drug-Naive", "Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")
# Reorder the rows of the matrix
task3_reordered_matrix <- task3_panscreen_toxicity_matrix[match(drug_order, rownames(task3_panscreen_toxicity_matrix)), ]
task3_reordered_matrix <- task3_reordered_matrix[ ,match(cell_order, colnames(task3_panscreen_toxicity_matrix))]

rownames(task3_reordered_matrix)[rownames(task3_reordered_matrix) == "YM155 (Sepantronium Bromide)"] <- "YM155"
rownames(task3_reordered_matrix)[rownames(task3_reordered_matrix) == "Carfilzomib (PR-171)"] <- "Carfilzomib"

heatmap1 <- Heatmap(task3_reordered_matrix, 
                    name = "% Viability \nin HTS", 
                    col = heatmap_colors, 
                    row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

# Convert to a grob (graphical object)
heatmap1_grob <- grid.grabExpr(draw(heatmap1))

# Convert to a ggplot object
heatmap1_plot <- ggplotify::as.ggplot(heatmap1_grob)

ggsave(file.path(figure_path, "hts_heatmap_dosesofinterest_paneladrugs.pdf"), plot = heatmap1_plot, width = 12, height = 8, units = "in")

#

# generate a heatmap of the follow-up panel a responses to the drugs at the doses of interest for the cell lines used in HTS

#filter the validation screen A data to remove the Screen column
task3_followup_panela_long = dplyr::select(followup_panela_long, -Screen)

#filter the validation screen A data to only include cell lines that were used in the main screens
task3_followup_panela_long = task3_followup_panela_long %>% filter(CellLine %in% c("Parental","D13","D4","C12","E2"))

#group Validation Screen A Data by Name, CellLine, Treatment, and Dose then calculating the average Toxicity for each group
task3_followup_panela_long_avg <- task3_followup_panela_long %>%
  group_by(Name, CellLine, Treatment, Dose) %>%
  summarise(PanelA_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

task3_followup_panela_long_avg = task3_followup_panela_long_avg %>% mutate(PanelA_Tox = (100 - PanelA_Tox))

#Filter only the doses used in screens 1 and 2
task3_followup_panela_long_avg <- task3_followup_panela_long_avg %>%
  mutate(Dose = round(as.numeric(Dose), digits = 2)) %>%
  # Filter for specific dose values
  filter(Dose %in% c(10, 1, 0.1, 0.01))

task3_followup_panela_heatmap_data <- task3_followup_panela_long_avg %>%
  # Join with the mapping dataframe
  inner_join(drug_dose_mapping, by = "Name") %>%
  # Filter rows where the Dose matches the DesiredDose from the mapping
  filter(Dose == DesiredDose) %>%
  # Optionally, remove the DesiredDose column if it's no longer needed
  select(-DesiredDose)

task3_panela_heatmap_data <- task3_followup_panela_heatmap_data %>% select(-Treatment) %>%
  pivot_wider(
    names_from = CellLine,   # Create columns from CellLine
    values_from = PanelA_Tox,  # Fill the cells with toxicity values
    values_fill = list(PanelA_Tox = NA)  # Fill missing values with NA
  ) %>% select(-Dose)

# Step 2: Convert the resulting dataframe to a matrix
task3_panela_toxicity_matrix <- as.matrix(task3_panela_heatmap_data[, -1])  # Exclude the first column (Name) for the matrix
rownames(task3_panela_toxicity_matrix) <- task3_panela_heatmap_data$Name  # Set row names as drug names
colnames(task3_panela_toxicity_matrix) <- c("Resistant Clone 4", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 1", "Drug-Naive")

# Now you can use the toxicity_matrix with ComplexHeatmap
drug_order = c("VLX1570", "Delanzomib (CEP-18770)", "Saracatinib (AZD0530)","EC330","Autophinib","IWP-O1","SR18662","Deguelin","MCB-613","NSC228155","YM155 (Sepantronium Bromide)",
               "Carfilzomib (PR-171)", "IMD 0354", "JH-RE-06", "Fluvastatin Sodium", "TPX-0005", "GMX1778 (CHS828)", "SB273005", "(+)-JQ1", "OTX015")

cell_order = c("Drug-Naive", "Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")
# Reorder the rows of the matrix
task3_panela_reordered_matrix <- task3_panela_toxicity_matrix[match(drug_order, rownames(task3_panela_toxicity_matrix)), ]
task3_panela_reordered_matrix <- task3_panela_reordered_matrix[ ,match(cell_order, colnames(task3_panela_toxicity_matrix))]

rownames(task3_panela_reordered_matrix)[rownames(task3_panela_reordered_matrix) == "YM155 (Sepantronium Bromide)"] <- "YM155"
rownames(task3_panela_reordered_matrix)[rownames(task3_panela_reordered_matrix) == "Carfilzomib (PR-171)"] <- "Carfilzomib"

heatmap2 <- Heatmap(task3_panela_reordered_matrix, 
                    name = "% Viability \nin Panel A", 
                    col = heatmap_colors, 
                    row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

# Convert to a grob (graphical object)
heatmap2_grob <- grid.grabExpr(draw(heatmap2, padding = unit(c(3, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmap2_plot <- ggplotify::as.ggplot(heatmap2_grob)

ggsave(file.path(figure_path, "panela_heatmap_dosesofinterest_paneladrugs.pdf"), plot = heatmap2_plot, width = 12, height = 8, units = "in")

#

#generate a heatmap of the follow-up panel a responses to the drugs that validated at the dose of interest for all cell lines 

drug_dose_mapping_panela <- data.frame(
  Name = sort(unique(followup_panela_long$Name)), 
  DesiredDose = c(1, 1, .01, 1, .01, 1, 1, .01, 1, 3.16, 1, 1, 1, 1, 10, .1, 1, 10, .1, .01) # The specific dose values for each drug, change out 6th position for EC330 
)

screena_heatmap_data <- dosenorm_followup_panela_long %>%
  # Join with the mapping dataframe
  inner_join(drug_dose_mapping_panela, by = "Name") %>% 
  # Filter rows where the Dose matches the DesiredDose from the mapping
  filter(abs(Dose - DesiredDose) < 0.005) %>%
  # Optionally, remove the DesiredDose column if it's no longer needed
  select(-DesiredDose) %>% select(-Dose) %>% select(-Screen)

controls_followup_panela_heatmap <- controls_followup_panela_long %>%
  pivot_longer(
    cols = starts_with("Toxicity_Rep"),
    names_to = c("Replicate"),
    names_pattern = "Toxicity_Rep(\\d)",
    values_to = "Toxicity"
  ) %>%
  mutate(Replicate = as.integer(Replicate)) # Ensure Replicate is an integer

controls_followup_panela_heatmap = controls_followup_panela_heatmap %>% filter(Name == "Dasatinib" & Dose == 0.1) %>% select(-Dose)

screena_heatmap_data = rbind(screena_heatmap_data, controls_followup_panela_heatmap)

screena_heatmap_data_avg <- screena_heatmap_data %>%
  group_by(Name, CellLine, Treatment) %>%
  summarise(PanelA_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

screena_heatmap_data_avg = screena_heatmap_data_avg %>% mutate(PanelA_Tox = (100 - PanelA_Tox))

screena_heatmap_data_avg <- screena_heatmap_data_avg %>% select(-Treatment) %>%
  pivot_wider(
    names_from = CellLine,   # Create columns from CellLine
    values_from = PanelA_Tox,  # Fill the cells with toxicity values
    values_fill = list(PanelA_Tox = NA)  # Fill missing values with NA
  )

# Step 2: Convert the resulting dataframe to a matrix
screena_heatmap_data_avg = screena_heatmap_data_avg %>% filter(Name %in% c("Dasatinib", "(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                           "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570"))

screena_toxicity_matrix <- as.matrix(screena_heatmap_data_avg[, -1])  # Exclude the first column (Name) for the matrix
rownames(screena_toxicity_matrix) <- screena_heatmap_data_avg$Name  # Set row names as drug names

# Reorder the matrix columns
reordered_screena_matrix <- screena_toxicity_matrix[, names(custom_names[1:14])]

# Create new labels for the columns
new_labels <- custom_names[names(custom_names[1:14])]

wrapped_labels <- sapply(new_labels, function(label) {
  if(grepl("Resistant Clone", label)) {
    gsub("Resistant Clone", "Resistant\nClone", label)
  } else if(grepl("Drug-Naive", label)) {
    "Drug-Naive"
  } else {
    # Default handling for other label formats
    gsub(" ", "\n", label)
  }
})

rownames(reordered_screena_matrix)[rownames(reordered_screena_matrix) == "YM155 (Sepantronium Bromide)"] <- "YM155"

heatmap <- Heatmap(reordered_screena_matrix, 
                   name = "% Viability", 
                   col = heatmap_colors, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE,
                   column_labels = wrapped_labels,
                   column_names_rot = 0,
                   column_names_centered = TRUE,
                   # Match your theme's font sizes
                   row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                   column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                   column_names_max_height = unit(6, "cm"),
                   bottom_annotation = HeatmapAnnotation(
                     foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
                   ),
                   # Legend text size
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                     labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

print(heatmap)

# Convert to a grob (graphical object)
heatmap_grob <- grid.grabExpr(draw(heatmap, padding = unit(c(3, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmap_plot <- ggplotify::as.ggplot(heatmap_grob)

ggsave(file.path(figure_path, "followuppanela_heatmap_dosesofinterest_validateddrugs.pdf"), plot = heatmap_plot, width = 24, height = 6, units = "in")

#make a version of the heatmap with only the drugs we used to predict combination sensitivity

reordered_screena_matrix_subset = reordered_screena_matrix[c(2:3,5:6,8:10,13:14,17),]

binned_heatmap_colors <- colorRamp2(c(0, 29.99, 30, 69.99, 70, 120), 
                             c("#C74444", "#C74444", "#E6E68A", "#E6E68A", "#6B9E6B", "#6B9E6B"))

heatmap <- Heatmap(reordered_screena_matrix_subset, 
                   name = "% Viability", 
                   col = binned_heatmap_colors, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE,
                   column_labels = wrapped_labels,
                   column_names_rot = 0,
                   column_names_centered = TRUE,
                   # Match your theme's font sizes
                   row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                   column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                   column_names_max_height = unit(6, "cm"),
                   bottom_annotation = HeatmapAnnotation(
                     foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
                   ),
                   # Legend text size
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                     labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

print(heatmap)

# Convert to a grob (graphical object)
heatmap_grob <- grid.grabExpr(draw(heatmap, padding = unit(c(3, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmap_plot <- ggplotify::as.ggplot(heatmap_grob)

ggsave(file.path(figure_path, "followuppanela_heatmap_dosesofinterest_binnedcombinationdrugs.pdf"), plot = heatmap_plot, width = 24, height = 6, units = "in")

reordered_screena_matrix_subset = reordered_screena_matrix_subset[,-1]

# Function to categorize and count
count_effects_pct <- function(viabilities) {
  n_total <- length(viabilities[!is.na(viabilities)])
  killed <- sum(viabilities < 30, na.rm = TRUE)
  partial <- sum(viabilities >= 30 & viabilities < 70, na.rm = TRUE)
  
  killed_pct <- round((killed / n_total) * 100, 1)
  partial_pct <- round((partial / n_total) * 100, 1)
  
  return(paste0(killed_pct, "% / ", partial_pct, "%"))  # Format: "killed% / partial%"
}

# Get drug names
drugs <- rownames(reordered_screena_matrix_subset)
n_drugs <- length(drugs)

# Create matrix to store results
result_matrix <- matrix("", nrow = n_drugs, ncol = n_drugs)
rownames(result_matrix) <- drugs
colnames(result_matrix) <- drugs

# Fill diagonal (single drug effects)
for (i in 1:n_drugs) {
  result_matrix[i, i] <- count_effects_pct(reordered_screena_matrix_subset[i, ])
}

# Fill lower triangle (drug combinations)
for (i in 2:n_drugs) {
  for (j in 1:(i-1)) {
    # Take minimum viability (best case scenario)
    combined_viability <- pmin(reordered_screena_matrix_subset[i, ], 
                               reordered_screena_matrix_subset[j, ])
    result_matrix[i, j] <- count_effects_pct(combined_viability)
  }
}

# Convert to data frame for plotting
result_df <- as.data.frame(result_matrix)

# Create theme with larger fonts
mytheme <- ttheme_default(
  base_size = 14,
  base_family = "Helvetica",
  core = list(
    fg_params = list(hjust = 0.5, vjust = 0.5, fontface = "plain"),
    bg_params = list(fill = c("white", "grey95"))
  ),
  colhead = list(
    fg_params = list(hjust = 0.5, vjust = 0.5, fontface = "bold"),
    bg_params = list(fill = "grey90")
  ),
  rowhead = list(
    fg_params = list(hjust = 1, vjust = 0.5, fontface = "bold"),
    bg_params = list(fill = "grey90")
  )
)

# Create and save
pdf(file.path(figure_path,"drug_combination_table.pdf"), width = 20, height = 20)
grid.table(result_df, theme = mytheme)
dev.off()

#Extract the cell line toxicity profiles (transpose matrix so cell lines are rows)
cell_line_profiles <- as.data.frame(t(reordered_screena_matrix))

#Calculate the distance matrix using Euclidean distance
dist_matrix <- dist(cell_line_profiles, method = "euclidean")

#Convert to a regular matrix for easier access/visualization
dist_matrix <- as.matrix(dist_matrix)

#Set row and column names to be the cell line names
rownames(dist_matrix) <- colnames(reordered_screena_matrix)
colnames(dist_matrix) <- colnames(reordered_screena_matrix)

# Step 6: Create the ComplexHeatmap
heatmap2 <- Heatmap(
  dist_matrix,
  name = "Distance",
  col = viridis(300),
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  column_labels = wrapped_labels,
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_labels = wrapped_labels,
  row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
  column_title_side = "top",
  row_title_side = "left",
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  bottom_annotation = HeatmapAnnotation(
    foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
  ),
  # Legend text size
  heatmap_legend_param = list(title = "Viability \nPercentage Points \nEuclidean Distance",
    title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

# Step 7: Draw the heatmap
print(heatmap2)

# Convert to a grob (graphical object)
heatmap2_grob <- grid.grabExpr(draw(heatmap2, padding = unit(c(3, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmap2_plot <- ggplotify::as.ggplot(heatmap2_grob)

ggsave(file.path(figure_path, "followuppanela_distancematrix.pdf"), plot = heatmap2_plot, width = 24, height = 9, units = "in")

# Calculate Spearman correlation matrix
# Transpose the data so cell lines are columns
spearman_cor <- cor(t(cell_line_profiles), method = "spearman")

#Convert to a regular matrix for easier access/visualization
spearman_cor <- as.matrix(spearman_cor)

#Set row and column names to be the cell line names
rownames(spearman_cor) <- colnames(reordered_screena_matrix)
colnames(spearman_cor) <- colnames(reordered_screena_matrix)

# Step 6: Create the ComplexHeatmap
heatmap3 <- Heatmap(
  spearman_cor,
  name = "Distance",
  col = viridis(300),
  clustering_distance_rows = "spearman",
  clustering_distance_columns = "spearman",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  column_labels = wrapped_labels,
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_labels = wrapped_labels,
  row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
  column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
  column_title_side = "top",
  row_title_side = "left",
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  bottom_annotation = HeatmapAnnotation(
    foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
  ),
  # Legend text size
  heatmap_legend_param = list(title = "Spearman\nCorrelation",
                              title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

# Step 7: Draw the heatmap
print(heatmap3)

# Convert to a grob (graphical object)
heatmap3_grob <- grid.grabExpr(draw(heatmap3, padding = unit(c(3, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmap3_plot <- ggplotify::as.ggplot(heatmap3_grob)

ggsave(file.path(figure_path, "followuppanela_spearmancorrelation.pdf"), plot = heatmap3_plot, width = 24, height = 9, units = "in")


#

# generate a heatmap of the follow-up panel a responses to the drugs that validated at an individual dose for all cell lines 

# Function to create a heatmap for a specific dose
create_dose_heatmap <- function(dose_value, all_drugs, cell_lines, custom_names, show_row_names = TRUE) {
  # Create drug-dose mapping for this specific dose
  drug_bydose_mapping <- data.frame(
    Name = unique(followup_panela_long$Name), 
    DesiredDose = rep(dose_value, times = length(unique(followup_panela_long$Name)))
  )
  
  # Filter data for this dose
  dose_heatmap_data <- dosenorm_followup_panela_long %>%
    # Join with the mapping dataframe
    inner_join(drug_bydose_mapping, by = "Name") %>% 
    # Filter rows where the Dose matches the DesiredDose from the mapping
    filter(Dose == DesiredDose) %>%
    # Remove unnecessary columns
    select(-DesiredDose, -Dose)
  
  # Check if we have any data for this dose
  if(nrow(dose_heatmap_data) == 0) {
    message(paste0("No data available for dose: ", dose_value))
    return(NULL)
  }
  
  # Calculate averages
  dose_heatmap_data_avg <- dose_heatmap_data %>%
    group_by(Name, CellLine, Treatment) %>%
    summarise(PanelA_Tox = mean(Toxicity, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>%
    select(-Treatment) %>%
    mutate(PanelA_Tox = (100 - PanelA_Tox)) %>%
    pivot_wider(
      names_from = CellLine,
      values_from = PanelA_Tox,
      values_fill = list(PanelA_Tox = NA)
    )
  
  # Filter for validated drugs - but keep this as a separate filtering step
  validated_drugs <- c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015", 
                       "SB273005", "SR18662", "TPX-0005", "YM155 (Sepantronium Bromide)",
                       "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", 
                       "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570")
  
  # Create a complete matrix with all drugs and all cell lines, filled with NA
  complete_matrix <- matrix(NA, 
                            nrow = length(all_drugs), 
                            ncol = length(cell_lines),
                            dimnames = list(all_drugs, cell_lines))
  
  # Fill in the values we have from the data
  if(nrow(dose_heatmap_data_avg) > 0) {
    for(drug in intersect(dose_heatmap_data_avg$Name, all_drugs)) {
      drug_row <- dose_heatmap_data_avg[dose_heatmap_data_avg$Name == drug, ]
      for(cell in intersect(names(drug_row)[-1], cell_lines)) {
        complete_matrix[drug, cell] <- drug_row[[cell]]
      }
    }
  }
  
  # Apply the reordering to columns 
  reordered_matrix <- complete_matrix[, cell_lines, drop = FALSE]
  
  # Create column labels
  new_labels <- custom_names[colnames(reordered_matrix)]
  
  # Shorten YM155 name if present
  if("YM155 (Sepantronium Bromide)" %in% rownames(reordered_matrix)) {
    rownames(reordered_matrix)[rownames(reordered_matrix) == "YM155 (Sepantronium Bromide)"] <- "YM155"
  }
  
  if("Carfilzomib (PR-171)" %in% rownames(reordered_matrix)) {
    rownames(reordered_matrix)[rownames(reordered_matrix) == "Carfilzomib (PR-171)"] <- "Carfilzomib"
  }
  
  # Shorten YM155 name if present
  if("GMX1778 (CHS828)" %in% rownames(reordered_matrix)) {
    rownames(reordered_matrix)[rownames(reordered_matrix) == "GMX1778 (CHS828)"] <- "GMX"
  }
  
  # Shorten YM155 name if present
  if("Fluvastatin Sodium" %in% rownames(reordered_matrix)) {
    rownames(reordered_matrix)[rownames(reordered_matrix) == "Fluvastatin Sodium"] <- "Fluv. Sodium"
  }
  
  # Create the heatmap with proper color mapping
  hm <- Heatmap(reordered_matrix, 
                name = paste0("Dose: ", dose_value), 
                col = heatmap_colors, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE,
                column_labels = new_labels,
                column_names_rot = 90,
                column_names_centered = TRUE,
                row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                column_names_max_height = unit(4, "cm"),
                # Configure row names display based on parameter
                show_row_names = show_row_names,
                # Show heatmap title
                column_title = paste0("Dose: ", round(dose_value, digits = 5), "\u00B5M"),
                column_title_gp = gpar(fontsize = 25.6, fontface = "bold"),
                # For the grid layout, we'll handle the legend separately
                show_heatmap_legend = FALSE)
  
  return(hm)
}

# Define the doses you want to visualize - use exact values from your data
doses <- unique(dosenorm_followup_panela_long$Dose)

doses <- doses[1:8]


# Define the validated drug order (same as in your original code)
drug_order <- c("VLX1570", "EC330", "Autophinib", "IWP-O1", "SR18662", "Deguelin", 
                "NSC228155", "YM155 (Sepantronium Bromide)", "Carfilzomib (PR-171)",
                "JH-RE-06", "Fluvastatin Sodium", "TPX-0005", "GMX1778 (CHS828)", 
                "SB273005", "(+)-JQ1", "OTX015")

# Get the cell line order based on custom_names
cell_line_order <- names(custom_names)

# Generate heatmaps for each dose and keep only non-NULL results
heatmap_list <- list()
for(i in seq_along(doses)) {
  # Calculate position in the grid
  show_rows <- (i == length(doses))
  
  # Create the heatmap with appropriate row names setting
  hm <- create_dose_heatmap(doses[i], drug_order, cell_line_order, custom_names, show_row_names = show_rows)
  
  if(!is.null(hm)) {
    heatmap_list <- c(heatmap_list, list(hm))
  }
}

# Check if we have any heatmaps
if(length(heatmap_list) == 0) {
  stop("No heatmaps could be created for any of the doses. Please check your data.")
}

# Create a common legend for all heatmaps
lgd <- Legend(col_fun = heatmap_colors, 
              title = "% Viability", 
              at = seq(0, 150, 25),  
              title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
              labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"))

# Create the grid of heatmaps
ht_list <- NULL
for(i in 1:length(heatmap_list)) {
  # Add to the list
  if(i == 1) {
    ht_list <- heatmap_list[[i]]
  } else {
    ht_list <- ht_list + heatmap_list[[i]]
  }
}

# Convert to a grob (graphical object)
heatmaps_grob <- grid.grabExpr(draw(ht_list, padding = unit(c(20, 2, 2, 2), "mm")))

# Convert to a ggplot object
heatmaps_plot <- ggplotify::as.ggplot(heatmaps_grob)

ggsave(file.path(figure_path, "followuppanela_heatmap_alldoses_validateddrugs.pdf"), plot = heatmaps_plot, width = 24, height = 10, units = "in")

legend_grob <- grid.grabExpr(draw(lgd))

legend_plot <- ggplotify::as.ggplot(legend_grob)

ggsave(file.path(figure_path, "followuppanela_legend_alldoses_validateddrugs.pdf"), plot = legend_plot, width = 2.56, height = 2.4, units = "in")

##

## Task 5 - Run PCA on the transcriptomic data from the bulk rna sequencing runs ----

#5a - Run PCA on the first bulk rna seq run, NOT averaged by replicate 

    task5a_data = run1_tpmdata %>% select(-sample) %>% filter(cell_line != "E9")
    
    #combine the cell line and replicate columns
    task5a_data <- task5a_data %>%
      unite("CellLine_Replicate", cell_line, replicate, sep = "_", remove = TRUE)
    
    task5a_data <- t(task5a_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5a_data) = as.character(unlist(task5a_data[1, ]))
    
    task5a_data = task5a_data[-1, ]
    
    task5a_data = t(task5a_data)
    
    task5a_data <- task5a_data[!rownames(task5a_data) %in% c("D11_Late2", "A15_Early"), ]
    
    task5a_data_numeric <- apply(task5a_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5a_data_numeric, 2, var) == 0
    
    #remove these columns
    task5a_data_numeric <- task5a_data_numeric[, !zero_var_cols]
    
    task5a_pca_result <- prcomp(task5a_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5a_pca_result))
    
    task5a_pca_data <- as.data.frame(task5a_pca_result$x)
    
    task5a_pca_scores <- as.data.frame(task5a_pca_result$x[,])
    task5a_pca_scores$cell_line <- rownames(task5a_data)
    
    task5a_pca_scores_mapped <- task5a_pca_scores
    
    # Apply the mapping to your data
    task5a_pca_scores_mapped$readable_label <- sapply(task5a_pca_scores$cell_line, create_readable_label)
    
    ggplot(task5a_pca_scores_mapped %>% filter(readable_label != ""), aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = readable_label),
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf,
                      size = 6) +
      labs(color = "Row Number",
           x = "Principal Component 1 - 8.419%",
           y = "Principal Component 2 - 6.470%",
           title = "PCA of the Transcriptomic Data in Run 1") +
      custom_figure_theme() +
      coord_cartesian(xlim = c(-50, 250), ylim = c(-200, 50))
    
    
    ggsave(file.path(figure_path, "pc1and2_transcriptomicdata_notaveraged.pdf"), width = 9.6, height = 9.6, units = "in")

#5b - Run PCA on the second bulk rna seq run, NOT averaged by replicate

    task5b_data = run2_tpmdata %>% select(-sample)
    
    #combine the cell line and replicate columns
    task5b_data <- task5b_data %>%
      unite("CellLine_Replicate", cell_line, replicate, sep = "_", remove = TRUE)
    
    task5b_data <- t(task5b_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5b_data) = as.character(unlist(task5b_data[1, ]))
    
    task5b_data = task5b_data[-1, ]
    
    task5b_data = t(task5b_data)
    
    task5b_data_numeric <- apply(task5b_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5b_data_numeric, 2, var) == 0
    
    #remove these columns
    task5b_data_numeric <- task5b_data_numeric[, !zero_var_cols]
    
    task5b_pca_result <- prcomp(task5b_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5b_pca_result))
    
    task5b_pca_data <- as.data.frame(task5b_pca_result$x)
    
    ggplot(task5b_pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      xlab("Principal Component 1") +
      ylab("Principal Component 2") +
      custom_theme()
    
    task5b_pca_scores <- as.data.frame(task5b_pca_result$x[,])
    task5b_pca_scores$cell_line <- rownames(task5b_data)
    
    task5b_pca_scores_mapped <- task5b_pca_scores
    
    # Apply the mapping to your data
    task5b_pca_scores_mapped$readable_label <- sapply(task5b_pca_scores$cell_line, create_readable_label)
    
    ggplot(task5b_pca_scores_mapped, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = readable_label), 
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf) +
      labs(color = "Row Number",
           x = "Principal Component 1 - 16.01%",
           y = "Principal Component 2 - 12.54%",
           title = "PCA Plot") +
      custom_theme()
    
    ggplot(task5b_pca_scores_mapped, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = case_when(
        grepl("Naive", readable_label) ~ "naive",
        grepl("#1[0-3]", readable_label) ~ "high_lines",
        TRUE ~ "other_lines"
      )), size = 3) + 
      geom_text_repel(aes(label = readable_label,
                          color = case_when(
                            grepl("Naive", readable_label) ~ "naive",
                            grepl("#1[0-3]", readable_label) ~ "high_lines",
                            TRUE ~ "other_lines"
                          )), 
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf) +
      scale_color_manual(values = c("high_lines" = "red", 
                                    "other_lines" = "black",
                                    "naive" = "#5E5E5E"),
                         breaks = c("high_lines", "other_lines", "naive"),
                         labels = c("Dabrafenib/Trametinib-Resistant", 
                                    "Vemurafenib-Resistant",
                                    "Drug-Naive")) +
      labs(color = "Primary Resistance",
           x = "Principal Component 1 - 16.01%",
           y = "Principal Component 2 - 12.54%",
           title = "PCA Plot") +
      custom_theme()
    
    ggplot(task5b_pca_scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), 
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf) +
      labs(color = "Row Number",
           x = "Principal Component 1 - 16.01%",
           y = "Principal Component 2 - 12.54%",
           title = "PCA Plot") +
      custom_theme()
    
#5c - Run PCA on the first bulk rna seq run, averaged by replicate
    
    task5c_data = run1_tpmdata %>% 
                  select(-sample) %>% 
                  filter(cell_line != "E9") %>% 
                  filter(!(cell_line == "A15" & replicate == "Early")) %>%
                  filter(!(cell_line == "D11" & replicate == "Late2"))
    
    #average by replicates
    task5c_data<- task5c_data %>%
      group_by(cell_line) %>%
      summarise(across(A1BG:last_col(), \(x) mean(x, na.rm = TRUE)))
    
    task5c_data <- t(task5c_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5c_data) = as.character(unlist(task5c_data[1, ]))
    
    task5c_data = task5c_data[-1, ]
    
    task5c_data = t(task5c_data)
    
    task5c_data_numeric <- apply(task5c_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5c_data_numeric, 2, var) == 0
    
    #remove these columns
    task5c_data_numeric <- task5c_data_numeric[, !zero_var_cols]
    
    task5c_pca_result <- prcomp(task5c_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5c_pca_result))
    
    task5c_pca_data <- as.data.frame(task5c_pca_result$x)
    
    task5c_pca_scores <- as.data.frame(task5c_pca_result$x[,])
    task5c_pca_scores$cell_line <- rownames(task5c_data)

#5d - Run PCA on the second bulk rna seq run, averaged by replicate
   
    task5d_data = run2_tpmdata %>% select(-sample)
    
    #average by replicates
    task5d_data<- task5d_data %>%
      group_by(cell_line) %>%
      summarise(across(A1BG:last_col(), \(x) mean(x, na.rm = TRUE)))
    
    task5d_data <- t(task5d_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5d_data) = as.character(unlist(task5d_data[1, ]))
    
    task5d_data = task5d_data[-1, ]
    
    task5d_data = t(task5d_data)
    
    task5d_data_numeric <- apply(task5d_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5d_data_numeric, 2, var) == 0
    
    #remove these columns
    task5d_data_numeric <- task5d_data_numeric[, !zero_var_cols]
    
    task5d_pca_result <- prcomp(task5d_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5d_pca_result))
    
    task5d_pca_data <- as.data.frame(task5d_pca_result$x)
    
    task5d_pca_scores <- as.data.frame(task5d_pca_result$x[,])
    task5d_pca_scores$cell_line <- rownames(task5d_data)
    
# 5e - Run PCA on the combined bulk rna seq runs, averaged by replicate
    
    #pull data from the first bulk RNA seq run
    task5e_run1_data = run1_tpmdata %>% select(-sample) %>% 
      filter(cell_line != "E9") %>%
      filter(!(cell_line == "A15" & replicate == "Early")) %>%
      filter(!(cell_line == "D11" & replicate == "Late2"))
    
    #pull data from the second bulk RNA seq run
    task5e_run2_data = run2_tpmdata %>% select(-sample)
    
    # Find common columns
    task5e_common_columns <- intersect(colnames(task5e_run1_data), colnames(task5e_run2_data))
    
    # Subset both data frames
    task5e_run1_subset <- task5e_run1_data[, task5e_common_columns]
    task5e_run2_subset <- task5e_run2_data[, task5e_common_columns]
    
    # Combine the subsetted data frames
    task5e_combined_data <- rbind(task5e_run1_subset, task5e_run2_subset)
    
    #average by replicates
    task5e_combined_data<- task5e_combined_data %>%
      group_by(cell_line) %>%
      summarise(across(A1BG:last_col(), \(x) mean(x, na.rm = TRUE)))
    
    task5e_data <- t(task5e_combined_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5e_data) = as.character(unlist(task5e_data[1, ]))
    
    task5e_data = task5e_data[-1, ]
    
    task5e_data = t(task5e_data)
    
    task5e_data_numeric <- apply(task5e_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5e_data_numeric, 2, var) == 0
    
    #remove these columns
    task5e_data_numeric <- task5e_data_numeric[, !zero_var_cols]
    
    task5e_pca_result <- prcomp(task5e_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5e_pca_result))
    
    task5e_pca_data <- as.data.frame(task5e_pca_result$x)
    
    task5e_pca_scores <- as.data.frame(task5e_pca_result$x[,])
    task5e_pca_scores$cell_line <- rownames(task5e_data)
    
    ggplot(task5e_pca_scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), 
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf) +
      labs(color = "Row Number",
           x = "Principal Component 1 - 14.140%",
           y = "Principal Component 2 - 7.875%",
           title = "PCA Plot") +
      custom_figure_theme()
    
    #Extract the cell line toxicity profiles (transpose matrix so cell lines are rows)
    cell_line_profiles_5e <- as.data.frame(task5e_data[c(1, 3, 8, 9, 11, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27), ])
    
    cell_line_profiles_5e_pca <- as.data.frame(task5e_pca_scores[,-1])
    cell_line_profiles_5e_pca <- cell_line_profiles_5e_pca %>% filter(cell_line %in% c("Parental", "Naive", "E2", "D13", "D4", "C12", "C14", "D5", "F9", "G12", "I11", "A2", "A11", "C3", "E8"))
    rownames(cell_line_profiles_5e_pca) <- cell_line_profiles_5e_pca$cell_line
    cell_line_profiles_5e_pca <- cell_line_profiles_5e_pca %>% select(-cell_line)
    
    #Calculate the distance matrix using Euclidean distance
    dist_matrix_5e <- dist(cell_line_profiles_5e, method = "euclidean")
    
    dist_matrix_5e_pca <- dist(cell_line_profiles_5e_pca, method = "euclidean")
    
    #Convert to a regular matrix for easier access/visualization
    dist_matrix_5e <- as.matrix(dist_matrix_5e)
    
    dist_matrix_5e_pca <- as.matrix(dist_matrix_5e_pca)
    
    #Set row and column names to be the cell line names
    rownames(dist_matrix_5e) <- rownames(task5e_data[c(1, 3, 8, 9, 11, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27), ])
    colnames(dist_matrix_5e) <- rownames(task5e_data[c(1, 3, 8, 9, 11, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27), ])
    
    rownames(dist_matrix_5e_pca) <- rownames(cell_line_profiles_5e_pca)
    colnames(dist_matrix_5e_pca) <- rownames(cell_line_profiles_5e_pca)
    
    
    custom_names_pca <- c("Parental" = "Drug-Naive Run 2", "Naive" = "Drug-Naive Run 1", "E2" = "Resistant Clone 1", "D13" = "Resistant Clone 2", "D4" = "Resistant Clone 3", "C12" = "Resistant Clone 4", "C14" = "Resistant Clone 5",
                      "D5" = "Resistant Clone 6", "F9" = "Resistant Clone 7", "G12" = "Resistant Clone 8", "I11" = "Resistant Clone 9", "A2" = "Resistant Clone 10", "A11" = "Resistant Clone 11",
                      "C3" = "Resistant Clone 12", "E8" = "Resistant Clone 13")
    
    # Get the current column/row names from your distance matrix
    current_names <- colnames(dist_matrix_5e_pca)
    
    # Create new labels based on the current names
    new_labels_pca <- custom_names_pca[current_names]
    
    wrapped_labels_pca <- sapply(new_labels_pca, function(label) {
      if(grepl("Resistant Clone", label)) {
        gsub("Resistant Clone", "Resistant\nClone", label)
      } else if(grepl("Drug-Naive Run 1", label)) {
        "Naive\nRun 1"
      } else if(grepl("Drug-Naive Run 2", label)) {
        "Naive\nRun 2"
      } else {
        # Default handling for other label formats
        gsub(" ", "\n", label)
      }
    })
    
    # Step 6: Create the ComplexHeatmap
    heatmap5e <- Heatmap(
      dist_matrix_5e,
      name = "Distance",
      col = viridis(300),
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      column_names_rot = 0,
      column_names_centered = TRUE,
      row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
      column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
      column_title_side = "top",
      row_title_side = "left",
      row_dend_width = unit(2, "cm"),
      column_dend_height = unit(2, "cm"),
      bottom_annotation = HeatmapAnnotation(
        foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
      ),
      # Legend text size
      heatmap_legend_param = list(title = "Transcriptomic \nEuclidean Distance",
                                  title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                                  labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))
    
    # Step 7: Draw the heatmap
    print(heatmap5e)
    
    # Step 6: Create the ComplexHeatmap
    heatmap5e_pca <- Heatmap(
      dist_matrix_5e_pca,
      name = "Distance",
      col = viridis(300),
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      column_names_rot = 0,
      column_names_centered = TRUE,
      column_labels = wrapped_labels_pca,
      row_labels = wrapped_labels_pca,
      row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
      column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
      column_title_side = "top",
      row_title_side = "left",
      row_dend_width = unit(2, "cm"),
      column_dend_height = unit(2, "cm"),
      bottom_annotation = HeatmapAnnotation(
        foo = anno_empty(border = FALSE, height = unit(0.5, "cm"))
      ),
      # Legend text size
      heatmap_legend_param = list(title = "PCA \nEuclidean \nDistance",
                                  title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                                  labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))
    
    # Step 7: Draw the heatmap
    print(heatmap5e_pca)
    
    # Convert to a grob (graphical object)
    heatmap5e_pca_grob <- grid.grabExpr(draw(heatmap5e_pca, padding = unit(c(3, 2, 2, 2), "mm")))
    
    # Convert to a ggplot object
    heatmap5e_pca_plot <- ggplotify::as.ggplot(heatmap5e_pca_grob)
    
    ggsave(file.path(figure_path, "transcriptomicpca_distancematrix.pdf"), plot = heatmap5e_pca_plot, width = 24, height = 9, units = "in")
    
# 5f - Run PCA on the combined bulk rna seq runs, NOT averaged by replicate
    
    #pull data from the first bulk RNA seq run
    task5f_run1_data = run1_tpmdata %>% select(-sample) %>% filter(cell_line != "E9")
    
    #combine the cell line and replicate columns
    task5f_run1_data <- task5f_run1_data %>%
      unite("sample.name", cell_line, replicate, sep = "_", remove = TRUE)
    
    #pull data from the second bulk RNA seq run
    task5f_run2_data = run2_tpmdata %>% select(-sample)
    
    #combine the cell line and replicate columns
    task5f_run2_data <- task5f_run2_data %>%
      unite("sample.name", cell_line, replicate, sep = "_", remove = TRUE)
    
    # Find common columns
    task5f_common_columns <- intersect(colnames(task5f_run1_data), colnames(task5f_run2_data))
    
    # Subset both data frames
    task5f_run1_subset <- task5f_run1_data[, task5f_common_columns]
    task5f_run2_subset <- task5f_run2_data[, task5f_common_columns]
    
    # Combine the subsetted data frames
    task5f_combined_data <- rbind(task5f_run1_subset, task5f_run2_subset)
    
    # Filter out outlier samples
    
    task5f_combined_data <- task5f_combined_data %>% filter(sample.name != "A15_Early") %>% filter(sample.name != "D11_Late2") 
    
    task5f_data <- t(task5f_combined_data)
    
    #use the first row which contains the cell lines as the column names of the new matrix
    colnames(task5f_data) = as.character(unlist(task5f_data[1, ]))
    
    task5f_data = task5f_data[-1, ]
    
    task5f_data = t(task5f_data)
    
    task5f_data_numeric <- apply(task5f_data, 2, as.numeric)
    
    #identify columns with zero variance
    zero_var_cols <- apply(task5f_data_numeric, 2, var) == 0
    
    #remove these columns
    task5f_data_numeric <- task5f_data_numeric[, !zero_var_cols]
    
    task5f_pca_result <- prcomp(task5f_data_numeric, center = TRUE, scale. = TRUE)
    print(summary(task5f_pca_result))
    
    task5f_pca_data <- as.data.frame(task5f_pca_result$x)
    
    ggplot(task5f_pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      xlab("Principal Component 1") +
      ylab("Principal Component 2") +
      custom_figure_theme()
    
    task5f_pca_scores <- as.data.frame(task5f_pca_result$x[, 1:4])
    task5f_pca_scores$cell_line <- rownames(task5f_data)
    
    ggplot(task5f_pca_scores, aes(x = PC3, y = PC4)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), 
                      box.padding = 0.5, 
                      point.padding = 0.5,
                      force = 2,
                      max.overlaps = Inf) +
      labs(color = "Row Number",
           x = "Principal Component 3 - 4.369%",
           y = "Principal Component 4 - 3.678%",
           title = "PCA Plot") +
      custom_figure_theme()
    
##

## Task 13 - PCA of screens 1 & 2 data, with 10 uM doses removed, averaged by replicates ----
    
    #remove the screen index column
    task13_line_by_drugpair = dplyr::select(panscreen_long, -c(Target, Pathway))
    
    #keep only the data of resistant colonies treated in combination with vemurafenib and the parental cells treated with DMSO only
    task13_line_by_drugpair = task13_line_by_drugpair %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                                                   (CellLine != "Parental" & Treatment != "DMSO"))
    
    task13_line_by_drugpair = task13_line_by_drugpair %>% filter(Dose != 10)
    
    #average replicates between screens (only parental has two)
    task13_line_by_drugpair <- task13_line_by_drugpair %>%
      group_by(Name, CellLine, Treatment, Dose) %>%
      summarise(PanScreen_Tox = mean(Toxicity, na.rm = TRUE)) %>%
      ungroup()
    
    #create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
    task13_line_by_drugpair_table = task13_line_by_drugpair %>%
      pivot_wider(
        names_from = c(Name, Dose),   # Create columns from CellLine
        values_from = PanScreen_Tox,  # Fill the cells with toxicity values
        values_fill = list(PanScreen_Tox = NA)  # Fill missing values with NA
      )
    
    #remove the treatment column
    task13_line_by_drugpair_table = dplyr::select(task13_line_by_drugpair_table, -Treatment)
    
    #transpose the table 
    task13_line_by_drugpair_table = t(task13_line_by_drugpair_table)
    
    #set the cell line row as the column names
    colnames(task13_line_by_drugpair_table) = as.character(unlist(task13_line_by_drugpair_table[1, ]))
    
    #remove the cell line row
    task13_line_by_drugpair_table = task13_line_by_drugpair_table[-1, ]
    
    #flip the table back
    task13_line_by_drugpair_table = t(task13_line_by_drugpair_table)
    
    #convert the table to numeric by column
    task13_line_by_drugpair_table_numeric = apply(task13_line_by_drugpair_table, 1, as.numeric)
    
    #flip the table back since the previous line transposed it
    task13_line_by_drugpair_table_numeric = t(task13_line_by_drugpair_table_numeric)
    
    #set the column names as the column names of the previous version of the table
    colnames(task13_line_by_drugpair_table_numeric) = colnames(task13_line_by_drugpair_table)
    
    #run pca
    task13_pca_result <- prcomp(task13_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)
    
    #print the results 
    print(summary(task13_pca_result))
    
    #store the values for each pc for each cell line as a data frame
    task13_pca_data <- as.data.frame(task13_pca_result$x)
    
    #store the loadings for each pc for each feature as a data frame
    task13_pca_loadings <- task13_pca_result$rotation
    
    #store pcs 1:4 in a new data frame
    task13_pca_scores <- as.data.frame(task13_pca_result$x[,])
    
    #add a column to the new dataframe with the cell line info
    task13_pca_scores$cell_line <- rownames(task13_line_by_drugpair_table_numeric)
    
    ##

    
    
## Task 10 - PCA based on drug toxicity data from Follow Up Panel A that validated & the control drugs in Follow Up Panel A, NOT averaged by replicates ----
    
    #Turn the toxicity_replicate columns into a toxicity column and a replicate column
    controls_followup_panela_task10 <- controls_followup_panela_long %>%
      pivot_longer(
        cols = starts_with("Toxicity_Rep"),
        names_to = c("Replicate"),
        names_pattern = "Toxicity_Rep(\\d)",
        values_to = "Toxicity"
      ) %>%
      mutate(Replicate = as.integer(Replicate)) # Ensure Replicate is an integer
    
    #remove the drugs that didn't replicate in this validation screen 
    full_followup_panela_long = followup_panela_long %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                            "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570"))
    
    #remove 10 uM doses
    full_followup_panela_long = full_followup_panela_long %>% filter(Dose != 10)
    
    #add the control panel a data to the main follow up panel a dataset
    full_followup_panela_task10 <- rbind((full_followup_panela_long %>% select(-Screen)), controls_followup_panela_task10)
    
    #create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
    full_line_by_drugpair_table = full_followup_panela_task10 %>%
      pivot_wider(
        names_from = c(Name, Dose),   # Create columns from CellLine
        values_from = Toxicity  # Fill the cells with toxicity values
      )
    
    #remove the treatment column
    full_line_by_drugpair_table = dplyr::select(full_line_by_drugpair_table, -Treatment)
    
    #combine the cell line and replicate columns
    full_line_by_drugpair_table <- full_line_by_drugpair_table %>%
      unite("CellLine_Replicate", CellLine, Replicate, sep = "_", remove = TRUE)
    
    #transpose the table 
    full_line_by_drugpair_table = t(full_line_by_drugpair_table)
    
    #set the cell line row as the column names
    colnames(full_line_by_drugpair_table) = as.character(unlist(full_line_by_drugpair_table[1, ]))
    
    #remove the cell line row
    full_line_by_drugpair_table = full_line_by_drugpair_table[-1, ]
    
    #flip the table back
    full_line_by_drugpair_table = t(full_line_by_drugpair_table)
    
    #convert the table to numeric by column
    full_line_by_drugpair_table_numeric = apply(full_line_by_drugpair_table, 1, as.numeric)
    
    #flip the table back since the previous line transposed it
    full_line_by_drugpair_table_numeric = t(full_line_by_drugpair_table_numeric)
    
    #set the column names as the column names of the previous version of the table
    colnames(full_line_by_drugpair_table_numeric) = colnames(full_line_by_drugpair_table)
    
    #run pca for the full table with validation screen a data and the corresponding control data
    full_pca_result <- prcomp(full_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)
    
    #print the results 
    print(summary(full_pca_result))
    
    #store the values for each pc for each cell line as a data frame
    full_pca_data <- as.data.frame(full_pca_result$x)
    
    #store the loadings for each pc for each feature as a data frame
    full_pca_loadings <- full_pca_result$rotation
    
    #plot the data for two pcs WITHOUT cell line labels on each point
    ggplot(full_pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      xlab("Principal Component 1") +
      ylab("Principal Component 2") +
      custom_theme()
    
    #store pcs 1:4 in a new data frame
    full_pca_scores <- as.data.frame(full_pca_result$x[, 1:4])
    
    #add a column to the new dataframe with the cell line info
    full_pca_scores$cell_line <- rownames(full_line_by_drugpair_table_numeric)
    
    #plot the data for two pcs WITH cell line labels on each point
    ggplot(full_pca_scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), max.overlaps = Inf) +
      custom_theme() +
      labs(color = "Row Number",
           x = "Principal Component 1",
           y = "Principal Component 2",
           title = "PCA of Follow-Up Panel A Data")
    
    ggplot(full_pca_scores, aes(x = PC3, y = PC4)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), max.overlaps = Inf) +
      custom_theme() +
      labs(color = "Row Number",
           x = "Principal Component 3",
           y = "Principal Component 4",
           title = "PCA of Follow-Up Panel A Data")
    
    ##
    
## Task 11 - PCA based on drug toxicity data from Follow Up Panel A & the control drugs in Follow Up Panel A, averaged by replicates ----
    
    #Turn the toxicity_replicate columns into a toxicity column and a replicate column
    avg_controls_followup_panela_task11 <- controls_followup_panela_long %>%
      pivot_longer(
        cols = starts_with("Toxicity_Rep"),
        names_to = c("Replicate"),
        names_pattern = "Toxicity_Rep(\\d)",
        values_to = "Toxicity"
      ) %>%
      mutate(Replicate = as.integer(Replicate)) # Ensure Replicate is an integer
    
    #remove the drugs that didn't replicate in this validation screen 
    avg_full_followup_panela_long = followup_panela_long %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                                "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570"))
    
    #add the control panel a data to the main follow up panel a dataset
    avg_full_followup_panela_task11 <- rbind((avg_full_followup_panela_long %>% select(-Screen)), avg_controls_followup_panela_task11)
    
    avg_full_followup_panela_task11$Dose <- as.numeric(avg_full_followup_panela_task11$Dose)
    
    avg_full_followup_panela_task11 <- avg_full_followup_panela_task11 %>% filter(Dose != 10)
    
    #average replicates within the follow up panel
    avg_full_line_by_drugpair <- avg_full_followup_panela_task11 %>%
      group_by(Name, CellLine, Treatment, Dose) %>%
      summarise(Toxicity = mean(Toxicity, na.rm = TRUE)) %>%
      ungroup()
    
    #create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
    avg_full_line_by_drugpair_table = avg_full_line_by_drugpair %>%
      pivot_wider(
        names_from = c(Name, Dose),   # Create columns from CellLine
        values_from = Toxicity  # Fill the cells with toxicity values
      )
    
    avg_full_line_by_drugpair_table <- avg_full_line_by_drugpair_table %>% select(CellLine, Treatment, '(+)-JQ1_1.001442', Autophinib_1.001442, 'Carfilzomib (PR-171)_0.010014',
                                                                                       Deguelin_1.001442, EC330_1, `Fluvastatin Sodium_1.001442`, `GMX1778 (CHS828)_0.010014`,
                                                                                       `IWP-O1_3.164557`, 'JH-RE-06_1.001442', NSC228155_1.001442, OTX015_1.001442, SB273005_1.001442,
                                                                                       SR18662_1.001442, `TPX-0005_3.164557`, VLX1570_0.100144, `YM155 (Sepantronium Bromide)_0.010029`)
    
    #remove the treatment column
    avg_full_line_by_drugpair_table = dplyr::select(avg_full_line_by_drugpair_table, -Treatment)
    
    #transpose the table 
    avg_full_line_by_drugpair_table = t(avg_full_line_by_drugpair_table)
    
    #set the cell line row as the column names
    colnames(avg_full_line_by_drugpair_table) = as.character(unlist(avg_full_line_by_drugpair_table[1, ]))
    
    #remove the cell line row
    avg_full_line_by_drugpair_table = avg_full_line_by_drugpair_table[-1, ]
    
    #flip the table back
    avg_full_line_by_drugpair_table = t(avg_full_line_by_drugpair_table)
    
    #convert the table to numeric by column
    avg_full_line_by_drugpair_table_numeric = apply(avg_full_line_by_drugpair_table, 1, as.numeric)
    
    #flip the table back since the previous line transposed it
    avg_full_line_by_drugpair_table_numeric = t(avg_full_line_by_drugpair_table_numeric)
    
    #set the column names as the column names of the previous version of the table
    colnames(avg_full_line_by_drugpair_table_numeric) = colnames(avg_full_line_by_drugpair_table)
    
    #run pca for the full table with validation screen a data and the corresponding control data
    avg_full_pca_result <- prcomp(avg_full_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)
    
    #print the results 
    print(summary(avg_full_pca_result))
    
    #store the values for each pc for each cell line as a data frame
    avg_full_pca_data <- as.data.frame(avg_full_pca_result$x)
    
    #store the loadings for each pc for each feature as a data frame
    avg_full_pca_loadings <- avg_full_pca_result$rotation
    
    #plot the data for two pcs WITHOUT cell line labels on each point
    ggplot(avg_full_pca_data, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      xlab("Principal Component 1") +
      ylab("Principal Component 2") +
      custom_figure_theme()
    
    #store pcs in a new data frame
    avg_full_pca_scores <- as.data.frame(avg_full_pca_result$x[,])
    
    #add a column to the new dataframe with the cell line info
    avg_full_pca_scores$cell_line <- rownames(avg_full_line_by_drugpair_table_numeric)
    
    #plot the data for two pcs WITH cell line labels on each point
    ggplot(avg_full_pca_scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) + 
      geom_text_repel(aes(label = cell_line), max.overlaps = Inf) +
      custom_figure_theme() +
      labs(color = "Row Number",
           x = "Principal Component 1",
           y = "Principal Component 2",
           title = "PCA Plot")
    
    ##
    
## Task 6 - Calculate and compare cell line distances in PCA space between transcriptomic and drug screening datasets ----
    
# Function to calculate pairwise Euclidean distances between cell lines in PCA space
  calculate_pca_distances <- function(pca_data, id_column = "cell_line") {
      # Extract cell line names
      cell_lines <- unique(pca_data[[id_column]])
      
      # Identify PCA component columns
      pca_cols <- grep("^PC", names(pca_data), value = TRUE)
      
      # Create a matrix for distances
      dist_matrix <- matrix(0, nrow = length(cell_lines), ncol = length(cell_lines))
      rownames(dist_matrix) <- cell_lines
      colnames(dist_matrix) <- cell_lines
      
      # Calculate Euclidean distances between all pairs of cell lines
      for (i in 1:length(cell_lines)) {
        for (j in 1:length(cell_lines)) {
          if (i != j) {
            # Extract PCA coordinates for cell line i
            cell_i_data <- pca_data[pca_data[[id_column]] == cell_lines[i], pca_cols]
            
            # Extract PCA coordinates for cell line j
            cell_j_data <- pca_data[pca_data[[id_column]] == cell_lines[j], pca_cols]
            
            # Calculate Euclidean distance
            dist_matrix[i, j] <- sqrt(sum((cell_i_data - cell_j_data)^2))
          }
        }
      }
      
      return(dist_matrix)
    }
    
# Function to compare distances between two datasets with statistical significance
  compare_distances <- function(dist_matrix1, dist_matrix2, labels = c("Transcriptomic", "Drug Screening")) {
      # Ensure both matrices have the same cell lines
      common_cell_lines <- intersect(rownames(dist_matrix1), rownames(dist_matrix2))
      
      if (length(common_cell_lines) < 2) {
        stop("Not enough common cell lines between datasets")
      }
      
      # Filter matrices to include only common cell lines
      dist_matrix1 <- dist_matrix1[common_cell_lines, common_cell_lines]
      dist_matrix2 <- dist_matrix2[common_cell_lines, common_cell_lines]
      
      # Convert distance matrices to vectors (excluding diagonal elements)
      dist_vec1 <- c()
      dist_vec2 <- c()
      pair_labels <- c()
      
      for (i in 1:(length(common_cell_lines) - 1)) {
        for (j in (i + 1):length(common_cell_lines)) {
          dist_vec1 <- c(dist_vec1, dist_matrix1[i, j])
          dist_vec2 <- c(dist_vec2, dist_matrix2[i, j])
          pair_labels <- c(pair_labels, paste(common_cell_lines[i], "-", common_cell_lines[j]))
        }
      }
      
      # Create a data frame for analysis and visualization
      comparison_df <- data.frame(
        pair = pair_labels,
        dist1 = dist_vec1,
        dist2 = dist_vec2,
        label1 = labels[1],
        label2 = labels[2]
      )
      
      # Calculate correlation between distances with p-values
      pearson_test <- cor.test(dist_vec1, dist_vec2, method = "pearson")
      spearman_test <- cor.test(dist_vec1, dist_vec2, method = "spearman")
      
      # For Mantel test, we need to ensure we're working with proper distance objects
      # Create full matrices with zeros on the diagonal
      n <- length(common_cell_lines)
      mat1 <- matrix(0, nrow = n, ncol = n)
      mat2 <- matrix(0, nrow = n, ncol = n)
      rownames(mat1) <- colnames(mat1) <- common_cell_lines
      rownames(mat2) <- colnames(mat2) <- common_cell_lines
      
      # Fill the matrices with the distance values
      index <- 1
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          mat1[i,j] <- mat1[j,i] <- dist_vec1[index]
          mat2[i,j] <- mat2[j,i] <- dist_vec2[index]
          index <- index + 1
        }
      }
      
      # Convert to dist objects
      dist1 <- as.dist(mat1)
      dist2 <- as.dist(mat2)
      
      # Perform Mantel test
      mantel_result <- mantel(dist1, dist2, method = "pearson", permutations = 999)
      
      
      pearson_cor <- pearson_test$estimate
      spearman_cor <- spearman_test$estimate
      pearson_p <- pearson_test$p.value
      spearman_p <- spearman_test$p.value
      mantel_cor <- mantel_result$statistic
      mantel_p <- mantel_result$signif
      
      # Determine significance based on conventional p < 0.05 threshold
      pearson_significant <- pearson_p < 0.05
      spearman_significant <- spearman_p < 0.05
      mantel_significant <- mantel_p < 0.05
      
      return(list(
        comparison_df = comparison_df,
        pearson_cor = pearson_cor,
        spearman_cor = spearman_cor,
        mantel_cor = mantel_cor,
        pearson_p = pearson_p,
        spearman_p = spearman_p,
        mantel_p = mantel_p,
        pearson_significant = pearson_significant,
        spearman_significant = spearman_significant,
        mantel_significant = mantel_significant
      ))
    }
  
  plot_distance_comparison <- function(comparison_results, 
                                       title = "Comparison of Cell Line \nDistances in PCA Space") {
    
    df <- comparison_results$comparison_df
    pearson <- comparison_results$pearson_cor
    spearman <- comparison_results$spearman_cor
    pearson_p <- comparison_results$pearson_p
    spearman_p <- comparison_results$spearman_p
    mantel <- comparison_results$mantel_cor
    mantel_p <- comparison_results$mantel_p
    
    # Format p-values for display
    pearson_p_text <- ifelse(pearson_p < 0.001, "p < 0.001", 
                             ifelse(pearson_p < 0.01, paste0("p = ", round(pearson_p, 3)),
                                    paste0("p = ", round(pearson_p, 3))))
    
    spearman_p_text <- ifelse(spearman_p < 0.001, "p < 0.001", 
                              ifelse(spearman_p < 0.01, paste0("p = ", round(spearman_p, 3)),
                                     paste0("p = ", round(spearman_p, 3))))
    
    mantel_p_text <- ifelse(mantel_p < 0.001, "p < 0.001", 
                              ifelse(mantel_p < 0.01, paste0("Mantel p = ", round(mantel_p, 3)),
                                     paste0("Mantel p = ", round(mantel_p, 3))))
    
    # Add significance asterisks
    pearson_sig <- ifelse(pearson_p < 0.001, "***", 
                          ifelse(pearson_p < 0.01, "**",
                                 ifelse(pearson_p < 0.05, "*", "")))
    
    spearman_sig <- ifelse(spearman_p < 0.001, "***", 
                           ifelse(spearman_p < 0.01, "**",
                                  ifelse(spearman_p < 0.05, "*", "")))
    
    mantel_sig <- ifelse(mantel_p < 0.001, "***", 
                           ifelse(mantel_p < 0.01, "**",
                                  ifelse(mantel_p < 0.05, "*", "")))
  
    #df$mapped_name <- custom_cellline_pairs_pcadist[df$pair]
    
    # Create scatterplot
    p <- ggplot(df, aes(x = dist1, y = dist2)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
      labs(
        title = title,
        subtitle = paste0("Pearson r = ", round(mantel, 3), " ", mantel_sig, " (", mantel_p_text, ")"),
        x = paste("Distance in", unique(df$label1), "PCA Space"),
        y = paste("Distance in", unique(df$label2), "PCA Space")
      ) +
     geom_text(data = df, label = df$pair) +
      custom_figure_theme()# +
    # coord_cartesian(xlim = c(160, 290), ylim = c(75, 150))
    
   # p <- p + geom_text(
   #    data = df,
   #    aes(label = df$mapped_name),
   #   nudge_y = max(df$dist2) * 0.01,  # Dynamic nudge based on data range
   #  size = 6,
   #   check_overlap = FALSE  # Ensure no labels are dropped
   # )
    
    # Calculate differences between distances
    df$diff <- abs(df$dist1 - df$dist2)
    
    # Get the most different pairs (up to 10, but handle case of fewer rows)
    n_diff <- min(10, nrow(df))
    most_different <- df[order(df$diff, decreasing = TRUE), ][1:n_diff, ]
    
    # Return the plot and data for most different pairs
    return(list(
      plot = p,
      most_different_pairs = most_different
    ))
  }
    
  
#Calculate distance matrices
trans_dist <- calculate_pca_distances(task5c_pca_scores)
drug_dist <- calculate_pca_distances(task13_pca_scores) 
  
#Compare and visualize
comparison <- compare_distances(trans_dist, drug_dist)
plot_results <- plot_distance_comparison(comparison)
print(plot_results$plot)

ggsave(file.path(figure_path, "pca_distancecomparisons.pdf"), width = 8, height = 8, units = "in")

trans_allclones_dist <- calculate_pca_distances(task5e_pca_scores %>% select(-PC1) %>% filter(cell_line != "Parental") %>% filter(cell_line != "Naive")) #task5e_pca_scores %>% select(-PC1)
drug_allclones_dist <- calculate_pca_distances(avg_full_pca_scores %>% filter(cell_line != "Parental"))

comparison_allclones <- compare_distances(trans_allclones_dist, drug_allclones_dist)
plot_results_allclones <- plot_distance_comparison(comparison_allclones)
print(plot_results_allclones$plot)

##

## Task 7 - Plot data from the in lab validation experiments ----

#create a data frame of the HTS data with only the drugs that have been used for validation
hts_validation_data = panscreen_long %>% filter(Name %in% c("(+)-JQ1", "Dasatinib", "Clofarabine", "Saracatinib (AZD0530)", "OTX015", "Epirubicin HCl", "S3I-201", "Camptothecin")) %>%
  select(-Target) %>% select(-Pathway)

#Mutate the Viability column in the lab data so it's in percent format rather than ratio format to match the HTS data
validation_data = lab_validation_data %>% mutate(Toxicity = 100 - (100 * Viability)) %>% select(-Viability) %>%
  filter(CellLine %in% c("Parental", "D13", "E2", "C12", "D4"))

#set up a dataframe with all of the naming inconsistencies between the HTS data and the lab data
name_mapping <- data.frame(NonStandard = c("JQ1", "Saracatinib"), Standard = c("(+)-JQ1", "Saracatinib (AZD0530)"))

#loop through the mapping data frame and update names in the lab data frame
for (i in 1:nrow(name_mapping)) {
  validation_data$Name[validation_data$Name == name_mapping$NonStandard[i]] <- name_mapping$Standard[i]
}

#Calculate the average for each cell line/dose/drug/treatment combination for each data frame
avg_hts_validation_data = hts_validation_data %>% group_by(CellLine, Name, Dose, Treatment) %>% summarize(AVG_Toxicity_HTS = mean(Toxicity, na.rm = TRUE), .groups = 'drop')

avg_lab_validation_data = validation_data %>% group_by(CellLine, Name, Dose, Treatment) %>% summarize(AVG_Toxicity_LAB = mean(Toxicity, na.rm = TRUE), .groups = 'drop')

#Use full_join to combine the two data frames while retaining all rows
all_validation_data = full_join(avg_lab_validation_data, avg_hts_validation_data, by = c("Name", "Dose", "Treatment", "CellLine"))

ggplot(all_validation_data, aes(x = (100 - AVG_Toxicity_HTS), y = (100 - AVG_Toxicity_LAB), color = CellLine)) +
  geom_point(aes(), size = 2) +
  facet_wrap(~Name, nrow = 2, ncol = 4) +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  geom_abline(slope = 1, intercept = 0, color = "red") +  # Adds the y=x line 
  coord_cartesian(xlim = c(0, 125), ylim = c(0, 125)) +
  labs(title = "Comparison of Viability Between HTS and Validation Screens",
       x = "High-Throughput Screen % Viability",
       y = "In-Lab Screen % Viability",
       color = "Cell Line") +
  custom_figure_theme() +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
  theme(legend.background = element_blank(),
        legend.position = (c(0.6,0.2)))

ggsave(file.path(figure_path, "labvalidationdata_vs_htsresults.pdf"), width = 14.2, height = 8, units = "in")

##


## Task 14 - Compare parental responses to drugs between the screens and responses to vem vs dmso in lines 1 and 2 for the threshold analysis----

task14_panscreen_parental = panscreen_long %>% filter(CellLine == "Parental") %>% filter(Treatment == "DMSO")
task14_panscreen_vemdmso = panscreen_long %>% filter(CellLine == ("D13") | CellLine == ("E2")) 

ggplot((task14_panscreen_vemdmso %>% pivot_wider(names_from = CellLine,
                                                 values_from = Toxicity) %>% filter(Treatment == "Vemurafenib")), aes(x = (100 - D13), y = (100 - E2), color = as.factor(Dose))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +  # Adds the y=x line
  scale_color_discrete(name = "Dose", 
                       breaks = names(custom_doses),
                       labels = custom_doses) +
  labs(x = str_wrap("% Viability of a Drug-Dose Condition on Resistant Clone 2", width = 30),
       y = str_wrap("% Viability of a Drug-Dose Condition on Resistant Clone 1", width = 30),
       color = "Dose") +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-5, 175), ylim = c(-5, 175)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85, 0.15),  # Adjust these values to move the legend
    legend.background = element_blank()
  )

ggsave(file.path(figure_path, "scatternogray_d13ande2_vemurafenibonly.pdf"), width = 7.68, height = 7.04, units = "in")

task14_panscreen_parental = task14_panscreen_parental %>% 
  pivot_wider(names_from = Screen,
              values_from = Toxicity)
task14_panscreen_vemdmso = task14_panscreen_vemdmso %>% 
  pivot_wider(names_from = Treatment,
              values_from = Toxicity)

task14_panscreen_parental$Screen1_Toxicity = task14_panscreen_parental$"1"
task14_panscreen_parental$Screen2_Toxicity = task14_panscreen_parental$"2"

#Plot the toxicity for each dose/pair in the parental cell line in combination with DMSO only, Screen 1 vs Screen 2

ggplot(task14_panscreen_parental, aes(x = (100 - Screen1_Toxicity), y = (100 - Screen2_Toxicity), color = as.factor(Dose))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +  # Adds the y=x line
  scale_color_discrete(name = "Dose", 
                       breaks = names(custom_doses),
                       labels = custom_doses) +
  labs(x = str_wrap("Screen A % Viability of a Drug- Dose Condition on Drug-Naive", width = 31),
       y = str_wrap("Screen B % Viability of a Drug- Dose Condition on Drug-Naive", width = 31),
       color = "Dose") +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-5, 175), ylim = c(-5, 175)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.89, 0.3),  # Adjust these values to move the legend
    legend.background = element_blank()
  )

ggsave(file.path(figure_path, "scatternogray_screen1and2_parentalonly.pdf"), width = 7.68, height = 7.04, units = "in")

ggplot(task14_panscreen_parental, aes(x = (100 - Screen1_Toxicity), y = (100 - Screen2_Toxicity), color = as.factor(Dose))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +
  geom_abline(slope = 1, intercept = 20, color = "blue", linetype = "dashed" ) +
  geom_abline(slope = 1, intercept = -20, color = "blue", linetype = "dashed" ) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 30, linetype = "dashed", color = "blue") +
  scale_color_discrete(name = "Dose", 
                       breaks = names(custom_doses),
                       labels = custom_doses) +
  labs(x = str_wrap("Screen A % Viability of a Drug-Dose Condition on Drug-Naive", width = 40),
       y = str_wrap("Screen B % Viability of a Drug-Dose Condition on Drug-Naive", width = 40),
       color = "Dose") +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-5, 175), ylim = c(-5, 175)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.55, 0.9),  # Adjust these values to move the legend
    legend.background = element_blank()
  )

ggsave(file.path(figure_path, "scatternogray_screen1and2_parentalonly_forthresholddiagram.pdf"), width = 9.6, height = 9.6, units = "in")

ggplot(task14_panscreen_vemdmso, aes(x = (100 - DMSO), y = (100 - Vemurafenib), color = as.factor(Dose))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +  # Adds the y=x line
  scale_color_discrete(name = "Dose", 
                       breaks = names(custom_doses),
                       labels = custom_doses) +
  labs(x = "% Viability of a Drug-Dose Condition \nwith DMSO in Screen A",
       y = "% Viability of a Drug-Dose \nCondition with Vem. in Screen A",
       color = "Dose") +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-5, 175), ylim = c(-5, 175)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85, 0.2),  # Adjust these values to move the legend
    legend.background = element_blank())

ggsave(file.path(figure_path, "scatternogray_dmsoandvem_screen1.pdf"), width = 7.68, height = 7.04, units = "in")

#create a new column in the main variable for this task that is the difference in toxicity values between the two screens for each drug/dose pair
task14_panscreen_parental = task14_panscreen_parental %>% mutate(Toxicity_Delta = Screen2_Toxicity - Screen1_Toxicity)
task14_panscreen_vemdmso = task14_panscreen_vemdmso %>% mutate(Toxicity_Delta = Vemurafenib - DMSO)

#calculate the mean and the standard deviation at each dose

mean_toxicity_delta_10_14 <- mean((task14_panscreen_parental %>% filter(Dose == 10))$Toxicity_Delta, na.rm = TRUE)
sd_toxicity_delta_10_14 <- sd((task14_panscreen_parental %>% filter(Dose == 10))$Toxicity_Delta, na.rm = TRUE)
mean_toxicity_delta_1_14 <- mean((task14_panscreen_parental %>% filter(Dose == 1))$Toxicity_Delta, na.rm = TRUE)
sd_toxicity_delta_1_14 <- sd((task14_panscreen_parental %>% filter(Dose == 1))$Toxicity_Delta, na.rm = TRUE)
mean_toxicity_delta_01_14 <- mean((task14_panscreen_parental %>% filter(Dose == .1))$Toxicity_Delta, na.rm = TRUE)
sd_toxicity_delta_01_14 <- sd((task14_panscreen_parental %>% filter(Dose == .1))$Toxicity_Delta, na.rm = TRUE)
mean_toxicity_delta_001_14 <- mean((task14_panscreen_parental %>% filter(Dose == .01))$Toxicity_Delta, na.rm = TRUE)
sd_toxicity_delta_001_14 <- sd((task14_panscreen_parental %>% filter(Dose == .01))$Toxicity_Delta, na.rm = TRUE)

mean_treatment_delta_10_14 <- mean((task14_panscreen_vemdmso %>% filter(Dose == 10))$Toxicity_Delta, na.rm = TRUE)
sd_treatment_delta_10_14 <- sd((task14_panscreen_vemdmso %>% filter(Dose == 10))$Toxicity_Delta, na.rm = TRUE)
mean_treatment_delta_1_14 <- mean((task14_panscreen_vemdmso %>% filter(Dose == 1))$Toxicity_Delta, na.rm = TRUE)
sd_treatment_delta_1_14 <- sd((task14_panscreen_vemdmso %>% filter(Dose == 1))$Toxicity_Delta, na.rm = TRUE)
mean_treatment_delta_01_14 <- mean((task14_panscreen_vemdmso %>% filter(Dose == .1))$Toxicity_Delta, na.rm = TRUE)
sd_treatment_delta_01_14 <- sd((task14_panscreen_vemdmso %>% filter(Dose == .1))$Toxicity_Delta, na.rm = TRUE)
mean_treatment_delta_001_14 <- mean((task14_panscreen_vemdmso %>% filter(Dose == .01))$Toxicity_Delta, na.rm = TRUE)
sd_treatment_delta_001_14 <- sd((task14_panscreen_vemdmso %>% filter(Dose == .01))$Toxicity_Delta, na.rm = TRUE)

# Calculate the lower and upper bounds (2 standard deviations from the mean) at each dose

lower_bound_10_14 <- mean_toxicity_delta_10_14 - 2 * sd_toxicity_delta_10_14
upper_bound_10_14 <- mean_toxicity_delta_10_14 + 2 * sd_toxicity_delta_10_14
lower_bound_1_14 <- mean_toxicity_delta_1_14 - 2 * sd_toxicity_delta_1_14
upper_bound_1_14 <- mean_toxicity_delta_1_14 + 2 * sd_toxicity_delta_1_14
lower_bound_01_14 <- mean_toxicity_delta_01_14 - 2 * sd_toxicity_delta_01_14
upper_bound_01_14 <- mean_toxicity_delta_01_14 + 2 * sd_toxicity_delta_01_14
lower_bound_001_14 <- mean_toxicity_delta_001_14 - 2 * sd_toxicity_delta_001_14
upper_bound_001_14 <- mean_toxicity_delta_001_14 + 2 * sd_toxicity_delta_001_14

lower_treatment_bound_10_14 <- mean_treatment_delta_10_14 - 2 * sd_treatment_delta_10_14
upper_treatment_bound_10_14 <- mean_treatment_delta_10_14 + 2 * sd_treatment_delta_10_14
lower_treatment_bound_1_14 <- mean_treatment_delta_1_14 - 2 * sd_treatment_delta_1_14
upper_treatment_bound_1_14 <- mean_treatment_delta_1_14 + 2 * sd_treatment_delta_1_14
lower_treatment_bound_01_14 <- mean_treatment_delta_01_14 - 2 * sd_treatment_delta_01_14
upper_treatment_bound_01_14 <- mean_treatment_delta_01_14 + 2 * sd_treatment_delta_01_14
lower_treatment_bound_001_14 <- mean_treatment_delta_001_14 - 2 * sd_treatment_delta_001_14
upper_treatment_bound_001_14 <- mean_treatment_delta_001_14 + 2 * sd_treatment_delta_001_14

#plot a histogram of the toxicity deltas with vertical lines showing the cutoff for what we are considering not significantly different

ggplot(task14_panscreen_parental, aes(x = Toxicity_Delta, fill = as.factor(Dose))) +
  geom_histogram(binwidth = 1, position = "dodge") +
  geom_vline(aes(xintercept = lower_bound_10_14), color = "#C77CFF", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_bound_10_14), color = "#C77CFF", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_bound_1_14), color = "#00BFC4", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_bound_1_14), color = "#00BFC4", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_bound_01_14), color = "#7CAE00", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_bound_01_14), color = "#7CAE00", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_bound_001_14), color = "#F8766D", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_bound_001_14), color = "#F8766D", linetype = "dashed", linewidth = 1) +
  scale_fill_discrete(name = "Dose", 
                      breaks = names(custom_doses),
                      labels = custom_doses) +
  labs(title = "Comparison of Parental \nViability Across Screens",
       subtitle = paste("Dashed lines represent ±2 SD at each dose"),
       x = "Screen A to B % Viability Difference",
       y = "Number of Drug-Dose Conditions") +
  custom_figure_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.8, 0.3),  # Adjust these values to move the legend
    legend.background = element_blank()
  )

ggsave(file.path(figure_path, "histogramwithsdlines_screen1and2_parentaldiff.pdf"), width = 8.096, height = 7.04, units = "in")

ggplot(task14_panscreen_vemdmso, aes(x = Toxicity_Delta, fill = as.factor(Dose))) +
  geom_histogram(binwidth = 1, position = "dodge") +
  geom_vline(aes(xintercept = lower_treatment_bound_10_14), color = "#C77CFF", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_treatment_bound_10_14), color = "#C77CFF", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_treatment_bound_1_14), color = "#00BFC4", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_treatment_bound_1_14), color = "#00BFC4", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_treatment_bound_01_14), color = "#7CAE00", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_treatment_bound_01_14), color = "#7CAE00", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = lower_treatment_bound_001_14), color = "#F8766D", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = upper_treatment_bound_001_14), color = "#F8766D", linetype = "dashed", linewidth = 1) +
  scale_fill_discrete(name = "Dose", 
                      breaks = names(custom_doses),
                      labels = custom_doses) +
  labs(title = "Vemurafenib vs DMSO \nViability in Screen A",
       subtitle = paste("Dashed lines represent ±2 SD at each dose"),
       x = "Change in Viability between Treatments",
       y = "Number of Drug-Dose Conditions") +
  custom_figure_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.8, 0.3),  # Adjust these values to move the legend
    legend.background = element_blank()
  )

ggsave(file.path(figure_path, "histogramwithsdlines_vemdmso_screen1diff.pdf"), width = 8.096, height = 7.04, units = "in")

# Create the dose_thresholds data frame
dose_thresholds_14 <- data.frame(
  Dose = c(10, 1, 0.1, 0.01),  # Doses in numeric format
  lower_bound = c(lower_bound_10_14, lower_bound_1_14, lower_bound_01_14, lower_bound_001_14),
  upper_bound = c(upper_bound_10_14, upper_bound_1_14, upper_bound_01_14, upper_bound_001_14)
)

dose_thresholds_treatment_14 <- data.frame(
  Dose = c(10, 1, 0.1, 0.01),  # Doses in numeric format
  lower_bound = c(lower_treatment_bound_10_14, lower_treatment_bound_1_14, lower_treatment_bound_01_14, lower_treatment_bound_001_14),
  upper_bound = c(upper_treatment_bound_10_14, upper_treatment_bound_1_14, upper_treatment_bound_01_14, upper_treatment_bound_001_14)
)

# Calculate counts outside the bounds for each dose
counts_outside_bounds_14 <- task14_panscreen_parental %>%
  group_by(Dose) %>%
  left_join(dose_thresholds_14, by = "Dose") %>%
  filter(Toxicity_Delta < lower_bound | Toxicity_Delta > upper_bound) %>%
  summarize(Count = n())

counts_outside_treatment_bounds_14 <- task14_panscreen_vemdmso %>%
  group_by(Dose) %>%
  left_join(dose_thresholds_treatment_14, by = "Dose") %>%
  filter(Toxicity_Delta < lower_bound | Toxicity_Delta > upper_bound) %>%
  summarize(Count = n())

counts_outside_treatment_bounds_14 = counts_outside_treatment_bounds_14 %>% mutate(Count = Count/2)

ggplot(counts_outside_bounds_14, aes(x = as.factor(Dose), y = Count, fill = as.factor(Dose))) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Values Outside \nAcross-Screen Bounds",
       subtitle = "Bounds: ±2 SD (calculated per dose)",
       x = "Dose",
       y = "Number of Drug-Dose Conditions") +
  scale_fill_discrete(name = "Dose", 
                      breaks = names(custom_doses),
                      labels = custom_doses) +
  scale_x_discrete(name = "Dose",
                   breaks = names(custom_doses),
                   labels = custom_doses) +
  # Update reference lines to use actual counts
  geom_segment(data = counts_outside_bounds_14,
               aes(x = as.numeric(factor(Dose)) - 0.45,
                   xend = as.numeric(factor(Dose)) + 0.45,
                   y = Count,
                   yend = Count),
               color = "blue",
               linetype = "dashed") +
  custom_figure_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",  # Adjust these values to move the legend
        legend.background = element_blank())

ggsave(file.path(figure_path, "countsoutsidebounds_screen1and2_parentalonly.pdf"), width = 7.168, height = 7.36, units = "in")

ggplot(counts_outside_treatment_bounds_14, aes(x = as.factor(Dose), y = Count, fill = as.factor(Dose))) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Values Outside \nWithin-Screen Bounds",
       subtitle = "Bounds: ±2 SD (calculated per dose)",
       x = "Dose",
       y = "Number of Drug-Dose Conditions") +
  scale_fill_discrete(name = "Dose", 
                      breaks = names(custom_doses),
                      labels = custom_doses) +
  scale_x_discrete(name = "Dose",
                   breaks = names(custom_doses),
                   labels = custom_doses) +
  # Update reference lines to use actual counts
  geom_segment(data = counts_outside_treatment_bounds_14,
               aes(x = as.numeric(factor(Dose)) - 0.45,
                   xend = as.numeric(factor(Dose)) + 0.45,
                   y = Count,
                   yend = Count),
               color = "blue",
               linetype = "dashed") +
  custom_figure_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",  # Adjust these values to move the legend
        legend.background = element_blank())

ggsave(file.path(figure_path, "countsoutsidebounds_vemdmso_screen1diff.pdf"), width = 7.168, height = 7.36, units = "in")

## calculate coefficient of variance 

task14_panscreen_parental_over70tox = task14_panscreen_parental %>% filter(Screen1_Toxicity >= 70 | Screen2_Toxicity >= 70)

  task14_panscreen_parental_over70tox_10 = task14_panscreen_parental_over70tox %>% filter(Dose == 10) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_over70tox_10_cov = (sd(task14_panscreen_parental_over70tox_10$Toxicity_Delta))/(mean(task14_panscreen_parental_over70tox_10$Toxicity_Delta))
  
  task14_panscreen_parental_over70tox_1 = task14_panscreen_parental_over70tox %>% filter(Dose == 1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_over70tox_1_cov = (sd(task14_panscreen_parental_over70tox_1$Toxicity_Delta))/(mean(task14_panscreen_parental_over70tox_1$Toxicity_Delta))
  
  task14_panscreen_parental_over70tox_01 = task14_panscreen_parental_over70tox %>% filter(Dose == 0.1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_over70tox_01_cov = (sd(task14_panscreen_parental_over70tox_01$Toxicity_Delta))/(mean(task14_panscreen_parental_over70tox_01$Toxicity_Delta))
  
  task14_panscreen_parental_over70tox_001 = task14_panscreen_parental_over70tox %>% filter(Dose == 0.01) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_over70tox_001_cov = (sd(task14_panscreen_parental_over70tox_001$Toxicity_Delta))/(mean(task14_panscreen_parental_over70tox_001$Toxicity_Delta))
  
task14_panscreen_parental_midtox = task14_panscreen_parental %>% filter((Screen1_Toxicity >= 30 & Screen1_Toxicity <= 70) | (Screen2_Toxicity >= 30 & Screen2_Toxicity <= 70))

  task14_panscreen_parental_midtox_10 = task14_panscreen_parental_midtox %>% filter(Dose == 10) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_midtox_10_cov = (sd(task14_panscreen_parental_midtox_10$Toxicity_Delta))/(mean(task14_panscreen_parental_midtox_10$Toxicity_Delta))
  
  task14_panscreen_parental_midtox_1 = task14_panscreen_parental_midtox %>% filter(Dose == 1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_midtox_1_cov = (sd(task14_panscreen_parental_midtox_1$Toxicity_Delta))/(mean(task14_panscreen_parental_midtox_1$Toxicity_Delta))
  
  task14_panscreen_parental_midtox_01 = task14_panscreen_parental_midtox %>% filter(Dose == 0.1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_midtox_01_cov = (sd(task14_panscreen_parental_midtox_01$Toxicity_Delta))/(mean(task14_panscreen_parental_midtox_01$Toxicity_Delta))
  
  task14_panscreen_parental_midtox_001 = task14_panscreen_parental_midtox %>% filter(Dose == 0.01) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_midtox_001_cov = (sd(task14_panscreen_parental_midtox_001$Toxicity_Delta))/(mean(task14_panscreen_parental_midtox_001$Toxicity_Delta))
  
task14_panscreen_parental_under30tox = task14_panscreen_parental %>% filter(Screen1_Toxicity <= 30 | Screen2_Toxicity <= 30)

  task14_panscreen_parental_under30tox_10 = task14_panscreen_parental_under30tox %>% filter(Dose == 10) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_under30tox_10_cov = (sd(task14_panscreen_parental_under30tox_10$Toxicity_Delta))/(mean(task14_panscreen_parental_under30tox_10$Toxicity_Delta))
  
  task14_panscreen_parental_under30tox_1 = task14_panscreen_parental_under30tox %>% filter(Dose == 1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_under30tox_1_cov = (sd(task14_panscreen_parental_under30tox_1$Toxicity_Delta))/(mean(task14_panscreen_parental_under30tox_1$Toxicity_Delta))
  
  task14_panscreen_parental_under30tox_01 = task14_panscreen_parental_under30tox %>% filter(Dose == 0.1) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_under30tox_01_cov = (sd(task14_panscreen_parental_under30tox_01$Toxicity_Delta))/(mean(task14_panscreen_parental_under30tox_01$Toxicity_Delta))
  
  task14_panscreen_parental_under30tox_001 = task14_panscreen_parental_under30tox %>% filter(Dose == 0.01) %>% mutate(Toxicity_Delta = abs(Toxicity_Delta))
  task14_panscreen_parental_under30tox_001_cov = (sd(task14_panscreen_parental_under30tox_001$Toxicity_Delta))/(mean(task14_panscreen_parental_under30tox_001$Toxicity_Delta))
  
cov_table = data.frame(
  `Viability Range` = c("Viability < 30%", "Viability < 30%", "Viability < 30%", "Viability < 30%", "30% < Viability < 70%", "30% < Viability < 70%", "30% < Viability < 70%", "30% < Viability < 70%", "Viability > 70%", "Viability > 70%", "Viability > 70%", "Viability > 70%"),
  `Dose (µM)` = c(10, 1, 0.1, 0.01, 10, 1, 0.1, 0.01, 10, 1, 0.1, 0.01),
  `Coefficient of Variation (%)` = round(c(task14_panscreen_parental_over70tox_10_cov, task14_panscreen_parental_over70tox_1_cov, task14_panscreen_parental_over70tox_01_cov, task14_panscreen_parental_over70tox_001_cov, task14_panscreen_parental_midtox_10_cov, task14_panscreen_parental_midtox_1_cov, task14_panscreen_parental_midtox_01_cov, task14_panscreen_parental_midtox_001_cov, task14_panscreen_parental_under30tox_10_cov, task14_panscreen_parental_under30tox_1_cov, task14_panscreen_parental_under30tox_01_cov, task14_panscreen_parental_under30tox_001_cov) * 100, 2),
  check.names = FALSE
  )

pdf("cov_table.pdf", width = 8, height = 6)
grid.newpage()
g <- tableGrob(cov_table, rows = NULL, theme = ttheme_default(rowhead=list(fg_params=list(cex=0))))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(g), l = 1, r = ncol(g))
g <- gtable_add_grob(g,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(g))
grid.draw(g)
dev.off()



##

## Task 15 - PCA of screens 1 & 2 data, NOT averaged by replicates, with 10 uM drug/dose pairs removed ----

#remove the screen index column
task15_line_by_drugpair = dplyr::select(panscreen_long, -c(Target, Pathway))

#keep only the data of resistant colonies treated in combination with vemurafenib and the parental cells treated with DMSO only
task15_line_by_drugpair = task15_line_by_drugpair %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                                               (CellLine != "Parental" & Treatment != "DMSO"))

#remove the drug/dose pairs at 10 uM
task15_line_by_drugpair = task15_line_by_drugpair %>% filter(Dose != 10)

#create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
task15_line_by_drugpair_table = task15_line_by_drugpair %>%
  pivot_wider(
    names_from = c(Name, Dose),   # Create columns from CellLine
    values_from = Toxicity  # Fill the cells with toxicity values
  )

#remove the treatment column
task15_line_by_drugpair_table = dplyr::select(task15_line_by_drugpair_table, -Treatment)

#combine the cell line and screen columns
task15_line_by_drugpair_table <- task15_line_by_drugpair_table %>%
  unite("CellLine_Screen", CellLine, Screen, sep = "_", remove = TRUE)

#transpose the table 
task15_line_by_drugpair_table = t(task15_line_by_drugpair_table)

#set the cell line row as the column names
colnames(task15_line_by_drugpair_table) = as.character(unlist(task15_line_by_drugpair_table[1, ]))

#remove the cell line row
task15_line_by_drugpair_table = task15_line_by_drugpair_table[-1, ]

#flip the table back
task15_line_by_drugpair_table = t(task15_line_by_drugpair_table)

#convert the table to numeric by column
task15_line_by_drugpair_table_numeric = apply(task15_line_by_drugpair_table, 1, as.numeric)

#flip the table back since the previous line transposed it
task15_line_by_drugpair_table_numeric = t(task15_line_by_drugpair_table_numeric)

#set the column names as the column names of the previous version of the table
colnames(task15_line_by_drugpair_table_numeric) = colnames(task15_line_by_drugpair_table)

#run pca
task15_pca_result <- prcomp(task15_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)

#print the results 
print(summary(task15_pca_result))

#store the values for each pc for each cell line as a data frame
task15_pca_data <- as.data.frame(task15_pca_result$x)

#store the loadings for each pc for each feature as a data frame
task15_pca_loadings <- task15_pca_result$rotation

#store pcs 1:4 in a new data frame
task15_pca_scores <- as.data.frame(task15_pca_result$x[, 1:4])

#add a column to the new dataframe with the cell line info
task15_pca_scores$cell_line <- rownames(task15_line_by_drugpair_table_numeric)

task15_pca_scores$custom_labels <- c("Resistant Clone 2 (A)", "Resistant Clone 1 (A)", "Drug-Naive (A)", "Resistant Clone 4 (B)", "Resistant Clone 3 (B)", "Drug-Naive (B)")  # Add your desired labels

#plot the data for two pcs WITH cell line labels on each point
ggplot(task15_pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) + 
  geom_text_repel(aes(label = custom_labels), max.overlaps = Inf, size = 6) +
  custom_figure_theme() +
  coord_cartesian(xlim = c(-85, 75), ylim = c(-85, 75)) + 
  labs(color = "Row Number",
       x = "Principal Component 1 - 32.68%",
       y = "Principal Component 2 - 23.46%",
       title = "PCA of Screen A and B Viability Data")

ggsave(file.path(figure_path, "pc1and2_htsdata_no10uM_notaveraged.pdf"), width = 9.6, height = 9.6, units = "in")


##

## Task 16 - PCA of combined screens 1/2 data and follow up panel a data, NOT averaged by replicates ----

#Turn the toxicity_replicate columns into a toxicity column and a replicate column
controls_followup_panela_task16 <- controls_followup_panela_long %>%
  pivot_longer(
    cols = starts_with("Toxicity_Rep"),
    names_to = c("Replicate"),
    names_pattern = "Toxicity_Rep(\\d)",
    values_to = "Toxicity"
  ) %>%
  mutate(Replicate = as.integer(Replicate)) # Ensure Replicate is an integer

#add a screen ID column to the controls data frame to match the main follow up panel a data
controls_followup_panela_task16 = controls_followup_panela_task16 %>% mutate(Screen = rep("A", length(controls_followup_panela_task16$Name)))

#remove the drugs that didn't replicate in this validation screen 
followup_panela_long_task16 = dosenorm_followup_panela_long %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                        "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570"))

#add the control panel a data to the main follow up panel a dataset
full_followup_panela_task16 <- rbind(followup_panela_long_task16, controls_followup_panela_task16)

#remove the target and pathway columns from the panscreen data
panscreen_long_task16 = panscreen_long %>% select(-c(Target, Pathway))

#Add a replicate column to the panscreen data
panscreen_long_task16 = panscreen_long_task16 %>% mutate(Replicate = panscreen_long_task16$Screen)

#filter the panscreen data to get only drug/dose/treatment pairs that exist in the follow up panel a data
panscreen_long_task16 = panscreen_long_task16 %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                                           (CellLine != "Parental" & Treatment != "DMSO"))

panscreen_long_task16 = panscreen_long_task16 %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                     "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570",
                                                                     "Dasatinib", "Epirubicin HCl", "Clofarabine", "S3I-201"))

#filter the follow up panel a data to get only drug/dose/treatment pairs that exist in the panscreen data
full_followup_panela_task16 = full_followup_panela_task16 %>% filter(Dose %in% c(10, 1, .1, .01))

#combine the cell line, screen, and replicate columns for the panscreen data
panscreen_long_task16 <- panscreen_long_task16 %>%
  unite("CellLine_Screen_Replicate", CellLine, Screen, Replicate, sep = "_", remove = TRUE)

#combine the cell line, screen, and replicate columns for the follow up panel a data
full_followup_panela_task16 <- full_followup_panela_task16 %>%
  unite("CellLine_Screen_Replicate", CellLine, Screen, Replicate, sep = "_", remove = TRUE)

#create a unique identifier for each drug-dose pair in both data frames
panscreen_long_task16 <- panscreen_long_task16 %>%
  unite("Name_Dose", Name, Dose, sep = "_", remove = FALSE)

full_followup_panela_task16 <- full_followup_panela_task16 %>%
  unite("Name_Dose", Name, Dose, sep = "_", remove = FALSE)

#Find the common drug-dose pairs
common_pairs_task16 <- intersect(panscreen_long_task16$Name_Dose, full_followup_panela_task16$Name_Dose)

#Filter both data frames to keep only the common pairs
panscreen_long_task16_filtered <- panscreen_long_task16 %>%
  filter(Name_Dose %in% common_pairs_task16)

full_followup_panela_task16_filtered <- full_followup_panela_task16 %>%
  filter(Name_Dose %in% common_pairs_task16)

#Combine the filtered data frames
combined_df_task16 <- bind_rows(panscreen_long_task16_filtered, full_followup_panela_task16_filtered)

#Remove the temporary Name_Dose column if no longer needed
combined_df_task16 <- combined_df_task16 %>%
  select(-Name_Dose)

#create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
task16_line_by_drugpair_table = combined_df_task16 %>%
  pivot_wider(
    names_from = c(Name, Dose),   # Create columns from CellLine
    values_from = Toxicity  # Fill the cells with toxicity values
  )

#remove the treatment column
task16_line_by_drugpair_table = dplyr::select(task16_line_by_drugpair_table, -Treatment)

#transpose the table 
task16_line_by_drugpair_table = t(task16_line_by_drugpair_table)

#set the cell line row as the column names
colnames(task16_line_by_drugpair_table) = as.character(unlist(task16_line_by_drugpair_table[1, ]))

#remove the cell line row
task16_line_by_drugpair_table = task16_line_by_drugpair_table[-1, ]

#flip the table back
task16_line_by_drugpair_table = t(task16_line_by_drugpair_table)

#convert the table to numeric by column
task16_line_by_drugpair_table_numeric = apply(task16_line_by_drugpair_table, 1, as.numeric)

#flip the table back since the previous line transposed it
task16_line_by_drugpair_table_numeric = t(task16_line_by_drugpair_table_numeric)

#set the column names as the column names of the previous version of the table
colnames(task16_line_by_drugpair_table_numeric) = colnames(task16_line_by_drugpair_table)

#run pca
task16_pca_result <- prcomp(task16_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)

#print the results 
print(summary(task16_pca_result))

#store the values for each pc for each cell line as a data frame
task16_pca_data <- as.data.frame(task16_pca_result$x)

#store the loadings for each pc for each feature as a data frame
task16_pca_loadings <- task16_pca_result$rotation

write.csv(task16_pca_loadings, "task16_pca_loadings.csv")

#plot the data for two pcs WITHOUT cell line labels on each point
ggplot(task16_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  custom_theme()

#store pcs 1:4 in a new data frame
task16_pca_scores <- as.data.frame(task16_pca_result$x[, 1:4])

#add a column to the new dataframe with the cell line info
task16_pca_scores$cell_line <- rownames(task16_line_by_drugpair_table_numeric)

#plot the data for two pcs WITH cell line labels on each point
ggplot(task16_pca_scores, aes(x = PC3, y = PC4)) +
  geom_point(size = 3) + 
  geom_text_repel(aes(label = cell_line), max.overlaps = Inf) +
  custom_theme() +
  labs(color = "Row Number",
       x = "Principal Component 3",
       y = "Principal Component 4",
       title = "PCA Plot")

##

## Task 17 - PCA of combined screens 1/2 data and follow up panel a data, NOT averaged by replicates, with 10 uM doses removed ----

#Turn the toxicity_replicate columns into a toxicity column and a replicate column
controls_followup_panela_task17 <- controls_followup_panela_long %>%
  pivot_longer(
    cols = starts_with("Toxicity_Rep"),
    names_to = c("Replicate"),
    names_pattern = "Toxicity_Rep(\\d)",
    values_to = "Toxicity"
  ) %>%
  mutate(Replicate = as.integer(Replicate)) # Ensure Replicate is an integer

#add a screen ID column to the controls data frame to match the main follow up panel a data
controls_followup_panela_task17 = controls_followup_panela_task17 %>% mutate(Screen = rep("A", length(controls_followup_panela_task17$Name)))

#remove the drugs that didn't replicate in this validation screen 
followup_panela_long_task17 = dosenorm_followup_panela_long %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                                   "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570"))

#add the control panel a data to the main follow up panel a dataset
full_followup_panela_task17 <- rbind(followup_panela_long_task17, controls_followup_panela_task17)

#remove the target and pathway columns from the panscreen data
panscreen_long_task17 = panscreen_long %>% select(-c(Target, Pathway))

#Add a replicate column to the panscreen data
panscreen_long_task17 = panscreen_long_task17 %>% mutate(Replicate = panscreen_long_task17$Screen)

#filter the panscreen data to get only drug/dose/treatment pairs that exist in the follow up panel a data
panscreen_long_task17 = panscreen_long_task17 %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                                           (CellLine != "Parental" & Treatment != "DMSO"))

panscreen_long_task17 = panscreen_long_task17 %>% filter(Name %in% c("(+)-JQ1", "Autophinib", "Fluvastatin Sodium", "OTX015","SB273005","SR18662","TPX-0005","YM155 (Sepantronium Bromide)",
                                                                     "Carfilzomib (PR-171)", "Deguelin", "EC330", "GMX1778 (CHS828)", "IWP-O1", "JH-RE-06", "NSC228155", "VLX1570",
                                                                     "Dasatinib", "Epirubicin HCl", "Clofarabine", "S3I-201"))

#filter the follow up panel a data to get only drug/dose/treatment pairs that exist in the panscreen data
full_followup_panela_task17 = full_followup_panela_task17 %>% filter(Dose %in% c(10, 1, .1, .01))

#combine the cell line, screen, and replicate columns for the panscreen data
panscreen_long_task17 <- panscreen_long_task17 %>%
  unite("CellLine_Screen_Replicate", CellLine, Screen, Replicate, sep = "_", remove = TRUE)

#combine the cell line, screen, and replicate columns for the follow up panel a data
full_followup_panela_task17 <- full_followup_panela_task17 %>%
  unite("CellLine_Screen_Replicate", CellLine, Screen, Replicate, sep = "_", remove = TRUE)

#create a unique identifier for each drug-dose pair in both data frames
panscreen_long_task17 <- panscreen_long_task17 %>%
  unite("Name_Dose", Name, Dose, sep = "_", remove = FALSE)

full_followup_panela_task17 <- full_followup_panela_task17 %>%
  unite("Name_Dose", Name, Dose, sep = "_", remove = FALSE)

#Find the common drug-dose pairs
common_pairs_task17 <- intersect(panscreen_long_task17$Name_Dose, full_followup_panela_task17$Name_Dose)

#Filter both data frames to keep only the common pairs
panscreen_long_task17_filtered <- panscreen_long_task17 %>%
  filter(Name_Dose %in% common_pairs_task17)

full_followup_panela_task17_filtered <- full_followup_panela_task17 %>%
  filter(Name_Dose %in% common_pairs_task17)

#Combine the filtered data frames
combined_df_task17 <- bind_rows(panscreen_long_task17_filtered, full_followup_panela_task17_filtered)

#Remove the temporary Name_Dose column if no longer needed
combined_df_task17 <- combined_df_task17 %>%
  select(-Name_Dose)

#remove the drug/dose pairs at 10 uM
combined_df_task17 = combined_df_task17 %>% filter(Dose != 10)

#create a table where rows are cell lines and columns are drug/dose columns (NOTE: excludes treatment)
task17_line_by_drugpair_table = combined_df_task17 %>%
  pivot_wider(
    names_from = c(Name, Dose),   # Create columns from CellLine
    values_from = Toxicity  # Fill the cells with toxicity values
  )

#remove the treatment column
task17_line_by_drugpair_table = dplyr::select(task17_line_by_drugpair_table, -Treatment)

#transpose the table 
task17_line_by_drugpair_table = t(task17_line_by_drugpair_table)

#set the cell line row as the column names
colnames(task17_line_by_drugpair_table) = as.character(unlist(task17_line_by_drugpair_table[1, ]))

#remove the cell line row
task17_line_by_drugpair_table = task17_line_by_drugpair_table[-1, ]

#flip the table back
task17_line_by_drugpair_table = t(task17_line_by_drugpair_table)

#convert the table to numeric by column
task17_line_by_drugpair_table_numeric = apply(task17_line_by_drugpair_table, 1, as.numeric)

#flip the table back since the previous line transposed it
task17_line_by_drugpair_table_numeric = t(task17_line_by_drugpair_table_numeric)

#set the column names as the column names of the previous version of the table
colnames(task17_line_by_drugpair_table_numeric) = colnames(task17_line_by_drugpair_table)

#run pca
task17_pca_result <- prcomp(task17_line_by_drugpair_table_numeric, center = TRUE, scale. = TRUE)

#print the results 
print(summary(task17_pca_result))

#store the values for each pc for each cell line as a data frame
task17_pca_data <- as.data.frame(task17_pca_result$x)

#store the loadings for each pc for each feature as a data frame
task17_pca_loadings <- task17_pca_result$rotation

write.csv(task17_pca_loadings, "task17_pca_loadings.csv")

#plot the data for two pcs WITHOUT cell line labels on each point
ggplot(task17_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  custom_theme()

#store pcs 1:4 in a new data frame
task17_pca_scores <- as.data.frame(task17_pca_result$x[, 1:4])

#add a column to the new dataframe with the cell line info
task17_pca_scores$cell_line <- rownames(task17_line_by_drugpair_table_numeric)

#plot the data for two pcs WITH cell line labels on each point
ggplot(task17_pca_scores, aes(x = PC3, y = PC4)) +
  geom_point(size = 3) + 
  geom_text_repel(aes(label = cell_line), max.overlaps = Inf) +
  custom_theme() +
  labs(color = "Row Number",
       x = "Principal Component 3",
       y = "Principal Component 4",
       title = "PCA Plot")

##

## Task 19 - Create pairwise scatterplots of each of the four resistant lines used in each HTS  ----

task19_data = panscreen_long %>% filter(CellLine != "Parental") %>% filter(Treatment == "Vemurafenib")

task19_data$Screen = NULL

task19_data = task19_data %>% spread(CellLine, Toxicity)

task19_data = task19_data %>% group_by(Name, Dose, Treatment, Target, Pathway) %>%
  summarize(
    D13 = mean(D13, na.rm = TRUE),
    E2 = mean(E2, na.rm = TRUE),
    D4 = mean(D4, na.rm = TRUE),
    C12 = mean(C12, na.rm = TRUE)
  ) %>%
  filter(!is.na(D13) | !is.na(E2) | !is.na(D4) | !is.na(C12))

#plot the scatter of screen 1 toxicity vs screen 2 toxicity with the shaded region representing the values that are inside three standard deviations, with those dots grayed out, faceted by dose
x_range_19 <- range(task19_data$D13, na.rm = TRUE)
extended_x_range_19 <- c(x_range_19[1] - diff(x_range_19), x_range_19[2] + diff(x_range_19))

#Create a combined plot for all six unique cell line pairs

# Function to calculate toxicity delta and count values outside bounds
calculate_counts_19 <- function(data, col1, col2, threshold_df) {
  data %>%
    mutate(Toxicity_Delta = !!sym(col1) - !!sym(col2)) %>%
    # Join with the threshold dataframe to get dose-specific bounds
    left_join(threshold_df, by = "Dose") %>%
    # Filter using the dose-specific bounds
    filter(Toxicity_Delta < lower_bound | Toxicity_Delta > upper_bound) %>%
    group_by(Dose) %>%
    summarize(Count = n()) %>%
    mutate(Comparison = paste(col1, "-", col2))
}

# Calculate counts for all combinations
cellline_combinations_19 <- list(
  c("D13", "E2"), c("D13", "D4"), c("D13", "C12"),
  c("E2", "D4"), c("E2", "C12"), c("D4", "C12")
)

all_counts_19 <- map_df(cellline_combinations_19, ~calculate_counts_19(task19_data, .x[1], .x[2], dose_thresholds_14))

ggplot(all_counts_19, aes(x = as.factor(Dose), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_segment(data = counts_outside_bounds_14,
               aes(x = as.numeric(factor(Dose)) - 0.45,
                   xend = as.numeric(factor(Dose)) + 0.45,
                   y = Count,
                   yend = Count),
               inherit.aes = FALSE,
               color = "blue",
               linetype = "dashed") +
  labs(title = str_wrap("Number of Drug-Dose Conditions Outside of Bounds when Comparing Between Cell Lines", width = 50),
       subtitle = str_wrap("Dashed Blue Lines: Number Expected from Across-Screen Threshold", width = 40),
       x = "Dose",
       y = "Number of Drug-Dose \nConditions Outside 2 SDs") +
  scale_fill_brewer(palette = "Set2", 
                    labels = custom_cellline_pairs) +
  scale_x_discrete(name = "Dose",
                   breaks = names(custom_doses),
                   labels = custom_doses) +
  scale_y_continuous(breaks = seq(0, max(all_counts_19$Count), by = 25)) +
  custom_figure_theme() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  coord_cartesian(ylim = c(0, max(all_counts_19$Count) * 1.1))

ggsave(file.path(figure_path,"countsoutsidebounds_allcelllinepairs.pdf"), width = 13.4, height = 7.68, units = "in")

all_counts_treatment_19 <- map_df(cellline_combinations_19, ~calculate_counts_19(task19_data, .x[1], .x[2], dose_thresholds_treatment_14))

ggplot(all_counts_treatment_19, aes(x = as.factor(Dose), y = Count, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_segment(data = counts_outside_treatment_bounds_14,
               aes(x = as.numeric(factor(Dose)) - 0.45,
                   xend = as.numeric(factor(Dose)) + 0.45,
                   y = Count,
                   yend = Count),
               inherit.aes = FALSE,
               color = "blue",
               linetype = "dashed") +
  labs(title = str_wrap("Number of Drug-Dose Conditions Outside of Bounds when Comparing Between Cell Lines", width = 50),
       subtitle = str_wrap("Dashed Blue Lines: Number Expected from Within-Screen Threshold", width = 40),
       x = "Dose",
       y = "Number of Drug-Dose \nConditions Outside 2 SDs") +
  scale_fill_brewer(palette = "Set2", 
                    labels = custom_cellline_pairs) +
  scale_x_discrete(name = "Dose",
                   breaks = names(custom_doses),
                   labels = custom_doses) +
  scale_y_continuous(breaks = seq(0, max(all_counts_treatment_19$Count), by = 25)) +
  custom_figure_theme() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank()) +
  coord_cartesian(ylim = c(0, max(all_counts_treatment_19$Count) * 1.1))

ggsave(file.path(figure_path, "countsoutsidetreatmentbounds_allcelllinepairs.pdf"), width = 13.4, height = 7.68, units = "in")


##

## Task 20 - Create grouped bar plots for dose/drug pairs that are unique to each cell line, shared, etc across ALL lines including parental ----

task20_panscreen = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                               (CellLine != "Parental" & Treatment == "Vemurafenib"))

task20_pathwaytotals = task20_panscreen %>% distinct(Name, Pathway) %>%  # Remove duplicate drug-pathway combinations
  count(Pathway, name = "total_drugs")

#transform the data frame into a wide format and summarize (only applicable for parental)
task20_panscreenwide = task20_panscreen %>% select(-Screen) %>% select(-Treatment) %>% pivot_wider(names_from = c(CellLine), values_from = Toxicity, values_fn = mean)

task20_panscreenwide <- task20_panscreenwide %>%
  mutate(Type = case_when(
    # Universal: All values are > 70
    D13 > 70 & E2 > 70 & Parental > 70 & C12 > 70 & D4 > 70 ~ "Universal",
    
    # Unique: Exactly one value > 70 and at least one value < 30
    (rowSums(across(c(D13, E2, Parental, C12, D4), ~ . > 70)) == 1) & 
      (rowSums(across(c(D13, E2, Parental, C12, D4), ~ . < 30)) >= 1) ~ "Unique",
    
    # Shared: At least one value > 70 and at least one value < 30
    (rowSums(across(c(D13, E2, Parental, C12, D4), ~ . > 70)) >= 1) & 
      (rowSums(across(c(D13, E2, Parental, C12, D4), ~ . < 30)) >= 1) ~ "Shared",
    
    # Default: Empty string if none of the conditions are met
    TRUE ~ "Other"
  ))

task20_panscreenwide = task20_panscreenwide %>% filter(Dose != 10) %>% filter(Type != "Other")

task20_panscreenwide_long <- task20_panscreenwide %>%
  pivot_longer(
    cols = c(D13, E2, Parental, C12, D4),
    names_to = "CellLine",
    values_to = "Toxicity"
  )

task20_panscreenwide_long = task20_panscreenwide_long %>% filter(Toxicity >= 70) %>% select(-Toxicity)

write_xlsx(task20_panscreenwide_long, file.path(figure_path, "list_uniqueshareduniversal_drugs.xlsx"))

complete_data_pathways <- task20_panscreenwide_long %>%
  select(Pathway, Type) %>%
  distinct() %>%
  crossing(CellLine = unique(task20_panscreenwide_long$CellLine)) %>%
  left_join(
    task20_panscreenwide_long %>% 
      group_by(Pathway, CellLine, Type) %>% 
      summarise(n = n_distinct(Name), .groups = 'drop'),
    by = c("Pathway", "CellLine", "Type")
  ) %>%
  mutate(n = replace_na(n, 0))

complete_data_pathways_normalized <- complete_data_pathways %>% left_join(task20_pathwaytotals, by = "Pathway") %>% mutate(Normalized_n = n/total_drugs)

ggplot(complete_data_pathways %>% 
         filter(Type == "Universal") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = unique(task20_panscreenwide_long$Pathway))),
       aes(x = Pathway, y = n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Universal Drugs",
       x = "Pathway",
       y = str_wrap("Number of Drugs", width = 10),
       fill = "Cell Line") +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = c(.9,.84),
        legend.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = -20, r = 0, b = 20, l = 0)))
                          
ggsave(file.path(figure_path, "allcommoncombos_bypathway_barplot_notfaceted.pdf"), width = 21.44, height = 5.48, units = "in", device = cairo_pdf)

ggplot(complete_data_pathways_normalized %>% 
         filter(Type == "Universal") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = unique(task20_panscreenwide_long$Pathway))),
       aes(x = Pathway, y = Normalized_n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Universal Drugs",
       x = "Pathway",
       y = str_wrap("Number of Drugs", width = 10),
       fill = "Cell Line") +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = c(.8,.84),
        legend.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = -20, r = 0, b = 20, l = 0)))

ggsave(file.path(figure_path, "allcommoncombos_bynormalizedpathway_barplot_notfaceted.pdf"), width = 21.44, height = 5.48, units = "in", device = cairo_pdf)

ggplot(complete_data_pathways %>%
         filter(Type == "Universal") %>%
         group_by(CellLine) %>%
         summarize(total = sum(n)) %>%
         mutate(percent = total/sum(total) * 100,
                label = paste0(round(percent, 1), "%"),
                CellLine = factor(CellLine, levels = names(custom_names)))) +
         geom_col(aes(x = "", y = percent, fill = CellLine)) +
         geom_text(aes(x = "", y = cumsum(percent) - percent/2, 
                       label = total),
                   size = 6, color = "white") +
         coord_polar(theta = "y", start = 0) +
         labs(title = str_wrap("Total per Cell Line", width = 9),
              x = NULL,
              y = NULL) +
         scale_fill_manual(values = custom_colors,
                           labels = custom_names,
                           breaks = names(custom_names)) +
         custom_figure_theme() +
         theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_blank(),    
               legend.position = "none")

ggsave(file.path(figure_path, "allcommoncombos_bycellline_piechart.pdf"), width = 3.2, height = 3.2, units = "in")

ggplot(complete_data_pathways %>% 
         filter(Type == "Shared") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = (unique(task20_panscreenwide_long$Pathway)))),
       aes(x = Pathway, y = n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Shared Drugs",
       x = "",
       y = str_wrap("Number of Drugs", width = 10)) +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path, "allsharedcombos_bypathway_barplot_notfaceted.pdf"), width = 21.44, height = 2.6, units = "in")

ggplot(complete_data_pathways_normalized %>% 
         filter(Type == "Shared") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = (unique(task20_panscreenwide_long$Pathway)))),
       aes(x = Pathway, y = Normalized_n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Shared Drugs",
       x = "",
       y = str_wrap("Number of Drugs", width = 10)) +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path, "allsharedcombos_bynormalizedpathway_barplot_notfaceted.pdf"), width = 21.44, height = 2.6, units = "in")

complete_data_pathways_shared = complete_data_pathways %>%
  filter(Type == "Shared") %>%
  group_by(CellLine) %>%
  summarize(total = sum(n)) %>%
  # Arrange by the order in custom_names
  arrange(match(CellLine, names(sort(custom_names)))) %>%
  mutate(CellLine = factor(CellLine, levels = unique(CellLine))) %>%
  # Calculate positions after arranging
  mutate(
    cs = rev(cumsum(rev(total))),
    pos = total/2 + lead(cs, 1),
    pos = if_else(is.na(pos), total/2, pos)
  )

ggplot(complete_data_pathways %>%
         filter(Type == "Shared") %>%
         group_by(CellLine) %>%
         summarize(total = sum(n)) %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names))), aes(x = "", y = total, fill = fct_reorder(CellLine, custom_names[CellLine]))) +
  geom_col(width = 1) +
  coord_polar(theta = "y", start = 0) +
  geom_text_repel(aes(y = pos, label = total),
                   data = complete_data_pathways_shared,
                   size = 6,
                   show.legend = FALSE,
                   nudge_x = 0.2,
                   color = "white",
                   segment.color = NA,
                   max.overlaps = Inf,
                   min.segment.length = 0,
                   force = 1,          # Adjust force of repulsion
                   direction = "y",    # Constrain movement to y direction
                   xlim = c(-0.5, 0.5)) +  # Remove label border
  labs(title = str_wrap("Total per Cell Line", width = 9),
       x = NULL,
       y = NULL) +
  scale_fill_manual(values = custom_colors,
                    labels = custom_names,
                    breaks = names(custom_names)) +
  custom_figure_theme() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),    
        legend.position = "none")

ggsave(file.path(figure_path,"allsharedcombos_bycellline_piechart.pdf"), width = 3.2, height = 3.2, units = "in")

ggplot(complete_data_pathways %>% 
         filter(Type == "Unique") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = (unique(task20_panscreenwide_long$Pathway)))),
       aes(x = Pathway, y = n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Breakdown of Unique Drugs",
       x = "",
       y = "                           Number of \nDrugs") +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        # Add these lines to control the appearance of the second x-axis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        plot.margin = margin(l = 0, r = 0, unit = "pt"),
        axis.title.y = element_text(angle = 90, hjust = 0.9)) + 
   scale_y_break(breaks = c(6, 30), ticklabels = c(30, 35))

ggsave(file.path(figure_path, "alluniquecombos_bypathway_barplot_notfaceted.pdf"), width = 21.76, height = 2.6, units = "in")

ggplot(complete_data_pathways_normalized %>% 
         filter(Type == "Unique") %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names)),
                Pathway = factor(Pathway, levels = (unique(task20_panscreenwide_long$Pathway)))),
       aes(x = Pathway, y = Normalized_n, fill = CellLine)) +
  geom_col(position = position_dodge2(width = 0.9),
           width = 0.8) +
  coord_cartesian(clip = "off") +
  labs(title = "Breakdown of Unique Drugs",
       x = "",
       y = "                           Number of \nDrugs") +
  scale_fill_manual(values = custom_colors, 
                    labels = custom_names, 
                    breaks = names(custom_names)) +
  scale_color_manual(values = custom_colors,  # Add color scale with same values
                     labels = custom_names, 
                     breaks = names(custom_names)) +
  scale_x_discrete(limits = (unique(task20_panscreenwide_long$Pathway))) +
  custom_figure_theme() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        # Add these lines to control the appearance of the second x-axis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        plot.margin = margin(l = 0, r = 0, unit = "pt"),
        axis.title.y = element_text(angle = 90, hjust = 0.9)) + 
  scale_y_break(breaks = c(0.05, 0.398), ticklabels = c(0.4))

ggsave(file.path(figure_path, "alluniquecombos_bynormalizedpathway_barplot_notfaceted.pdf"), width = 21.76, height = 2.6, units = "in")

complete_data_pathways_unique = complete_data_pathways %>%
  filter(Type == "Unique") %>%
  group_by(CellLine) %>%
  summarize(total = sum(n)) %>%
  # Arrange by the order in custom_names
  arrange(match(CellLine, names(sort(custom_names)))) %>%
  mutate(CellLine = factor(CellLine, levels = unique(CellLine))) %>%
  # Calculate positions after arranging
  mutate(
    cs = rev(cumsum(rev(total))),
    pos = total/2 + lead(cs, 1),
    pos = if_else(is.na(pos), total/2, pos)
  )

ggplot(complete_data_pathways %>%
         filter(Type == "Unique") %>%
         group_by(CellLine) %>%
         summarize(total = sum(n)) %>%
         mutate(CellLine = factor(CellLine, levels = names(custom_names))), aes(x = "", y = total, fill = fct_reorder(CellLine, custom_names[CellLine]))) +
  geom_col(width = 1) +
  coord_polar(theta = "y", start = 0) +
  geom_text_repel(aes(y = pos, label = total),
                  data = complete_data_pathways_unique,
                  size = 6,
                  show.legend = FALSE,
                  nudge_x = 0.2,
                  color = "white",
                  segment.color = NA,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  force = 0,          # Adjust force of repulsion
                  direction = "y",    # Constrain movement to y direction
                  xlim = c(-0.5, 0.5)) +  # Remove label border
  labs(title = str_wrap("Total per Cell Line", width = 9),
       x = NULL,
       y = NULL) +
  scale_fill_manual(values = custom_colors,
                    labels = custom_names,
                    breaks = names(custom_names)) +
  custom_figure_theme() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),    
        legend.position = "none")

ggsave(file.path(figure_path,"alluniquecombos_bycellline_piechart.pdf"), width = 3.2, height = 3.2, units = "in")

#Create stacked bar plots that are sorted by the number of cell lines for which a drug_dose pair kills
transformed_data = plot_data %>% filter(!grepl("_10$", NameDose))
transformed_data = transformed_data %>%
  group_by(NameDose) %>%
  summarise(Count = n_distinct(CellLine))

ggplot(transformed_data, aes(x = Count)) +
  geom_bar() +
  labs(title = "Distribution of Drug/Dose Pairs with >70% Toxicity",
       x = "Number of Resistant Lines Killed",
       y = "Number of Drug/Dose Pairs") +
  custom_theme()

# Step 6 - make a stacked bar plot for all the drug/dose pairs that have > 70% toxicity for each resistant line and <30% toxicity for at least one other line

common_combinations_6 <- Reduce(intersect, sets_maxmin)
all_combinations_6 <- unique(unlist(sets_maxmin))

#Identify unique combinations for each cell line
unique_combinations_6 <- lapply(names(sets_maxmin), function(cell_line) {
  other_combinations_6 <- unlist(sets_maxmin[names(sets_maxmin) != cell_line])
  setdiff(sets_maxmin[[cell_line]], other_combinations_6)
})

names(unique_combinations_6) <- names(sets_maxmin)

#Identify shared combinations
shared_combinations_6 <- lapply(names(sets_maxmin), function(cell_line) {
  setdiff(sets_maxmin[[cell_line]], c(unique_combinations_6[[cell_line]], common_combinations_6))
})

names(shared_combinations_6) <- names(sets_maxmin)

#Create data for plotting
plot_data_6 <- list()
for (cell_line in names(sets_maxmin)) {
  unique_data_6 <- if(length(unique_combinations_6[[cell_line]]) > 0) 
    data.frame(NameDose = unique_combinations_6[[cell_line]], CellLine = cell_line, Type = "Unique") 
  else NULL
  shared_data_6 <- if(length(shared_combinations_6[[cell_line]]) > 0)
    data.frame(NameDose = shared_combinations_6[[cell_line]], CellLine = cell_line, Type = "Shared")
  else NULL
  common_data_6 <- if(length(common_combinations_6) > 0) 
    data.frame(NameDose = common_combinations_6, CellLine = cell_line, Type = "Common") 
  else NULL
  
  # Combine the data
  combined_data_6 <- rbind(unique_data_6, shared_data_6, common_data_6)
  if (!is.null(combined_data_6)) {
    plot_data_6[[cell_line]] <- combined_data_6
  }
}
plot_data_6 <- do.call(rbind, plot_data_6)

plot_data_6$Type <- factor(plot_data_6$Type, levels = c("Unique", "Shared", "Common"))

#Create stacked bar plots
ggplot((plot_data_6 %>%
          filter(!grepl("_10$", NameDose)) %>%
          mutate(CellLine = factor(CellLine, levels = names(custom_names)))), aes(x = CellLine, fill = Type)) +
  geom_bar() +
  custom_theme() + 
  labs(title = "Distribution of Drug-Dose Pairs with >70% Toxicity & <30% Toxicity",
       x = "Cell Line",
       y = "Number of Drug/Dose Pairs") +
  scale_x_discrete(name = "Cell Line",
                   breaks = names(custom_names),
                   labels = custom_names) +
  scale_fill_brewer(palette = "Set1")

ggplot((plot_data_6 %>%
          filter(!grepl("_10$", NameDose)) %>%
          mutate(CellLine = factor(CellLine, levels = names(custom_names)))), aes(x = CellLine, fill = Type)) +
  geom_bar() +
  geom_text(stat = "count",
            aes(label = Type),
            position = position_stack(vjust = 0.5),
            size = 7,
            color = "white") +
  custom_figure_theme() + 
  labs(x = "Cell Line",
       y = "Number of Toxic Drug/Dose Pairs") +
  scale_x_discrete(name = "Cell Line",
                   breaks = names(custom_names),
                   labels = function(x) str_wrap(custom_names[x], width = 10)) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none") +
  annotate("text", x = 3.5, y = 23, 
           label = str_wrap("Shared pairs are toxic to 2-4 cell lines, but not all", width = 30),
           color = "black", size = 7)

ggsave("plots_for_draftfigures/drugdosepairs_uniqueshared_barplot.pdf", width = 8, height = 7.264, units = "in")

plot_data_6_pathways <- plot_data_6 %>%
  separate(NameDose, into = c("Name", "Dose"), sep = "_", remove = TRUE)

# Create a named vector for lookup
pathway_lookup <- setNames(panscreen_long$Pathway, panscreen_long$Name)

# Add Pathway column based on the Name
plot_data_6_pathways <- plot_data_6_pathways %>%
  mutate(Pathway = pathway_lookup[Name])

pathway_lookup <- setNames(panscreen_long$Target, panscreen_long$Name)

plot_data_6_pathways <- plot_data_6_pathways %>%
  mutate(Target = pathway_lookup[Name])

ggplot((plot_data_6_pathways %>% filter(Dose != 10) %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), aes(x = CellLine, fill = CellLine)) +
  geom_bar(position = "dodge") +
  facet_wrap(~Pathway) +
  labs(title = "Pathway Breakdown of Shared and Unique Drug-Dose Pairs",
       x = "Pathway",
       y = "Number of Drug/Dose Pairs") +
  custom_figure_theme() +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  theme(axis.text.x = element_blank())

ggsave("plots_for_draftfigures/alllines_pathways_sharedunique_barplot.pdf", width = 16, height = 6.4, units = "in")

ggplot((plot_data_6_pathways %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), aes(x = CellLine, fill = CellLine)) +
  geom_bar(position = "dodge") +
  facet_wrap(~Pathway) +
  labs(title = "Pathway Breakdown of Shared and Unique Drug-Dose Pairs",
       x = "Pathway",
       y = "Number of Drug/Dose Pairs") +
  custom_figure_theme() +
  scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  theme(axis.text.x = element_blank())

ggsave("plots_for_draftfigures/alllines_alldoses_pathways_sharedunique_barplot.pdf", width = 16, height = 6.4, units = "in")

#create a data table of all of the dose pairs with >70% & <30% toxicity at 1 um or below
d13_1uM = task20_results[["D13"]] %>% filter(Dose != 10)
e2_1uM = task20_results[["E2"]] %>% filter(Dose != 10)
c12_1uM = task20_results[["C12"]] %>% filter(Dose != 10)
d4_1uM = task20_results[["D4"]] %>% filter(Dose != 10)
parental_1uM = task20_results[["Parental"]] %>% filter(Dose != 10)

d13_1uM_flipped = task20_results_flipped[["D13"]] %>% filter(Dose != 10)
e2_1uM_flipped = task20_results_flipped[["E2"]] %>% filter(Dose != 10)
c12_1uM_flipped = task20_results_flipped[["C12"]] %>% filter(Dose != 10)
d4_1uM_flipped = task20_results_flipped[["D4"]] %>% filter(Dose != 10)
parental_1uM_flipped = task20_results_flipped[["Parental"]] #%>% filter(Dose != 10)

#create a data table of all of the dose pairs with >70% at 1 um or below
d13_1uM_70 = task20_results_70only[["D13"]] %>% filter(Dose != 10)
e2_1uM_70 = task20_results_70only[["E2"]] %>% filter(Dose != 10)
c12_1uM_70 = task20_results_70only[["C12"]] %>% filter(Dose != 10)
d4_1uM_70 = task20_results_70only[["D4"]] %>% filter(Dose != 10)
parental_1uM_70 = task20_results_70only[["Parental"]] %>% filter(Dose != 10)

## Calculate the difference between the parental response and the average resistant line response

task20_panscreenwide_resavg = task20_panscreenwide %>% mutate(ResistantAVG = ((D13 + C12 + D4 + E2)/4))

task20_panscreenwide_resavg = task20_panscreenwide_resavg %>% mutate(ResDelta = ResistantAVG - Parental)

task20_panscreenwide_resavg_filtered = task20_panscreenwide_resavg %>% filter(Dose != 10)

##

## Task 21 - Create rank order bar plots sorted by toxicity for each cell line ----

#filter the data to only capture parental cells treated without and resistant cells treated with vemurafenib
task21_panscreen = panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                               (CellLine != "Parental" & Treatment == "Vemurafenib")) %>% filter(Dose != 10)

task21_panscreen_avg = task21_panscreen %>% group_by(Name, Target, Pathway, CellLine, Treatment, Dose) %>%
                                            summarize(Toxicity = mean(Toxicity, na.rm = TRUE))

#transform the data frame into a wide format and summarize (only applicable for parental)
task21_panscreenwide = task21_panscreen %>% select(-Screen) %>% select(-Treatment) %>% pivot_wider(names_from = c(CellLine), values_from = Toxicity, values_fn = mean)

#filter out drug/dose pairs where not everything is killed
task21_panscreenwide_filtered = task21_panscreenwide %>% filter(!(D13 >= 50 & E2 >= 50 & C12 >= 50 & D4 >= 50 & Parental >= 50))

#filter out drug/dose pairs where everything is killed
task21_panscreenwide_filteredout = task21_panscreenwide %>% filter(D13 >= 50 & E2 >= 50 & C12 >= 50 & D4 >= 50 & Parental >= 50)

filteredout_pathways <- as.data.frame(table(task21_panscreenwide_filteredout$Pathway))
colnames(filteredout_pathways) <- c("Pathway", "Count")

# Sort by frequency in descending order
filteredout_pathways <- filteredout_pathways[order(-filteredout_pathways$Count),]

ggplot(filteredout_pathways, aes(x = reorder(Pathway, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  custom_figure_theme() +
  labs(title = "Frequency of Pathways in Drug-Dose Conditions \nwith < 50% Viability for all Cell Lines",
       x = "Pathway",
       y = "Number of Occurrences") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_path, "pathways_rankorder_filteredout.pdf"), width = 12.8, height = 9.6, units = "in", device = cairo_pdf)

# Calculate total unique drugs per pathway
pathway_totals <- task21_panscreenwide_filtered %>%
  distinct(Name, Pathway) %>%  # Remove duplicate drug-pathway combinations
  count(Pathway, name = "total_drugs")

# Now modify the get_top_pathways function to use this
get_top_pathways <- function(df, cell_line, n_pathways = 10) {
  df %>%
    filter(.data[[cell_line]] >= 70) %>%
    #distinct(Name, Pathway) %>%  # Count each drug only once per pathway
    count(Pathway) %>%
    # Join with pathway totals
    left_join(pathway_totals, by = "Pathway") %>%
    # Calculate normalized count (as percentage)
    mutate(normalized_count = (n / total_drugs) * 100) %>%
    # Arrange by normalized count instead of raw count
    arrange(desc(normalized_count)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= n_pathways) %>%
    mutate(
      Cell_Line = cell_line,
      Rank = rank
    ) %>%
    select(Cell_Line, Rank, Pathway, Count = n, Total_Drugs = total_drugs, Normalized_Percent = normalized_count)
}

# Apply function to each cell line
cell_lines <- c("D13", "E2", "Parental", "C12", "D4")
results <- do.call(rbind, lapply(cell_lines, function(cl) {
  get_top_pathways(task21_panscreenwide_filtered, cl, n_pathways = 8)
}))

# View results
print(results)

###

# Perform hypergeometric enrichment test
enrichment_analysis <- results %>%
  group_by(Cell_Line) %>%
  mutate(
    # Total drugs selected for this cell line
    Total_Selected = sum(Count),
    # Total drugs across all pathways (assuming no overlap)
    Total_Pool = sum(Total_Drugs),
    
    # Expected count if random
    Expected = (Total_Drugs / Total_Pool) * Total_Selected,
    
    # Hypergeometric test p-value
    # phyper(q, m, n, k, lower.tail = FALSE) where:
    # q = observed - 1 (for "more extreme")
    # m = pathway size (Total_Drugs for this pathway)
    # n = total drugs NOT in pathway
    # k = total drugs selected
    p_value = phyper(
      q = Count - 1,
      m = Total_Drugs,
      n = Total_Pool - Total_Drugs,
      k = Total_Selected,
      lower.tail = FALSE
    ),
    
    # Fold enrichment
    Fold_Enrichment = Count / Expected,
    
    # Enrichment direction
    Enrichment = case_when(
      Count > Expected ~ "Over-represented",
      Count < Expected ~ "Under-represented",
      TRUE ~ "As expected"
    )
  ) %>%
  ungroup() %>%
  group_by(Cell_Line) %>%
  mutate(
    # Adjust p-values for multiple testing (Benjamini-Hochberg)
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(Cell_Line, p_adj)

write_xlsx(enrichment_analysis, file.path(figure_path, "pathwaycounts_totals_normalizedpercents.xlsx"))

###

# Create a data frame for the legend
legend_df <- data.frame(
  Pathway = unique(results$Pathway),
  x = rep(1, length(unique(results$Pathway))),
  y = 1:length(unique(results$Pathway))
)

legend_df <- rbind(legend_df, (data.frame(Pathway = "Non-Top Pathways", x = 1, y = (length(unique(results$Pathway)) + 1))))

ggplot(legend_df, aes(x=0, y=0, fill=Pathway)) +
  geom_col() +
  scale_fill_manual(values = pathway_color_mapping,
                    labels = function(x) ifelse(x == "NF-\\u03BAB", "NF-κB", x)) +
  custom_figure_theme() +
  theme(
    legend.position = "left",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave(file.path(figure_path, "rankorderplot_pathwaylegend.pdf"), width = 4, height = 6.4, units = "in", device = cairo_pdf)

#plot the filtered data with the top 8 pathways colored for a single cell line
ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "Parental") %>% pull(Pathway)), Pathway, "Other"))), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), Parental), y = (100 - Parental), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "Parental") %>% arrange(desc(Count)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "Parental") %>%
                                 arrange(desc(Count)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  labs(title = "Drug-Naive Cells",
       x = "Drug-Dose Condition",
       y = "% Viability",
       fill = "Pathway") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_parental.pdf"), width = 11.84, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "Parental") %>% pull(Pathway)), Pathway, "Other")) %>% filter(Parental >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), Parental), y = (100 - Parental), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "Parental") %>% arrange(desc(Count)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "Parental") %>%
                                 arrange(desc(Count)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Drug-Naive Cells",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")
  
ggsave(file.path(figure_path,"rankorderplot_zoomin_allabove50removed_parental.pdf"), width = 11.84, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "Parental") %>% pull(Pathway)), Pathway, "Other")) %>% filter(Parental >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), Parental), y = (100 - Parental), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "Parental") %>% arrange(desc(Count)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "Parental") %>%
                                 arrange(desc(Count)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Drug-Naive Cells",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_text(size = 19.2, angle = 90, hjust = 1, vjust = 0.5),  
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_zoomin_allabove50removed_parental_withdrugnames.pdf"), width = 24, height = 9.6, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "D13") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D13 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D13), y = (100 - D13), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "D13") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "D13") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 2",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path, "rankorderplot_allabove50removed_d13.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "D13") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D13 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D13), y = (100 - D13), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "D13") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "D13") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 2",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_text(size = 19.2, angle = 90, hjust = 1, vjust = 0.5),  
        legend.position = "none")

ggsave(file.path(figure_path, "rankorderplot_allabove50removed_d13_withdrugnames.pdf"), width = 24, height = 7.36, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "E2") %>% pull(Pathway)), Pathway, "Other")) %>% filter(E2 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), E2), y = (100 - E2), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "E2") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "E2") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 1",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_e2.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "E2") %>% pull(Pathway)), Pathway, "Other")) %>% filter(E2 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), E2), y = (100 - E2), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "E2") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "E2") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 1",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_text(size = 19.2, angle = 90, hjust = 1, vjust = 0.5),  
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_e2_withdrugnames.pdf"), width = 24, height = 7.36, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "C12") %>% pull(Pathway)), Pathway, "Other")) %>% filter(C12 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), C12), y = (100 - C12), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "C12") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "C12") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 4",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_c12.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "C12") %>% pull(Pathway)), Pathway, "Other")) %>% filter(C12 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), C12), y = (100 - C12), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "C12") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "C12") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 4",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_text(size = 19.2, angle = 90, hjust = 1, vjust = 0.5),  
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_c12_withdrugnames.pdf"), width = 24, height = 7.36, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "D4") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D4 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D4), y = (100 - D4), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "D4") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "D4") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 3",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none") 

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_d4.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results %>% filter(Cell_Line == "D4") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D4 >= 70)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D4), y = (100 - D4), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results %>% filter(Cell_Line == "D4") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results %>% 
                                 filter(Cell_Line == "D4") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 3",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_text(size = 19.2, angle = 90, hjust = 1, vjust = 0.5),  
        legend.position = "none") 

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_d4_withdrugnames.pdf"), width = 24, height = 7.36, units = "in")

task21_panscreenwide_filtered_forheatmap = task21_panscreenwide_filtered %>% filter(D13 >= 70 | E2 >= 70 | C12 >= 70 | D4 >= 70 | Parental >= 70)

task21_panscreenwide_filtered_forheatmap = task21_panscreenwide_filtered_forheatmap %>% select(-Target)

task21_panscreenwide_filtered_forheatmap = task21_panscreenwide_filtered_forheatmap %>% mutate(DrugDoseCondition = paste(Name, Dose, sep = "_"))

task21_panscreenwide_filtered_forheatmap = task21_panscreenwide_filtered_forheatmap %>% select(-Dose) %>% select(-Name) 

task21_panscreenwide_filtered_forheatmap = task21_panscreenwide_filtered_forheatmap %>% mutate(D13 = 100 - D13) %>% 
                                                                                        mutate(E2 = 100 - E2) %>%
                                                                                        mutate(C12 = 100 - C12) %>% 
                                                                                        mutate(D4 = 100 - D4) %>%
                                                                                        mutate(Parental = 100 - Parental)

task21_heatmap_matrix <- as.matrix(task21_panscreenwide_filtered_forheatmap[, -7])  # Exclude the first column (Name) for the matrix
rownames(task21_heatmap_matrix) <- task21_panscreenwide_filtered_forheatmap$DrugDoseCondition  # Set row names as drug names
colnames(task21_heatmap_matrix) <- c("Pathway", "Resistant Clone 2", "Resistant Clone 1", "Drug-Naive", "Resistant Clone 4", "Resistant Clone 3")

heatmap_pathways <- task21_heatmap_matrix[, "Pathway"]

task21_heatmap_matrix_data <- task21_heatmap_matrix[,-1]

task21_heatmap_matrix_data <- apply(task21_heatmap_matrix_data, 2, as.numeric)

cell_order = c("Drug-Naive", "Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")
# Reorder the rows of the matrix
task21_heatmap_matrix_data <- task21_heatmap_matrix_data[ ,match(cell_order, colnames(task21_heatmap_matrix_data))]

heatmap21 <- Heatmap(task21_heatmap_matrix_data, 
                    name = "% Viability \nin HTS", 
                    col = (colorRamp2(seq(0, 125, length.out = 125), viridis(125))), 
                    show_row_names = FALSE,
                    column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    cluster_rows = TRUE, cluster_columns = FALSE,
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")),
                    left_annotation = rowAnnotation(
                      Pathway = heatmap_pathways,
                      col = list(Pathway = pathway_color_mapping),
                      annotation_name_gp = gpar(fontsize = 25.6, fontfamily = "Helvetica"),
                      annotation_legend_param = list(
                        Pathway = list(
                          title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                          labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))))

# Convert to a grob (graphical object)
heatmap21_grob <- grid.grabExpr(draw(heatmap21))

# Convert to a ggplot object
heatmap21_plot <- ggplotify::as.ggplot(heatmap21_grob)

ggsave(file.path(figure_path, "hts_heatmap_rankorderplot_drugs.pdf"), plot = heatmap21_plot, width = 24, height = 12, units = "in")

pathway_order <- order(task21_heatmap_matrix[, "Pathway"])
task21_heatmap_matrix <- task21_heatmap_matrix[pathway_order, ]

heatmap_pathways <- task21_heatmap_matrix[, "Pathway"]

task21_heatmap_matrix <- task21_heatmap_matrix[,-1]

task21_heatmap_matrix <- apply(task21_heatmap_matrix, 2, as.numeric)

cell_order = c("Drug-Naive", "Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")
# Reorder the rows of the matrix
task21_heatmap_matrix <- task21_heatmap_matrix[ ,match(cell_order, colnames(task21_heatmap_matrix))]


heatmap21b <- Heatmap(task21_heatmap_matrix_sorted, 
                     name = "% Viability \nin HTS", 
                     col = (colorRamp2(seq(0, 125, length.out = 125), viridis(125))), 
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                     cluster_rows = FALSE, cluster_columns = FALSE,
                     heatmap_legend_param = list(
                       title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                       labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")),
                     left_annotation = rowAnnotation(
                       Pathway = heatmap_pathways,
                       col = list(Pathway = pathway_color_mapping),
                       annotation_name_gp = gpar(fontsize = 25.6, fontfamily = "Helvetica"),
                       annotation_legend_param = list(
                         Pathway = list(
                           title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                           labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))))

# Convert to a grob (graphical object)
heatmap21b_grob <- grid.grabExpr(draw(heatmap21b))

# Convert to a ggplot object
heatmap21b_plot <- ggplotify::as.ggplot(heatmap21b_grob)

ggsave(file.path(figure_path, "hts_heatmap_rankorderplot_drugs_unclustered.pdf"), plot = heatmap21b_plot, width = 24, height = 12, units = "in")

# Now modify the get_top_pathways function to use this
get_top_pathways_new <- function(df, cell_line, n_pathways = 10) {
  df %>%
    filter(.data[[cell_line]] >= 90) %>%
    #distinct(Name, Pathway) %>%  # Count each drug only once per pathway
    count(Pathway) %>%
    # Join with pathway totals
    left_join(pathway_totals, by = "Pathway") %>%
    # Calculate normalized count (as percentage)
    mutate(normalized_count = (n / total_drugs) * 100) %>%
    # Arrange by normalized count instead of raw count
    arrange(desc(normalized_count)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= n_pathways) %>%
    mutate(
      Cell_Line = cell_line,
      Rank = rank
    ) %>%
    select(Cell_Line, Rank, Pathway, Count = n, Total_Drugs = total_drugs, Normalized_Percent = normalized_count)
}

# Apply function to each cell line
cell_lines <- c("D13", "E2", "Parental", "C12", "D4")
results_new <- do.call(rbind, lapply(cell_lines, function(cl) {
  get_top_pathways_new(task21_panscreenwide_filtered, cl, n_pathways = 8)
}))

# View results
print(results_new)

# Create a data frame for the legend
legend_df_new <- data.frame(
  Pathway = unique(results_new$Pathway),
  x = rep(1, length(unique(results_new$Pathway))),
  y = 1:length(unique(results_new$Pathway))
)

legend_df_new <- rbind(legend_df_new, (data.frame(Pathway = "Non-Top Pathways", x = 1, y = (length(unique(results_new$Pathway)) + 1))))

ggplot(legend_df_new, aes(x=0, y=0, fill=Pathway)) +
  geom_col() +
  scale_fill_manual(values = pathway_color_mapping,
                    labels = function(x) ifelse(x == "NF-\\u03BAB", "NF-κB", x)) +
  custom_figure_theme() +
  theme(
    legend.position = "left",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave(file.path(figure_path, "rankorderplot_pathwaylegend_new.pdf"), width = 4, height = 6.4, units = "in", device = cairo_pdf)

#plot the filtered data with the top 8 pathways colored for a single cell line
ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results_new %>% filter(Cell_Line == "Parental") %>% pull(Pathway)), Pathway, "Other")) %>% filter(Parental >= 90)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), Parental), y = (100 - Parental), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results_new %>% filter(Cell_Line == "Parental") %>% arrange(desc(Count)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results_new %>% 
                                 filter(Cell_Line == "Parental") %>%
                                 arrange(desc(Count)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Drug-Naive Cells",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_zoomin_allabove50removed_parental_new.pdf"), width = 11.84, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results_new %>% filter(Cell_Line == "D13") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D13 >= 90)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D13), y = (100 - D13), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results_new %>% filter(Cell_Line == "D13") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results_new %>% 
                                 filter(Cell_Line == "D13") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 2",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path, "rankorderplot_allabove50removed_d13_new.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results_new %>% filter(Cell_Line == "E2") %>% pull(Pathway)), Pathway, "Other")) %>% filter(E2 >= 90)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), E2), y = (100 - E2), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results_new %>% filter(Cell_Line == "E2") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results_new %>% 
                                 filter(Cell_Line == "E2") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 1",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_e2_new.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results_new %>% filter(Cell_Line == "C12") %>% pull(Pathway)), Pathway, "Other")) %>% filter(C12 >= 90)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), C12), y = (100 - C12), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results_new %>% filter(Cell_Line == "C12") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results_new %>% 
                                 filter(Cell_Line == "C12") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 4",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_c12_new.pdf"), width = 9.6, height = 3.2, units = "in")

ggplot((task21_panscreenwide_filtered %>% 
          mutate(Pathway_Color = ifelse(Pathway %in% (results_new %>% filter(Cell_Line == "D4") %>% pull(Pathway)), Pathway, "Other")) %>% filter(D4 >= 90)), 
       aes(x = reorder(paste(Name, Dose, sep = "_"), D4), y = (100 - D4), fill = Pathway_Color)) +
  geom_col() +
  scale_fill_manual(values = c(pathway_color_mapping[results_new %>% filter(Cell_Line == "D4") %>% arrange(desc(Normalized_Percent)) %>% pull(Pathway)], "Other" = "gray80"),
                    name = "Pathway (Ranked)",
                    breaks = c(results_new %>% 
                                 filter(Cell_Line == "D4") %>%
                                 arrange(desc(Normalized_Percent)) %>%
                                 pull(Pathway),
                               "Other")) +
  custom_figure_theme() +
  coord_cartesian(ylim = c(-5, 30)) +
  labs(title = "Resistant Clone 3",
       x = "Drug-Dose Condition",
       y = "% Viability") +
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none") 

ggsave(file.path(figure_path,"rankorderplot_allabove50removed_d4_new.pdf"), width = 9.6, height = 3.2, units = "in")

##

## Task 22 - Plot Control Drug Data from the follow up panel a vs from the original screens (NEEDS TO BE CLEANED UP) ----

#turn the two toxicity columns into a toxicity column and a replicate column
controls_followup_panela_long_rep <- controls_followup_panela_long %>%
  pivot_longer(cols = starts_with("Toxicity"),
               names_to = "Replicate",
               names_prefix = "Toxicity_",
               values_to = "Toxicity") %>%
  mutate(Replicate = as.numeric(gsub("Rep", "", Replicate)))

#filter the control follow up panel a data to only include cell lines that were used in the main screens
controls_followup_panela_long_filtered = controls_followup_panela_long_rep %>% filter(CellLine %in% c("Parental","D13","D4","C12","E2"))

#filter the panscreen data to only include cell line/treatment pairs also used in the validation follow up panel a
panscreen_long_filtered_5 <- panscreen_long %>%
  filter((CellLine == "Parental" & Treatment == "DMSO") | 
           (CellLine != "Parental" & Treatment == "Vemurafenib")
  )

panscreen_controls_long_filtered <- panscreen_controls_long %>%
  filter((CellLine == "Parental" & Treatment == "DMSO") | 
           (CellLine != "Parental" & Treatment == "Vemurafenib")
  )

panscreen_controls_long_filtered = panscreen_controls_long_filtered %>%
  pivot_longer(cols = starts_with("Toxicity"),
               names_to = "Replicate",
               names_prefix = "Toxicity_",
               values_to = "Toxicity") %>%
  mutate(Replicate = as.numeric(gsub("Rep", "", Replicate)))

#filter the pan screen data to only include drugs used in the control validation follow up panel a and remove the Target and Pathway columns
panscreen_long_filtered_5 <- panscreen_long_filtered_5 %>%
  filter(Name %in% controls_followup_panela_long_filtered$Name) %>% select(-Target) %>% select(-Pathway)

#group Validation follow up panel a Data by Name, CellLine, Treatment, and Dose then calculating the average Toxicity for each group
avg_controls_followup_panela_long_filtered <- controls_followup_panela_long_filtered %>%
  group_by(Name, CellLine, Treatment, Dose) %>%
  summarise(ScreenA_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

#group Pan Screen Data by Name, CellLine, Treatment, and Dose then calculating the average Toxicity for each group
avg_panscreen_long_filtered_5 <- panscreen_long_filtered_5 %>%
  group_by(Name, CellLine, Treatment, Dose) %>%
  summarise(PanScreen_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

avg_panscreen_controls_long_filtered = panscreen_controls_long_filtered %>%
  group_by(Drug, CellLine, Treatment, Dose) %>%
  summarise(PanScreen_Tox = mean(Toxicity, na.rm = TRUE)) %>%
  ungroup()

#Round the doses in the Validation follow up panel a data and filter only the ones used in screens 1 and 2
avg_controls_followup_panela_long_filtered <- avg_controls_followup_panela_long_filtered %>%
  # Convert Dose to numeric and round to 2 decimal places
  mutate(Dose = round(as.numeric(Dose), digits = 2)) %>%
  # Filter for specific dose values
  filter(Dose %in% c(10, 1, 0.1, 0.01))

#rename the drug column in the follow up panel a control data to match the pan screen data so it will be properly used to join the toxicity columns
avg_panscreen_controls_long_filtered = avg_panscreen_controls_long_filtered %>% rename(Name = Drug)

#Add the validation follow up panel a toxicity data to the pan screen data
controls_panscreen_followup_panela_overlap <- left_join(avg_panscreen_long_filtered_5, 
                                                        avg_controls_followup_panela_long_filtered %>% 
                                                          select(CellLine, Treatment, Dose, Name, ScreenA_Tox),
                                                        by = c("CellLine", "Treatment", "Dose", "Name"))

controls_only_overlap <- left_join(avg_panscreen_controls_long_filtered, 
                                   avg_controls_followup_panela_long_filtered %>% 
                                     select(CellLine, Treatment, Dose, Name, ScreenA_Tox),
                                   by = c("CellLine", "Treatment", "Dose", "Name"))

ggplot(controls_only_overlap, aes(x = PanScreen_Tox, y = ScreenA_Tox, color = CellLine)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed" ) +  # Adds the y=x line 
  facet_wrap(~Dose, labeller = labeller(Dose = as_labeller(custom_doses))) +
  labs(title = "Comparison of Toxicity Across Screens",
       x = "Toxicity in Screens 1 & 2",
       y = "Toxicity in Follow Up Panel A") +
  scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
  xlim(0, NA) +
  custom_theme() +
  coord_fixed()

##


## Task 23 - Make a heatmap of the panscreen data for all drugs targeting the PI3K/Akt/mTOR pathway or integrin pathway ----

panscreen_23 <- panscreen_long %>% filter((CellLine == "Parental" & Treatment == "DMSO") | 
                                              (CellLine != "Parental" & Treatment == "Vemurafenib"))

panscreen_23 = panscreen_23 %>% mutate(Viability = 100 - Toxicity) %>% select(-Toxicity)

#transform the data frame into a wide format and summarize (only applicable for parental)
panscreen_wide_23 = panscreen_23 %>% select(-Screen) %>% select(-Treatment) %>% pivot_wider(names_from = c(CellLine), values_from = Viability, values_fn = mean)

panscreen_wide_23_pi3k = panscreen_wide_23 %>% filter(grepl("PI3K/Akt/mTOR", Pathway) | grepl("PI3K", Target))

panscreen_wide_23_pi3k = panscreen_wide_23_pi3k %>% select(-Target) %>% select(-Pathway)

task23_heatmapdata <- panscreen_wide_23_pi3k %>% unite("Name_Dose", Name, Dose, sep = "_", remove = TRUE)

task23_heatmapdata <- task23_heatmapdata %>% filter(D13 < 90 | E2 < 90 | D4 < 90 | C12 < 90 | Parental < 90)

# Step 2: Convert the resulting dataframe to a matrix
panscreen_wide_23_matrix <- as.matrix(task23_heatmapdata[, -1])  # Exclude the first column (Name) for the matrix
rownames(panscreen_wide_23_matrix) <- task23_heatmapdata$Name_Dose  # Set row names as drug names
colnames(panscreen_wide_23_matrix) <- c("Resistant Clone 2", "Resistant Clone 1", "Drug-Naive", "Resistant Clone 4", "Resistant Clone 3")

cell_order = c("Drug-Naive","Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")

task23_reordered_matrix <- panscreen_wide_23_matrix[ ,match(cell_order, colnames(panscreen_wide_23_matrix))]

heatmap23 <- Heatmap(task23_reordered_matrix, 
                    name = "% Viability \n in HTS of \n PI3K Drugs", 
                    col = (colorRamp2(seq(0, 125, length.out = 125), viridis(125))), 
                    show_row_names = FALSE,
                    column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                    cluster_rows = TRUE, cluster_columns = FALSE,
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

print(heatmap23)

# Convert to a grob (graphical object)
heatmap23_grob <- grid.grabExpr(draw(heatmap23))

# Convert to a ggplot object
heatmap23_plot <- ggplotify::as.ggplot(heatmap23_grob)

ggsave(file.path(figure_path, "hts_heatmap_pi3kdrugs.pdf"), plot = heatmap23_plot, width = 13.76, height = 9.6, units = "in")

panscreen_wide_23_int = panscreen_wide_23 %>% filter(grepl("Integrin", Target))

panscreen_wide_23_int = panscreen_wide_23_int %>% select(-Target) %>% select(-Pathway)

task23b_heatmapdata <- panscreen_wide_23_int %>% unite("Name_Dose", Name, Dose, sep = "_", remove = TRUE)

task23b_heatmapdata <- task23b_heatmapdata %>% filter(D13 < 90 | E2 < 90 | D4 < 90 | C12 < 90 | Parental < 90)

# Step 2: Convert the resulting dataframe to a matrix
panscreen_wide_23b_matrix <- as.matrix(task23b_heatmapdata[, -1])  # Exclude the first column (Name) for the matrix
rownames(panscreen_wide_23b_matrix) <- task23b_heatmapdata$Name_Dose  # Set row names as drug names
colnames(panscreen_wide_23b_matrix) <- c("Resistant Clone 2", "Resistant Clone 1", "Drug-Naive", "Resistant Clone 4", "Resistant Clone 3")

cell_order = c("Drug-Naive" ,"Resistant Clone 1", "Resistant Clone 2", "Resistant Clone 3", "Resistant Clone 4")

task23b_reordered_matrix <- panscreen_wide_23b_matrix[ ,match(cell_order, colnames(panscreen_wide_23b_matrix))]

heatmap23b <- Heatmap(task23b_reordered_matrix, 
                     name = "% Viability \n in HTS of \nIntegrin Drugs", 
                     col = (colorRamp2(seq(0, 125, length.out = 125), viridis(125))), 
                     row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                     column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
                     cluster_rows = TRUE, cluster_columns = FALSE,
                     heatmap_legend_param = list(
                       title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                       labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

print(heatmap23b)

# Convert to a grob (graphical object)
heatmap23b_grob <- grid.grabExpr(draw(heatmap23b))

# Convert to a ggplot object
heatmap23b_plot <- ggplotify::as.ggplot(heatmap23b_grob)

ggsave(file.path(figure_path, "hts_heatmap_integrindrugs.pdf"), plot = heatmap23b_plot, width = 13.76, height = 9.6, units = "in")


## Task 24 - Plot data from the in lab validation experiments done for revisions ----

revision_round1_validation_data_adjusted = revision_round1_validation_data %>% mutate(Viability = Viability*100)

  selected_drug <- c("SR18662")
  
  selected_data = revision_round1_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)

  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 1 µM", selected_drug),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)
  
  ggsave(file.path(figure_path, "sr18662viability_1uM_revisions_round1.pdf"), width = 9.6, height = 9.6, units = "in")

  selected_drug <- c("Saracatinib (AZD 0530)")
  
  selected_data = revision_round1_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)

  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 10 µM", selected_drug),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)
  
  ggsave(file.path(figure_path, "saracatinibviability_10uM_revisions_round1.pdf"), width = 9.6, height = 9.6, units = "in")

revision_round2_validation_data_adjusted = revision_round2_validation_data %>% mutate(Viability = Viability*100)

  selected_drug <- c("Deguelin")
  
  selected_data = revision_round2_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)
  
  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("%.1f", after_stat(y))),
                 vjust = -0.5, size = 4, color = "black") +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 1 µM", selected_drug),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)

  ggsave(file.path(figure_path, "deguelinviability_1uM_revisions_round2.pdf"), width = 7.68, height = 7.68, units = "in")
  
  selected_drug <- c("SB273005")
  
  selected_data = revision_round2_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)
  
  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("%.1f", after_stat(y))),
                 vjust = -0.5, size = 4, color = "black") +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 0.1 µM", selected_drug),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)

  ggsave(file.path(figure_path, "sb273005viability_0_1uM_revisions_round2.pdf"), width = 7.68, height = 7.68, units = "in")
  
  selected_drug <- c("Deguelin + SB273005")
  
  selected_data = revision_round2_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)

  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    stat_summary(fun = mean, geom = "text", 
                 aes(label = sprintf("%.1f", after_stat(y))),
                 vjust = -0.5, size = 4, color = "black") +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 1 µM Deguelin + \n0.1 µM SB273005"),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)
  
  ggsave(file.path(figure_path, "deguelinsb273005viability_revisions_round2.pdf"), width = 7.68, height = 8.2, units = "in")
  
  selected_drug <- c("IWP-O1")
  
  selected_data = revision_round2_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)
  
  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 3.16 µM", selected_drug),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)
  
  ggsave(file.path(figure_path, "iwpo1viability_3_16uM_revisions_round2.pdf"), width = 7.68, height = 7.68, units = "in")
  
  selected_drug <- c("IWP-O1 + SB273005")
  
  selected_data = revision_round2_validation_data_adjusted %>% filter((CellLine == "Parental" & Treatment == "DMSO") | (CellLine != "Parental" & Treatment == "Vemurafenib"))
  
  selected_data <- selected_data %>% filter(Name %in% selected_drug)
  
  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)
  
  #plot the data for only the resistant colonies for a single drug in combination with vemurafenib only
  ggplot((selected_data %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = CellLine, y = Viability, fill = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "bar", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    scale_fill_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to 3.16 µM IWP-O1 + \n0.1 µM SB273005"),
         x = "",
         y = "% Viability") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") +
    scale_x_discrete(labels = custom_names)+
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    geom_point(data = selected_data,
               aes(x = CellLine, y = Viability), 
               size = 3, color = "black", alpha = 0.5)
  
  ggsave(file.path(figure_path, "iwpo1sb273005viability_revisions_round2.pdf"), width = 7.68, height = 8.2, units = "in")

selected_data = jq1_coculture_data %>% mutate(Viability = Viability*100)

  selected_data$Dose <- as.numeric(selected_data$Dose)
  selected_data$Treatment <- as.factor(selected_data$Treatment)
  
  ggplot((selected_data %>% filter(Treatment == "Vemurafenib") %>% mutate(CellLine = factor(CellLine, levels = names(custom_names)))), 
         aes(x = Dose, y = Viability, color = CellLine)) +
    # Plot lines for the average of the points
    stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1) +
    scale_color_manual(values = custom_colors, labels = custom_names, breaks = names(custom_names)) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(title = paste("Response to (+)-JQ1"),
         x = "Dose",
         y = "% Viability",
         color = "Cell Line") +
    custom_figure_theme() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.background = element_blank(),
      legend.position = c(0.8, 0.8)) +
    scale_x_log10(breaks = as.numeric(names(custom_doses)), labels = custom_doses) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125))
  
  ggsave(file.path(figure_path, "jq1cocultureviability.pdf"), width = 9, height = 9, units = "in")

