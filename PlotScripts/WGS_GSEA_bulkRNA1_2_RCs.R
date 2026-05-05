library(ggplot2)
library(tidyverse)
library(readxl)
library(here)
library(ggrepel)
library(fgsea)
library(gridGraphics)
library(grid)
library(scales)
library(RColorBrewer)
library(circlize)
library(viridis)
library(gridExtra)
library(magrittr)
library(Cairo)
library(VennDiagram)
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

## Task 0a - Load the data from the first bulk rna seq run and put it into a gene x cell line matrix ----

current_path <- here()

parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "BulkRNASeq_Run1")

figure_path <- file.path(parent_path, "Plots")

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

##

## Task 0b - Load the data from the second bulk rna seq run ----

current_path <- here()

parent_path <- dirname(current_path)

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

## Task 0c - Load the mutation data table ----

current_path <- here()

parent_path <- dirname(current_path)

new_path <- file.path(parent_path, "ExtractedData", "WholeGenomeSequencing")

mutations_table <- read_tsv(file.path(new_path, "cadd_phred_over_15.tsv"), 
                            skip = 6,
                            na = c("", "NA", "N/A"))

## Task 0d - Load custom theme and other custom functions ----

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

custom_colors <- c("Naive" = "#5E5E5E", "Parental" = "#5E5E5E", "D13" = "#FF7E79", "E2" = "#4294F8", "C12" = "#6539F8", "A15" = "#3E8D27", "E2 + D13" = "#8E1356", "C14" = "#953877", "D5" = "#6ADCC1", "D4" = "#9E6D15",
                   "A2" = "#EF8BF9", "A11" = "#CC5C76", "C3" = "#F57946", "E8" = "#F9AD2A", "F9" = "#1D457F", "I11" = "#625A94", "G12" = "#0013F0")

custom_names <- c("Naive" = "Drug-Naive Run 1", "Parental" = "Drug-Naive Run 2", "E2" = "Resistant Line #1", "D13" = "Resistant Line #2", "D4" = "Resistant Line #3", "C12" = "Resistant Line #4", "C14" = "Resistant Line #5",
                  "D5" = "Resistant Line #6", "F9" = "Resistant Line #7", "G12" = "Resistant Line #8", "I11" = "Resistant Line #9", "A2" = "Resistant Line #10", "A11" = "Resistant Line #11",
                  "C3" = "Resistant Line #12", "E8" = "Resistant Line #13")


keytarget_genes_of_interest <- c("KLF5", "SRC", "ITGAV", "ITGB3", "ITGB5", "BRD4", "BRD2", "BRD3", "BRDT", "PORCN", "BIRC5", "NAMPT", "NDUFS7", "NDUFS8", "NDUFV2", "NDUFS3", "NDUFS2", "NDUFV1", "NDUFS1",
                                 "NDUFS6", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "NDUFA12", "NDUFS4", "NDUFA9", "NDUFAB1", "NDUFA2", "NDUFA1", "NDUFB3", "NDUFA5", "NDUFA6", "NDUFA6", "NDUFA11",
                                 "NDUFB11", "NDUFS5", "NDUFB4", "NDUFA13", "NDUFB7", "NDUFA8", "NDUFB9", "NDUFB10", "NDUFB8", "NDUFC2", "NDUFB2", "NDUFA7", "NDUFA3", "NDUFB5", "NDUFB1", "NDUFC1", "NDUFA10",
                                 "NDUFV3", "NDUFB6", "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "STAT3", "JAK2", "YES", "FYN", "FGR", "LCK", "HCK", "BLK", "LYN", "FRK", "NRF2", "GSK3", "PRKAA1", "PRKAA2",
                                 "PRKAB1", "PRKAB2", "PRKAG1", "PRKAG2", "PRKAG3", "PIK3C3", "PIK3R4", "LIF", "NADPH", "REV1", "REV7", "EGFR", "HER2", "HER3", "HER4", "ALK", "ROS1", "TRKA", "TRKB", "TRKC")

##

## Task 6 - GSEA of Resistant Programs in Boe et al. ----

#Transform to long format
long_data <- run1_tpmdata %>%
  pivot_longer(
    cols = -c(sample, cell_line, replicate),  # Exclude these columns from pivoting
    names_to = "gene_name",                   # Column names become gene names
    values_to = "tpm"                         # Values become TPM values
  )

#Group by cell_line and gene_name, then calculate average TPM per gene across replicates
average_by_cell_line <- long_data %>%
  group_by(cell_line, gene_name) %>%
  summarize(avg_tpm = mean(tpm, na.rm = TRUE), .groups = "drop")

#Reshape back to wide format with genes as rows and cell lines as columns
gene_by_cell_line <- average_by_cell_line %>%
  pivot_wider(
    id_cols = gene_name,                      # Genes will be the row identifier
    names_from = cell_line,                   # Cell lines will become columns
    values_from = avg_tpm                     # Values will be the average TPM
  )

data = gene_by_cell_line %>%
  column_to_rownames("gene_name")

# Calculate log2 fold change for each cell line
ranked_lists_logfc <- lapply(colnames(data), function(cell_line) {
  other_cell_lines <- setdiff(colnames(data), cell_line)
  
  # Calculate the mean expression across all other cell lines
  mean_others <- rowMeans(data[, other_cell_lines])
  
  # Compute the log2 fold change
  logfc <- log2((data[, cell_line] + 1) / (mean_others + 1))
  
  # Rank the genes based on the log2 fold change
  logfc_ranked <- sort(logfc, decreasing = TRUE)
  
  return(logfc_ranked)
})

# Assuming 'data' is the original data frame with cell line names as column names
names(ranked_lists_logfc) <- colnames(data)

## Specify the gene set you want to use

metaprograms <- readRDS("metaprograms_list.rds") #ryan's metaprograms
# Assuming 'metaprograms' is a list of 13 gene sets

# Initialize a nested list to store fgsea results for each gene set and cell line
fgsea_results_nested_list <- list()

# Loop over each gene set in metaprograms
for (gene_set_name in names(metaprograms)) {
  # Create the gene set list for the current gene set
  current_gene_set <- list(metaprograms[[gene_set_name]])
  names(current_gene_set) <- gene_set_name
  
  # Initialize a list to store fgsea results for each cell line for the current gene set
  fgsea_results_list <- list()
  
  # Loop over each cell line
  for (i in seq_along(ranked_lists_logfc)) {
    # Running fgsea
    rankedGenesNoTies <- ranked_lists_logfc[[i]] + runif(length(ranked_lists_logfc[[i]]), min=-1e-10, max=1e-10)
    
    # Run fgsea without the nperm argument to use fgseaMultilevel
    fgsea_results_list[[i]] <- fgsea(pathways = current_gene_set, stats = rankedGenesNoTies, maxSize = 500)
    
    # Clean up
    rm(rankedGenesNoTies)
  } 
  
  # Name the results list with cell line names
  names(fgsea_results_list) <- colnames(data)
  
  # Store the results for the current gene set
  fgsea_results_nested_list[[gene_set_name]] <- fgsea_results_list
}

# Initialize an empty matrix to store NES values
nes_matrix <- matrix(nrow = length(fgsea_results_nested_list), ncol = length(fgsea_results_nested_list[[1]]), 
                     dimnames = list(names(fgsea_results_nested_list), names(fgsea_results_nested_list[[1]])))

# Populate the NES matrix
for (gene_set_name in names(fgsea_results_nested_list)) {
  for (cell_line_name in names(fgsea_results_nested_list[[gene_set_name]])) {
    # Extract NES
    nes_matrix[gene_set_name, cell_line_name] <- fgsea_results_nested_list[[gene_set_name]][[cell_line_name]]$NES[fgsea_results_nested_list[[gene_set_name]][[cell_line_name]]$pathway == gene_set_name]
  }
}

# Draw the heatmap
row_labels <- c("Ribosome", "Cell Division #1", "Cell Division #2", "Interferon", "Smooth Muscle", "EMT #1", "EMT #2", "Mitochondrial", "ECM", "Neural Crest", "Cytoskeleton", "Complement", "Melanocyte")

# Create a named vector for your labels, where names match the column names exactly
custom_names_sorted <- c(
  "A15" = "Additional Clone",
  "B13" = "Additional Clone", 
  "B4" = "Additional Clone",
  "B5" = "Additional Clone",
  "B9" = "Additional Clone",
  "C12" = "Resistant Clone 4",
  "C14" = "Additional Clone",
  "C16" = "Additional Clone",
  "C4" = "Additional Clone",
  "C5" = "Additional Clone",
  "C6" = "Additional Clone",
  "D1" = "Additional Clone",
  "D11" = "Additional Clone",
  "D13" = "Resistant Clone 2",
  "D4" = "Resistant Clone 3",
  "D5" = "Additional Clone",
  "D9" = "Additional Clone",
  "E2" = "Resistant Clone 1",
  "E9" = "Additional Clone",
  "Naive" = "Drug-Naive"
)

nes_heatmap <- Heatmap(nes_matrix,
        name = "Normalized \nEnrichment \nScore",
        col = viridis(128),
        row_labels = row_labels,
        column_labels = custom_names_sorted[colnames(nes_matrix)],
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
        column_names_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica"),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
          labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

# Convert to a grob (graphical object)
nes_heatmap_grob <- grid.grabExpr(draw(nes_heatmap))

# Convert to a ggplot object
nes_heatmap_plot <- ggplotify::as.ggplot(nes_heatmap_grob)

ggsave(file.path(figure_path, "nes_heatmap_resistanceprograms.pdf"), plot = nes_heatmap_plot, width = 22.4, height = 8, units = "in")

##




## Task 1 - Look at deleterious mutations list for presence of key target genes ----

mutations_table_simplified = mutations_table %>% select(Gene, Samples, Phred...14, cDNA_change)

mutations_table_simplified_filtered = mutations_table_simplified %>% filter(Gene %in% keytarget_genes_of_interest)

mutations_table_simplified_filtered = mutations_table_simplified_filtered %>% mutate(GeneChange = paste0(Gene, " (", cDNA_change, ")")) %>%
  select(-cDNA_change) %>% select(-Gene)

mutations_table_long <- mutations_table_simplified_filtered %>%
  mutate(row_id = row_number()) %>%
  separate_rows(Samples, Phred...14, sep = ";", convert = TRUE) %>%
  group_by(GeneChange, row_id) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  select(-row_id)

mutations_table_long = mutations_table_long %>% select(-position) %>% mutate(Samples = sub("\\..*", "", Samples))

mutations_table_long = mutations_table_long %>% rename(CellLine = Samples,
                                                       CADDScore = Phred...14)

mutations_table_wide = mutations_table_long %>% pivot_wider(names_from = CellLine,
                                                            values_from = CADDScore,
                                                            values_fill = 0)

mutations_table_wide = mutations_table_wide %>%
  column_to_rownames(var = "GeneChange")

mutation_colors = colorRamp2(
  c(0, seq(15, 25, length.out = 11)), 
  c("white", viridis(11))
)

mutations_matrix = as.matrix(mutations_table_wide)

heatmap_mutations <- Heatmap(mutations_matrix, 
                   name = "CADD\nScore", 
                   col = mutation_colors, 
                   cluster_rows = FALSE, 
                   cluster_columns = FALSE,
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
                     at = c(15, seq(16, 25, by = 1)),
                     title_gp = gpar(fontsize = 25.6, fontface = "bold", fontfamily = "Helvetica"),
                     labels_gp = gpar(fontsize = 19.2, fontfamily = "Helvetica")))

print(heatmap_mutations)

# Convert to a grob (graphical object)
heatmap_grob_mutations <- grid.grabExpr(draw(heatmap_mutations, padding = unit(c(3, 2, 2, 10), "mm"), heatmap_legend_side = "left"))

# Convert to a ggplot object
heatmap_plot_mutations <- ggplotify::as.ggplot(heatmap_grob_mutations)

ggsave(file.path(figure_path, "heatmap_mutationsofinterest_phredscore.pdf"), plot = heatmap_plot_mutations, width = 24, height = 8, units = "in")


## Task 2 - Look at expression of key target genes ----



#pull data from the first bulk RNA seq run
task2_run1_data = run1_tpmdata %>% select(-sample) %>% 
  filter(cell_line != "E9") %>%
  filter(!(cell_line == "A15" & replicate == "Early")) %>%
  filter(!(cell_line == "D11" & replicate == "Late2"))

#pull data from the second bulk RNA seq run
task2_run2_data = run2_tpmdata %>% select(-sample)

# Find common columns
task2_common_columns <- intersect(colnames(task2_run1_data), colnames(task2_run2_data))

# Subset both data frames
task2_run1_subset <- task2_run1_data[, task2_common_columns]
task2_run2_subset <- task2_run2_data[, task2_common_columns]

# Combine the subsetted data frames
task2_combined_data <- rbind(task2_run1_subset, task2_run2_subset)

#average by replicates
task2_combined_data_avg <- task2_combined_data %>%
  group_by(cell_line) %>%
  summarise(across(A1BG:last_col(), \(x) mean(x, na.rm = TRUE)))

task2_combined_data_avg = task2_combined_data_avg %>% filter(cell_line %in% c("Parental", "Naive", "E2", "D13", "D4", "C12", "C14", "D5", "F9", "G12", "I11", "A11", "A2", "C3", "E8"))

task2_combined_data_table <- t(task2_combined_data_avg)

#use the first row which contains the cell lines as the column names of the new matrix
colnames(task2_combined_data_table) = as.character(unlist(task2_combined_data_table[1, ]))

task2_combined_data_table = task2_combined_data_table[-1, ]

# Filter the matrix to only include genes in your list
filtered_task2_data_table <- task2_combined_data_table[rownames(task2_combined_data_table) %in% keytarget_genes_of_interest, ]

# Convert matrix to data frame and reshape to long format
plot_data <- filtered_task2_data_table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, 
               names_to = "Cell_Line", 
               values_to = "Expression") %>%
  mutate(Expression = as.numeric(Expression))

# Create the small multiples plot
ggplot(plot_data %>% mutate(Cell_Line = factor(Cell_Line, levels = names(custom_names))), aes(x = Cell_Line, y = Expression, color = Cell_Line)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = custom_colors, labels = custom_names) +
  scale_y_continuous(limits = c(0, NA)) +
  custom_figure_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 19.2, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = custom_names) +
  labs(
    x = "Cell Line",
    y = "Expression (TPM)",
    color = "Cell Line"
  )

ggsave(file.path(figure_path, "expressionofkeytargetgenes_allfollowuppanellines.pdf"), width = 24, height = 31, units = "in")


