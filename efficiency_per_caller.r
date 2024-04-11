# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

cell_lines_induced <- c("KIF5BE21d", "TFGC7_1d", "V1C1_1_d", "V3aA3_1d")
cell_lines_known <-c('C2924hc1', 'C824hc1', 'C924hc1', 'YO_24hc1')
cell_lines_public <-c("TCGA-78-7163","TCGA-86-A4P8","TCGA-50-8460","TCGA-67-6215","TCGA-05-4397","TCGA-97-8176","TCGA-97-7941","TCGA-05-4433","TCGA-05-4382","TCGA-55-7284","TCGA-50-8457","TCGA-50-6594","TCGA-50-5072","TCGA-55-A490","TCGA-55-8094","TCGA-97-8175","TCGA-L9-A5IP","TCGA-67-6217","TCGA-99-7458","TCGA-50-6593","TCGA-93-A4JN","TCGA-99-8028","TCGA-38-A44F","TCGA-55-A494","TCGA-53-A4EZ","TCGA-99-8025","TCGA-55-8089","TCGA-44-7661","TCGA-97-8174","TCGA-49-AAQV","TCGA-97-8179","TCGA-L9-A50W","TCGA-MN-A4N1","TCGA-MN-A4N5","TCGA-97-7937","TCGA-99-AA5R","TCGA-73-A9RS","TCGA-55-A48Z","TCGA-55-7913","TCGA-86-A4P8","TCGA-55-7914","TCGA-55-8513")


# Create a list of expected fusions for each sample
expected_fusions_induced <- c("KIF5B:ALK", "TFG:ALK", "EML4:ALK", "EML4:ALK")
expected_fusions_known <-rep("EML4:ALK", length(cell_lines_known)) # Bcs it is same for all known fusions
expected_fusions_public <- rep("EML4:ALK", length(cell_lines_public))  


STARfusion_induced <- read.table("~/data/starfusion_output/STARfusion_out_induced_ALK_data.txt", header = TRUE, sep = "\t")
arriba_induced <- read.table("~/data/arriba_output/arriba_induced_ALK_data.txt", header = TRUE, sep = "\t")
fcirc_induced <- read.table("~/data/fcirc_output/fcirc_induced_ALK_data.txt", header = TRUE, sep = "\t")
infusion_induced <- read.table("~/data/infusion_output/infusion_induced_ALK_data.txt", header = TRUE, sep = "\t")
fusioncatcher_induced <- read.table("~/data/fusioncatcher_output/fusioncatcher_induced_ALK_data.txt", header = TRUE, sep = "\t")
jaffa_d_induced <- read.table("~/data/jaffa_output/jaffa_d_induced_ALK_data.txt", header = TRUE, sep = "\t")
squid_induced<-read.table("~/data/squid_output/squid_induced_ALK_data.txt", header = TRUE, sep = "\t") 
#change all column names in gene1 and gene2
colnames(STARfusion_induced)[1:2] <- c("gene1", "gene2")
colnames(arriba_induced)[1:2] <- c("gene1", "gene2")
colnames(fcirc_induced)[1:2] <- c("gene1", "gene2")
colnames(infusion_induced)[1:2] <- c("gene1", "gene2")
colnames(fusioncatcher_induced)[1:2] <- c("gene1", "gene2")

squid_induced$FusedGenes <- as.character(squid_induced$FusedGenes)
squid_induced$gene1 <- sapply(strsplit(squid_induced$FusedGenes, ":"), `[`, 1)
squid_induced$gene2 <- sapply(strsplit(squid_induced$FusedGenes, ":"), `[`, 2)


#Genes in star has something like KIF5B^ENSG00000170759.11 (so i am getting rid of it)
STARfusion_induced$gene1 <- gsub("\\^.*", "", STARfusion_induced$gene1)
STARfusion_induced$gene2 <- gsub("\\^.*", "", STARfusion_induced$gene2)

# Reading the data files
STARfusion_known <- read.table("~/data/starfusion_output/STARfusion_out_known_ALK_data.txt", header = TRUE, sep = "\t")
arriba_known <- read.table("~/data/arriba_output/arriba_known_ALK_data.txt", header = TRUE, sep = "\t")
fcirc_known <- read.table("~/data/fcirc_output/fcirc_known_ALK_data.txt", header = TRUE, sep = "\t")
infusion_known <- read.table("~/data/infusion_output/infusion_known_ALK_data.txt", header = TRUE, sep = "\t")
fusioncatcher_known <- read.table("~/data/fusioncatcher_output/fusioncatcher_known_ALK_data.txt", header = TRUE, sep = "\t")
jaffa_d_known <- read.table("~/data/jaffa_output/jaffa_d_known_ALK_data.txt", header = TRUE, sep = "\t")
squid_known <- read.table("~/data/squid_output/squid_known_ALK_data.txt", header = TRUE, sep = "\t") 

# Renaming the columns
colnames(STARfusion_known)[1:2] <- c("gene1", "gene2")
colnames(arriba_known)[1:2] <- c("gene1", "gene2")
colnames(fcirc_known)[1:2] <- c("gene1", "gene2")
colnames(infusion_known)[1:2] <- c("gene1", "gene2")
colnames(fusioncatcher_known)[1:2] <- c("gene1", "gene2")

# Formatting the 'FusedGenes' column for the 'squid' dataset
squid_known$FusedGenes <- as.character(squid_known$FusedGenes)
squid_known$gene1 <- sapply(strsplit(squid_known$FusedGenes, ":"), `[`, 1)
squid_known$gene2 <- sapply(strsplit(squid_known$FusedGenes, ":"), `[`, 2)

# Removing additional identifiers from the 'gene1' and 'gene2' columns for the 'STARfusion' dataset
STARfusion_known$gene1 <- gsub("\\^.*", "", STARfusion_known$gene1)
STARfusion_known$gene2 <- gsub("\\^.*", "", STARfusion_known$gene2)

STARfusion_public <- read.table("~/data/starfusion_output/starfusion_public_ALK_data.txt", header = TRUE, sep = "\t")
arriba_public <- read.table("~/data/arriba_output/arriba_public_ALK_data.txt", header = TRUE, sep = "\t")
fcirc_public <- read.table("~/data/fcirc_output/fcirc_public_ALK_data.txt", header = TRUE, sep = "\t")
infusion_public <- read.table("~/data/infusion_output/infusion_public_ALK_data.txt", header = TRUE, sep = "\t")
fusioncatcher_public <- read.table("~/data/fusioncatcher_output/fusioncatcher_public_ALK_data.txt", header = TRUE, sep = "\t")
jaffa_d_public <- read.table("~/data/jaffa_output/jaffa_d_public_ALK_data.txt", header = TRUE, sep = "\t")


#change all column names in gene1 and gene2
colnames(STARfusion_public)[1:2] <- c("gene1", "gene2")
colnames(arriba_public)[1:2] <- c("gene1", "gene2")
colnames(fcirc_public)[1:2] <- c("gene1", "gene2")
colnames(infusion_public)[1:2] <- c("gene1", "gene2")
colnames(fusioncatcher_public)[1:2] <- c("gene1", "gene2")


#Genes in star has something like KIF5B^ENSG00000170759.11 (so i am getting rid of it)
STARfusion_public$gene1 <- gsub("\\^.*", "", STARfusion_public$gene1)
STARfusion_public$gene2 <- gsub("\\^.*", "", STARfusion_public$gene2) 




############# Changing the values in character ############

# Convert 'gene1' and 'gene2' columns to character in STARfusion_public
STARfusion_public$gene1 <- as.character(STARfusion_public$gene1)
STARfusion_public$gene2 <- as.character(STARfusion_public$gene2)

# Convert 'gene1' and 'gene2' columns to character in arriba_public
arriba_public$gene1 <- as.character(arriba_public$gene1)
arriba_public$gene2 <- as.character(arriba_public$gene2)

# Convert 'gene1' and 'gene2' columns to character in fcirc_public
fcirc_public$gene1 <- as.character(fcirc_public$gene1)
fcirc_public$gene2 <- as.character(fcirc_public$gene2)

# Convert 'gene1' and 'gene2' columns to character in infusion_public
infusion_public$gene1 <- as.character(infusion_public$gene1)
infusion_public$gene2 <- as.character(infusion_public$gene2)

# Convert 'gene1' and 'gene2' columns to character in fusioncatcher_public
fusioncatcher_public$gene1 <- as.character(fusioncatcher_public$gene1)
fusioncatcher_public$gene2 <- as.character(fusioncatcher_public$gene2)

# Convert 'gene1' and 'gene2' columns to character in jaffa_d_public
jaffa_d_public$gene1 <- as.character(jaffa_d_public$gene1)
jaffa_d_public$gene2 <- as.character(jaffa_d_public$gene2)
# Convert 'cell_line' column to character in STARfusion_public
STARfusion_public$cell_line <- as.character(STARfusion_public$cell_line)

# Convert 'cell_line' column to character in arriba_public
arriba_public$cell_line <- as.character(arriba_public$cell_line)

# Convert 'cell_line' column to character in fcirc_public
fcirc_public$cell_line <- as.character(fcirc_public$cell_line)

# Convert 'cell_line' column to character in infusion_public
infusion_public$cell_line <- as.character(infusion_public$cell_line)

# Convert 'cell_line' column to character in fusioncatcher_public
fusioncatcher_public$cell_line <- as.character(fusioncatcher_public$cell_line)

# Convert 'cell_line' column to character in jaffa_d_public
jaffa_d_public$cell_line <- as.character(jaffa_d_public$cell_line)

# Convert 'sub_dir' column to character in STARfusion_public
STARfusion_public$sub_dir <- as.character(STARfusion_public$sub_dir)

# Convert 'sub_dir' column to character in arriba_public
arriba_public$sub_dir <- as.character(arriba_public$sub_dir)

# Convert 'sub_dir' column to character in fcirc_public
fcirc_public$sub_dir <- as.character(fcirc_public$sub_dir)

# Convert 'sub_dir' column to character in infusion_public
infusion_public$sub_dir <- as.character(infusion_public$sub_dir)

# Convert 'sub_dir' column to character in fusioncatcher_public
fusioncatcher_public$sub_dir <- as.character(fusioncatcher_public$sub_dir)

# Convert 'sub_dir' column to character in jaffa_d_public
jaffa_d_public$sub_dir <- as.character(jaffa_d_public$sub_dir)

###########################################################################




# Loop through each dataset and add the relevant columns
datasets <- list(STARfusion_induced, arriba_induced, fcirc_induced, infusion_induced, fusioncatcher_induced, jaffa_d_induced, squid_induced,
                 STARfusion_known, arriba_known, fcirc_known, infusion_known, fusioncatcher_known, jaffa_d_known, squid_known,
                 STARfusion_public, arriba_public, fcirc_public, infusion_public, fusioncatcher_public, jaffa_d_public)

dataset_names <- c("STARfusion_induced", "arriba_induced", "fcirc_induced", "infusion_induced", "fusioncatcher_induced", "jaffa_d_induced", "squid_induced",
                   "STARfusion_known", "arriba_known", "fcirc_known", "infusion_known", "fusioncatcher_known", "jaffa_d_known", "squid_known",
                   "STARfusion_public", "arriba_public", "fcirc_public", "infusion_public", "fusioncatcher_public", "jaffa_d_public")

# Initialize empty all_data
all_data <- data.frame()

# Loop through each dataset and add the relevant columns
for (i in seq_along(datasets)) {
  if (nrow(datasets[[i]]) > 0) {
    current_dataset <- datasets[[i]]
    current_name <- dataset_names[i]
    
    # Add fusion, caller, and dataset_type columns
    current_dataset$fusion <- paste(current_dataset$gene1, current_dataset$gene2, sep = ":")
    current_dataset$caller <- gsub('(_induced|_known|_public)', '', current_name)
    current_dataset$dataset_type <- ifelse(grepl('_induced', current_name), 'induced', ifelse(grepl('_known', current_name), 'known', 'public'))
    
    # Filter based on expected fusions for the dataset type
    expected_fusions_current <- switch(current_dataset$dataset_type[1],
                                       "induced" = expected_fusions_induced,
                                       "known" = expected_fusions_known,
                                       "public" = expected_fusions_public)
    
    current_dataset <- current_dataset %>% filter(fusion %in% expected_fusions_current)
    
    # Add to all_data
    all_data <- bind_rows(all_data, current_dataset)
  }
}

# Keep only the relevant columns
all_data <- all_data %>% select(caller, dataset_type, cell_line, fusion)

# Initialize lists for cell_lines and expected_fusions
cell_lines_list <- list(
  "induced" = c("KIF5BE21d", "TFGC7_1d", "V1C1_1_d", "V3aA3_1d"),
  "known" = c('C2924hc1', 'C824hc1', 'C924hc1', 'YO_24hc1'),
  "public" = c("TCGA-78-7163","TCGA-86-A4P8","TCGA-50-8460","TCGA-67-6215","TCGA-05-4397","TCGA-97-8176","TCGA-97-7941","TCGA-05-4433","TCGA-05-4382","TCGA-55-7284","TCGA-50-8457","TCGA-50-6594","TCGA-50-5072","TCGA-55-A490","TCGA-55-8094","TCGA-97-8175","TCGA-L9-A5IP","TCGA-67-6217","TCGA-99-7458","TCGA-50-6593","TCGA-93-A4JN","TCGA-99-8028","TCGA-38-A44F","TCGA-55-A494","TCGA-53-A4EZ","TCGA-99-8025","TCGA-55-8089","TCGA-44-7661","TCGA-97-8174","TCGA-49-AAQV","TCGA-97-8179","TCGA-L9-A50W","TCGA-MN-A4N1","TCGA-MN-A4N5","TCGA-97-7937","TCGA-99-AA5R","TCGA-73-A9RS","TCGA-55-A48Z","TCGA-55-7913","TCGA-86-A4P8","TCGA-55-7914","TCGA-55-8513")
)

expected_fusions_list <- list(
  "induced" = c("KIF5B:ALK", "TFG:ALK", "EML4:ALK", "EML4:ALK"),
  "known" = rep("EML4:ALK", length(cell_lines_list[["known"]])),
  "public" = rep("EML4:ALK", length(cell_lines_list[["public"]]))
)
# Create an empty dataframe for the heatmap data
heatmap_data <- data.frame()

# Loop through each type of dataset (induced, known, public)
for (type in c("induced", "known", "public")) {
  # Filter all_data by dataset_type
  sub_data <- subset(all_data, dataset_type == type)
  
  # Get the corresponding cell lines and expected fusions from the lists
  current_cell_lines <- cell_lines_list[[type]]
  
  # Additional check to make sure only specified cell lines are included
  sub_data <- sub_data %>% filter(cell_line %in% current_cell_lines)
  
  current_expected_fusions <- expected_fusions_list[[type]]
  
  # Loop through each cell line
  for (cl in current_cell_lines) {
    # Loop through each expected fusion
    for (ef in current_expected_fusions) {
      # Loop through each caller
      for (caller in unique(sub_data$caller)) {
        # Check if the caller detected the expected fusion for this cell line
        detection_status <- ifelse(any(sub_data$caller == caller & sub_data$cell_line == cl & sub_data$fusion == ef), 1, 0)
        
        # Append to heatmap_data
        heatmap_data <- rbind(heatmap_data, data.frame(caller, type, cl, ef, detection_status))
      }
    }
  }
}

# Rename the columns for clarity
colnames(heatmap_data) <- c("caller", "dataset_type", "cell_line", "expected_fusion", "detection_status")


#Plot heat map
ggplot(heatmap_data, aes(x = caller, y = cell_line)) +
  geom_tile(aes(fill = factor(detection_status)), width = 0.9, height = 0.9) +
  scale_fill_manual(name = "Detection", values = c("0" = NA, "1" = '#4CAF50')) +
  facet_grid(dataset_type ~ ., scales = "free_y") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5))  # Adjust the size of font 
