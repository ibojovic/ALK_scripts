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

# Define the cell lines and expected fusions
cell_lines_induced <- c("KIF5BE21d", "TFGC7_1d", "V1C1_1_d", "V3aA3_1d")
cell_lines_known <- c('C2924hc1', 'C824hc1', 'C924hc1', 'YO_24hc1')
cell_lines_public <- c("TCGA-78-7163","TCGA-86-A4P8","TCGA-50-8460","TCGA-67-6215","TCGA-05-4397","TCGA-97-8176","TCGA-97-7941","TCGA-05-4433","TCGA-05-4382","TCGA-55-7284","TCGA-50-8457","TCGA-50-6594","TCGA-50-5072","TCGA-55-A490","TCGA-55-8094","TCGA-97-8175","TCGA-L9-A5IP","TCGA-67-6217","TCGA-99-7458","TCGA-50-6593","TCGA-93-A4JN","TCGA-99-8028","TCGA-38-A44F","TCGA-55-A494","TCGA-53-A4EZ","TCGA-99-8025","TCGA-55-8089","TCGA-44-7661","TCGA-97-8174","TCGA-49-AAQV","TCGA-97-8179","TCGA-L9-A50W","TCGA-MN-A4N1","TCGA-MN-A4N5","TCGA-97-7937","TCGA-99-AA5R","TCGA-73-A9RS","TCGA-55-A48Z","TCGA-55-7913","TCGA-86-A4P8","TCGA-55-7914","TCGA-55-8513")

expected_fusions_induced <- c("KIF5B:ALK", "TFG:ALK", "EML4:ALK", "EML4:ALK")
expected_fusions_known <- rep("EML4:ALK", length(cell_lines_known))
expected_fusions_public <- rep("EML4:ALK", length(cell_lines_public))

# Create a data frame named 'expected_fusions'
expected_fusions <- data.frame(
  cell_line = c(cell_lines_induced, cell_lines_known, cell_lines_public),
  expected_fusion = c(expected_fusions_induced, expected_fusions_known, expected_fusions_public)
)

# Create a function to determine the category based on the suffix of cell lines
get_category <- function(cell_line) {
  if (cell_line %in% cell_lines_induced) {
    return("induced")
  } else if (cell_line %in% cell_lines_known) {
    return("known")
  } else if (cell_line %in% cell_lines_public) {
    return("public")
  } else {
    return(NA)
  }
}

# Add the 'category' column based on the cell lines
expected_fusions$category <- sapply(expected_fusions$cell_line, get_category)


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

###############################################################################
###############################################################################
###############################################################################
# Select the necessary columns from each data frame
selected_columns <- c("gene1", "gene2", "cell_line")

# Add a 'fusion_caller' column to each data frame
STARfusion_induced_selected <- STARfusion_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "STAR-fusion")
arriba_induced_selected <- arriba_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Arriba")
fcirc_induced_selected <- fcirc_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "F-circ")
infusion_induced_selected <- infusion_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "InFusion")
fusioncatcher_induced_selected <- fusioncatcher_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "FusionCatcher")
jaffa_d_induced_selected <- jaffa_d_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Jaffa Direct")
squid_induced_selected <- squid_induced %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Squid")

STARfusion_known_selected <- STARfusion_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "STAR-fusion")
arriba_known_selected <- arriba_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Arriba")
fcirc_known_selected <- fcirc_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "F-circ")
infusion_known_selected <- infusion_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "InFusion")
fusioncatcher_known_selected <- fusioncatcher_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "FusionCatcher")
jaffa_d_known_selected <- jaffa_d_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Jaffa Direct")
squid_known_selected <- squid_known %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Squid")

STARfusion_public_selected <- STARfusion_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "STAR-fusion")
arriba_public_selected <- arriba_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Arriba")
fcirc_public_selected <- fcirc_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "F-circ")
infusion_public_selected <- infusion_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "InFusion")
fusioncatcher_public_selected <- fusioncatcher_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "FusionCatcher")
jaffa_d_public_selected <- jaffa_d_public %>% select(all_of(selected_columns)) %>% mutate(fusion_caller = "Jaffa Direct")

# Combine the selected columns from all data frames
combined_data <- bind_rows(STARfusion_induced_selected, arriba_induced_selected, fcirc_induced_selected, infusion_induced_selected, fusioncatcher_induced_selected, jaffa_d_induced_selected, squid_induced_selected, STARfusion_known_selected, arriba_known_selected, fcirc_known_selected, infusion_known_selected, fusioncatcher_known_selected, jaffa_d_known_selected, squid_known_selected, STARfusion_public_selected, arriba_public_selected, fcirc_public_selected, infusion_public_selected, fusioncatcher_public_selected, jaffa_d_public_selected)

# Add 'fusion_caller' and 'category' columns to each data frame
add_columns <- function(df, fusion_caller, category) {
  df %>% select(all_of(selected_columns)) %>% 
    mutate(fusion_caller = fusion_caller, category = category)
}

# Apply function to add columns
STARfusion_induced_selected <- add_columns(STARfusion_induced, "STAR-fusion", "induced")
arriba_induced_selected <- add_columns(arriba_induced, "Arriba", "induced")
fcirc_induced_selected <- add_columns(fcirc_induced, "F-circ", "induced")
infusion_induced_selected <- add_columns(infusion_induced, "InFusion", "induced")
fusioncatcher_induced_selected <- add_columns(fusioncatcher_induced, "FusionCatcher", "induced")
jaffa_d_induced_selected <- add_columns(jaffa_d_induced, "Jaffa Direct", "induced")
squid_induced_selected <- add_columns(squid_induced, "Squid", "induced")

STARfusion_known_selected <- add_columns(STARfusion_known, "STAR-fusion", "known")
arriba_known_selected <- add_columns(arriba_known, "Arriba", "known")
fcirc_known_selected <- add_columns(fcirc_known, "F-circ", "known")
infusion_known_selected <- add_columns(infusion_known, "InFusion", "known")
fusioncatcher_known_selected <- add_columns(fusioncatcher_known, "FusionCatcher", "known")
jaffa_d_known_selected <- add_columns(jaffa_d_known, "Jaffa Direct", "known")
squid_known_selected <- add_columns(squid_known, "Squid", "known")

STARfusion_public_selected <- add_columns(STARfusion_public, "STAR-fusion", "public")
arriba_public_selected <- add_columns(arriba_public, "Arriba", "public")
fcirc_public_selected <- add_columns(fcirc_public, "F-circ", "public")
infusion_public_selected <- add_columns(infusion_public, "InFusion", "public")
fusioncatcher_public_selected <- add_columns(fusioncatcher_public, "FusionCatcher", "public")
jaffa_d_public_selected <- add_columns(jaffa_d_public, "Jaffa Direct", "public")

# Combine the selected columns from all data frames
combined_data <- bind_rows(
  STARfusion_induced_selected, arriba_induced_selected, fcirc_induced_selected, 
  infusion_induced_selected, fusioncatcher_induced_selected, jaffa_d_induced_selected, 
  squid_induced_selected, STARfusion_known_selected, arriba_known_selected, fcirc_known_selected, 
  infusion_known_selected, fusioncatcher_known_selected, jaffa_d_known_selected, 
  squid_known_selected, STARfusion_public_selected, arriba_public_selected, fcirc_public_selected, 
  infusion_public_selected, fusioncatcher_public_selected, jaffa_d_public_selected
)

################################################################################
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Your data setup goes here
# combined_data and expected_fusions should be loaded here

# Count duplicates for each fusion_caller and category
duplicate_counts <- combined_data %>%
  group_by(fusion_caller, category, cell_line, gene1, gene2) %>%
  tally() %>%
  filter(n > 1) %>%
  summarise(total_duplicates = sum(n - 1), .groups = 'drop')

# Remove duplicates from combined_data
combined_data_unique <- combined_data %>%
  distinct(fusion_caller, category, cell_line, gene1, gene2, .keep_all = TRUE)

# Generate full grid for each category separately
full_grids <- lapply(unique(expected_fusions$category), function(cat) {
  expand_grid(
    fusion_caller = unique(combined_data$fusion_caller),
    cell_line = unique(expected_fusions$cell_line[expected_fusions$category == cat]),
    category = cat
  )
})

# Combine them back into a single data frame
full_grid <- bind_rows(full_grids)

# Merge the unique combined_data with expected_fusions
merged_data <- left_join(combined_data_unique, expected_fusions, by = c("cell_line", "category"))

# Create a new column to indicate whether the fusion is expected or not
merged_data$is_expected <- with(merged_data, paste0(gene1, ":", gene2) == expected_fusion)

# Find the "not called fusions" by anti-joining the full grid with the merged_data
not_called_fusions <- anti_join(full_grid, 
                                merged_data %>% select(fusion_caller, cell_line, category),
                                by = c("fusion_caller", "cell_line", "category"))

# Label them as "not expected"
not_called_fusions$is_expected <- FALSE

# Combine the merged_data with not_called_fusions
complete_data <- bind_rows(
  merged_data %>% select(fusion_caller, category, cell_line, is_expected),
  not_called_fusions
)

# Add the duplicates as unexpected
duplicate_counts$is_expected <- FALSE
complete_data <- bind_rows(
  complete_data,
  duplicate_counts %>% select(fusion_caller, category, is_expected)
)

# Group by fusion_caller, category, and is_expected, then calculate the count and proportion
summary_data <- complete_data %>% 
  group_by(fusion_caller, category, is_expected) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(fusion_caller, category) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Generate the plot
p <- ggplot(summary_data, aes(x = fusion_caller, y = proportion, fill = as.factor(is_expected))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "lightgreen", "FALSE" = "lightcoral")) +
  facet_wrap(~ category) +
  labs(title = "Proportion of Expected vs Unexpected Fusions",
       x = "Fusion Caller",
       y = "Proportion",
       fill = "Is Expected") +
  theme_minimal() +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.01)), vjust = -0.5)

print(p)

