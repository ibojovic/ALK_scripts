# Load the necessary library
library(UpSetR)
library(tidyr)

# Define the cell lines to iterate through
cell_lines <- c("C2924hc1", "C824hc1", "C924hc1", "YO_24hc1")


#STAR-fusion
# Loop through each cell line and extract the required rows
STARfusion <- data.frame()
for (cell_line in cell_lines) {
  # Read the fusions data for the current cell line into a data frame
  fusions <- read.delim(paste0("~/temp/data/fusion_out/STARfusion_out/known_fusions/", cell_line, "/star-fusion.fusion_predictions.tsv"))
  
  # Extract the required columns from the data frame and add a new column for the cell line
  fusion_data <- data.frame(fusions[, c("LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint", "FFPM", "annots")])
  fusion_data$cell_line <- cell_line
  
  # Append the fusion data for the current cell line to the STARfusion data frame
  STARfusion <- rbind(STARfusion, fusion_data)
}



# Arriba
arriba <- data.frame()
for (cell_line in cell_lines) {
  # Read the fusions data for the current cell line into a data frame
  fusions <- read.delim(paste0("~/data/arriba_output/known_fusions/", cell_line, "/fusions.tsv"))
  
  # Extract the required columns from the data frame and add a new column for the cell line
  fusion_data <- data.frame(fusions[, c("X.gene1", "gene2", "breakpoint1", "breakpoint2", "type", "confidence", "coverage1", "coverage2")])
  fusion_data$cell_line <- cell_line
  
  # Append the fusion data for the current cell line to the arriba data frame
  arriba <- rbind(arriba, fusion_data)
}

# Fcirc
fcirc <- data.frame()

for (cell_line in cell_lines) {
  fusions <- read.delim(paste0("data/fcirc_output/known_fusions/", cell_line, "/fusion_results.tsv"))
  fusion_data <- data.frame(fusions[, c("X5.Gene", "X3.Gene", "X5.Gene_BreakPoint_Pos", "X3.Gene_BreakPoint_Pos", "BreakpointReads_Count", "P.Value")])
  fusion_data$cell_line <- cell_line
  fcirc <- rbind(fcirc, fusion_data)
}

# Infusion
infusion <- data.frame()

for (cell_line in cell_lines) {
  fusions <- read.delim(paste0("~/data/infusion_output/known_fusions/", cell_line, "/infusion_output/fusions.detailed.full.txt"))
  fusion_data <- data.frame(fusions[, c("gene_1", "gene_2", "break_pos1", "break_pos2", "pap_rate")])
  fusion_data$cell_line <- cell_line
  infusion <- rbind(infusion, fusion_data)
}

# FusionCatcher
fusioncatcher <- data.frame()

for (cell_line in cell_lines) {
  fusions <- read.delim(paste0("~/temp/data/fusion_out/fusioncatcher_out/known_fusions/", cell_line, "/final-list_candidate-fusion-genes.txt"))
  fusion_data <- data.frame(fusions[, c("Gene_1_symbol.5end_fusion_partner.", "Gene_2_symbol.3end_fusion_partner.", "Fusion_point_for_gene_1.5end_fusion_partner.", "Fusion_point_for_gene_2.3end_fusion_partner.", "Counts_of_common_mapping_reads")])
  fusion_data$cell_line <- cell_line
  fusioncatcher <- rbind(fusioncatcher, fusion_data)
}

# Jaffa Direct
jaffa_d <- data.frame()

for (cell_line in cell_lines) {
  fusions <- read.delim(paste0("~/data/jaffa_output/direct/known_fusions/", cell_line, "/jaffa_results.csv"), sep = ",")
  fusion_data <- data.frame(fusions[, c("fusion.genes", "chrom1", "base1", "chrom2", "base2", "classification", "known")])
  fusion_data <- separate(fusion_data, fusion.genes, into = c("gene1", "gene2"), sep = ":", remove = FALSE)
  fusion_data$cell_line <- cell_line
  jaffa_d <- rbind(jaffa_d, fusion_data)
}

# SQUID
squid <- data.frame()

for (cell_line in cell_lines) {
  fusions <- read.delim(paste0("~/temp/data/fusion_out/squid_out/known_fusions/", cell_line, "/", cell_line, "_squid_annotated.txt"))
  fusion_data <- data.frame(fusions[, c("FusedGenes", "start1", "end1", "start2", "end2", "score")])
  fusion_data$cell_line <- cell_line
  squid <- rbind(squid, fusion_data)
}

#change all column names in gene1 and gene2
colnames(STARfusion)[1:2] <- c("gene1", "gene2")
colnames(arriba)[1:2] <- c("gene1", "gene2")
colnames(fcirc)[1:2] <- c("gene1", "gene2")
colnames(infusion)[1:2] <- c("gene1", "gene2")
colnames(fusioncatcher)[1:2] <- c("gene1", "gene2")

squid$FusedGenes <- as.character(squid$FusedGenes)
squid$gene1 <- sapply(strsplit(squid$FusedGenes, ":"), `[`, 1)
squid$gene2 <- sapply(strsplit(squid$FusedGenes, ":"), `[`, 2)


#Genes in star has something like KIF5B^ENSG00000170759.11 (so i am getting rid of it)
STARfusion$gene1 <- gsub("\\^.*", "", STARfusion$gene1)
STARfusion$gene2 <- gsub("\\^.*", "", STARfusion$gene2)

colnames(STARfusion)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(fcirc)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(infusion)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(fusioncatcher)[3:4] <- c("breakpoint1", "breakpoint2")
names(jaffa_d)[5]='breakpoint1'
names(jaffa_d)[7]='breakpoint2'
names(squid)[5]='breakpoint1'
names(squid)[3]='breakpoint2'

arriba$breakpoint1 <- as.character(arriba$breakpoint1)
arriba$breakpoint2 <- as.character(arriba$breakpoint2)
arriba$breakpoint1 <- sapply(strsplit(arriba$breakpoint1, ":"), function(x) x[[2]])
arriba$breakpoint2 <- sapply(strsplit(arriba$breakpoint2, ":"), function(x) x[[2]])

STARfusion$breakpoint1 <- as.character(STARfusion$breakpoint1)
STARfusion$breakpoint2 <- as.character(STARfusion$breakpoint2)
STARfusion$breakpoint1 <- sub("^[^:]+:(.+):.+$", "\\1", STARfusion$breakpoint1)
STARfusion$breakpoint2 <- sub("^[^:]+:(.+):.+$", "\\1", STARfusion$breakpoint2)

fusioncatcher$breakpoint1 <- as.character(fusioncatcher$breakpoint1)
fusioncatcher$breakpoint2 <- as.character(fusioncatcher$breakpoint2)
fusioncatcher$breakpoint1 <- sub("^[^:]+:(.+):.+$", "\\1", fusioncatcher$breakpoint1)
fusioncatcher$breakpoint2 <- sub("^[^:]+:(.+):.+$", "\\1", fusioncatcher$breakpoint2)  


# Convert breakpoint1 column to character in each table
STARfusion$breakpoint1 <- as.character(STARfusion$breakpoint1)
arriba$breakpoint1 <- as.character(arriba$breakpoint1)
fcirc$breakpoint1 <- as.character(fcirc$breakpoint1)
infusion$breakpoint1 <- as.character(infusion$breakpoint1)
fusioncatcher$breakpoint1 <- as.character(fusioncatcher$breakpoint1)
jaffa_d$breakpoint1 <- as.character(jaffa_d$breakpoint1)
squid$breakpoint1 <- as.character(squid$breakpoint1)

# Convert breakpoint1 column to character in each table
STARfusion$breakpoint2 <- as.character(STARfusion$breakpoint2)
arriba$breakpoint2 <- as.character(arriba$breakpoint2)
fcirc$breakpoint2 <- as.character(fcirc$breakpoint2)
infusion$breakpoint2 <- as.character(infusion$breakpoint2)
fusioncatcher$breakpoint2 <- as.character(fusioncatcher$breakpoint2)
jaffa_d$breakpoint2 <- as.character(jaffa_d$breakpoint2)
squid$breakpoint2 <- as.character(squid$breakpoint2)  

# Add the tool name as a new column to each data frame
STARfusion$tool <- "STARfusion"
arriba$tool <- "Arriba"
fcirc$tool <- "Fcirc"
infusion$tool <- "Infusion"
fusioncatcher$tool <- "FusionCatcher"
jaffa_d$tool <- "Jaffa Direct"
squid$tool <- "SQUID"

# Select only the necessary columns
STARfusion <- STARfusion[, c("gene1", "gene2", "cell_line", "tool")]
arriba <- arriba[, c("gene1", "gene2", "cell_line", "tool")]
fcirc <- fcirc[, c("gene1", "gene2", "cell_line", "tool")]
infusion <- infusion[, c("gene1", "gene2", "cell_line", "tool")]
fusioncatcher <- fusioncatcher[, c("gene1", "gene2", "cell_line", "tool")]
jaffa_d <- jaffa_d[, c("gene1", "gene2", "cell_line", "tool")]
squid <- squid[, c("gene1", "gene2", "cell_line", "tool")]

# Deduplicate rows in each dataframe based on your columns
fcirc <- unique(fcirc[, c("gene1", "gene2", "cell_line")])
STARfusion <- unique(STARfusion[, c("gene1", "gene2", "cell_line")])
arriba <- unique(arriba[, c("gene1", "gene2", "cell_line")])
fusioncatcher <- unique(fusioncatcher[, c("gene1", "gene2", "cell_line")])
jaffa_d <- unique(jaffa_d[, c("gene1", "gene2", "cell_line")])
squid <- unique(squid[, c("gene1", "gene2", "cell_line")])
infusion <- unique(infusion[, c("gene1", "gene2", "cell_line")])


# Combine all the dataframes into a list
list_of_dfs <- list(fcirc = fcirc, STARfusion = STARfusion, arriba = arriba, 
                    fusioncatcher = fusioncatcher, jaffa_d = jaffa_d, squid = squid, infusion =infusion)
# First, create a unique identifier for each dataframe as before
list_of_dfs_with_id <- lapply(list_of_dfs, function(df) {
  df$unique_id <- paste(df$gene1, df$gene2, df$cell_line, sep = "_")
  return(df$unique_id)
})

# Get all unique identifiers
all_unique_ids <- unique(unlist(list_of_dfs_with_id))

# Create a binary matrix
binary_matrix <- matrix(0, nrow = length(all_unique_ids), ncol = length(list_of_dfs_with_id))
rownames(binary_matrix) <- all_unique_ids
colnames(binary_matrix) <- names(list_of_dfs_with_id)

# Fill in the binary matrix based on the presence of unique ids in each dataframe
for (i in seq_along(list_of_dfs_with_id)) {
  binary_matrix[rownames(binary_matrix) %in% list_of_dfs_with_id[[i]], i] <- 1
}

# Convert matrix to data.frame and then plot with UpSetR
binary_df <- as.data.frame(binary_matrix)

library(UpSetR)
upset(binary_df, sets = names(list_of_dfs_with_id), order.by = "degree")
