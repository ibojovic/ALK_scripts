# INDUCED

library(ggplot2)
library(reshape2)
#Import all databases 

STARfusion <- read.table("~/data/starfusion_output/STARfusion_out_induced_ALK_data.txt", header = TRUE, sep = "\t")
arriba <- read.table("~/data/arriba_output/arriba_induced_ALK_data.txt", header = TRUE, sep = "\t")
fcirc <- read.table("~/data/fcirc_output/fcirc_induced_ALK_data.txt", header = TRUE, sep = "\t")
infusion <- read.table("~/data/infusion_output/infusion_induced_ALK_data.txt", header = TRUE, sep = "\t")
fusioncatcher <- read.table("~/data/fusioncatcher_output/fusioncatcher_induced_ALK_data.txt", header = TRUE, sep = "\t")
jaffa_d <- read.table("~/data/jaffa_output/jaffa_d_induced_ALK_data.txt", header = TRUE, sep = "\t")
squid<-read.table("~/data/squid_output/squid_induced_ALK_data.txt", header = TRUE, sep = "\t") 

starfusion_palette <- c("#0E515A", "#007880", "#00A9A9", "#32C5D2", "#80CED7")
fcirc_palette <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")
infusion_palette <- c("#006633", "#339966", "#66cccc", "#99cc99", "#ccff66")
fusioncatcher_palette <- c("#4B0082", "#483D8B", "#6A5ACD", "#7B68EE", "#9370DB")
arriba_palette <- c("#0072BB", "#0099D5", "#8AC7DE", "#D1E7F5", "#F8D7B5", "#E18D96", "#B53920")
jaffa_palette <- c("#FFD700", "#F0E68C", "#FFDAB9", "#FFA500", "#FFFF00")
squid_palette<- c("#0E1447", "#2E4B8C", "#508CD7", "#8CC7FF", "#C5E3FF")

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

# create a dataframe with cell lines, expected fusions, and breakpoints
cell_lines <- c("KIF5BE21d", "TFGC7_1d", "V1C1_1_d", "V3aA3_1d")
# Create a list of expected fusions for each sample
expected_fusions <- c("KIF5B:ALK", "TFG:ALK", "EML4:ALK", "EML4:ALK")
gene1 <- c("KIF5B", "TFG", "EML4", "EML4")
gene2 <- c("ALK", "ALK", "ALK", "ALK")
breakpoint1 <- c(32017143, 100728858, 42295516, 42264731)
breakpoint2 <- c(29223528, 29223528, 29223528, 29223528)

# Select only the required columns from each table
STARfusion <- STARfusion[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
arriba <- arriba[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
fcirc <- fcirc[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
infusion <- infusion[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
fusioncatcher <- fusioncatcher[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
jaffa_d <- jaffa_d[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
squid <- squid[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]

library(ggplot2)
library(dplyr)

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

library(ggplot2)
library(dplyr)
library(tidyr)

# Create a list of table names
table_names <- c("STARfusion", "arriba", "fcirc", "infusion", "fusioncatcher", "jaffa_d", "squid")

# Combine all tables into a single dataframe with an additional column for table names
combined_df <- bind_rows(
  STARfusion %>% mutate(table_name = "STARfusion"),
  arriba %>% mutate(table_name = "arriba"),
  fcirc %>% mutate(table_name = "fcirc"),
  infusion %>% mutate(table_name = "infusion"),
  fusioncatcher %>% mutate(table_name = "fusioncatcher"),
  jaffa_d %>% mutate(table_name = "jaffa_d"),
  squid %>% mutate(table_name = "squid")
)

# Count the occurrences of each unique cell_line value for each table
counts <- combined_df %>%
  group_by(table_name, cell_line) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the counts for each cell_line value within each table
ggplot(counts, aes(x = table_name, y = count, fill = cell_line)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Table", y = "Count", fill = "Cell Line") +
  ggtitle("Counts of Called Fusions for each Fusion Caller") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



######## LOGARITHMIC SCALE #######

# library(ggplot2)
# library(dplyr)
# library(scales)
# 
# # Create a list of table names
# table_names <- c("STARfusion", "arriba", "fcirc", "infusion", "fusioncatcher", "jaffa_d", "squid")
# 
# # Combine all tables into a single dataframe with an additional column for table names
# combined_df <- bind_rows(
#   STARfusion %>% mutate(table_name = "STARfusion"),
#   arriba %>% mutate(table_name = "arriba"),
#   fcirc %>% mutate(table_name = "fcirc"),
#   infusion %>% mutate(table_name = "infusion"),
#   fusioncatcher %>% mutate(table_name = "fusioncatcher"),
#   jaffa_d %>% mutate(table_name = "jaffa_d"),
#   squid %>% mutate(table_name = "squid")
# )
# 
# # Count the occurrences of each unique cell_line value for each table
# counts <- combined_df %>%
#   group_by(table_name, cell_line) %>%
#   summarise(count = n()) %>%
#   ungroup()
# 
# # Apply logarithmic scaling to the count values
# counts$log_count <- log10(counts$count)
# 
# # Plot the counts for each cell_line value within each table
# ggplot(counts, aes(x = table_name, y = log_count, fill = cell_line)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Table", y = "Logarithmic Count", fill = "Cell Line") +
#   ggtitle("Counts of Called Fusions for each Fusion Caller") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_y_log10()  # Apply logarithmic scaling to the y-axis
