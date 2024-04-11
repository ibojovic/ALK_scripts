# Known

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

#Import all databases 

STARfusion <- read.table("~/data/starfusion_output/STARfusion_out_public_ALK_data.txt", header = TRUE, sep = "\t")
arriba <- read.table("~/data/arriba_output/arriba_public_ALK_data.txt", header = TRUE, sep = "\t")
fcirc <- read.table("~/data/fcirc_output/fcirc_public_ALK_data.txt", header = TRUE, sep = "\t")
infusion <- read.table("~/data/infusion_output/infusion_public_ALK_data.txt", header = TRUE, sep = "\t")
fusioncatcher <- read.table("~/data/fusioncatcher_output/fusioncatcher_public_ALK_data.txt", header = TRUE, sep = "\t")
jaffa_d <- read.table("~/data/jaffa_output/jaffa_d_public_ALK_data.txt", header = TRUE, sep = "\t")



#change all column names in gene1 and gene2
colnames(STARfusion)[1:2] <- c("gene1", "gene2")
colnames(arriba)[1:2] <- c("gene1", "gene2")
colnames(fcirc)[1:2] <- c("gene1", "gene2")
colnames(infusion)[1:2] <- c("gene1", "gene2")
colnames(fusioncatcher)[1:2] <- c("gene1", "gene2")


#Genes in star has something like KIF5B^ENSG00000170759.11 (so i am getting rid of it)
STARfusion$gene1 <- gsub("\\^.*", "", STARfusion$gene1)
STARfusion$gene2 <- gsub("\\^.*", "", STARfusion$gene2)

colnames(STARfusion)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(fcirc)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(infusion)[3:4] <- c("breakpoint1", "breakpoint2")
colnames(fusioncatcher)[3:4] <- c("breakpoint1", "breakpoint2")
names(jaffa_d)[5]='breakpoint1'
names(jaffa_d)[7]='breakpoint2'


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


# Select only the required columns from each table
STARfusion <- STARfusion[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
arriba <- arriba[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
fcirc <- fcirc[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
infusion <- infusion[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
fusioncatcher <- fusioncatcher[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]
jaffa_d <- jaffa_d[, c("gene1", "gene2", "breakpoint1", "breakpoint2", "cell_line")]


# Convert breakpoint1 column to character in each table
STARfusion$breakpoint1 <- as.character(STARfusion$breakpoint1)
arriba$breakpoint1 <- as.character(arriba$breakpoint1)
fcirc$breakpoint1 <- as.character(fcirc$breakpoint1)
infusion$breakpoint1 <- as.character(infusion$breakpoint1)
fusioncatcher$breakpoint1 <- as.character(fusioncatcher$breakpoint1)
jaffa_d$breakpoint1 <- as.character(jaffa_d$breakpoint1)


# Convert breakpoint1 column to character in each table
STARfusion$breakpoint2 <- as.character(STARfusion$breakpoint2)
arriba$breakpoint2 <- as.character(arriba$breakpoint2)
fcirc$breakpoint2 <- as.character(fcirc$breakpoint2)
infusion$breakpoint2 <- as.character(infusion$breakpoint2)
fusioncatcher$breakpoint2 <- as.character(fusioncatcher$breakpoint2)
jaffa_d$breakpoint2 <- as.character(jaffa_d$breakpoint2)
STARfusion$gene1 <- as.character(STARfusion$gene1)
STARfusion$gene2 <- as.character(STARfusion$gene2)
arriba$gene1 <- as.character(arriba$gene1)
arriba$gene2 <- as.character(arriba$gene2)
fcirc$gene1 <- as.character(fcirc$gene1)
fcirc$gene2 <- as.character(fcirc$gene2)
infusion$gene1 <- as.character(infusion$gene1)
infusion$gene2 <- as.character(infusion$gene2)
fusioncatcher$gene1 <- as.character(fusioncatcher$gene1)
fusioncatcher$gene2 <- as.character(fusioncatcher$gene2)
jaffa_d$gene1 <- as.character(jaffa_d$gene1)
jaffa_d$gene2 <- as.character(jaffa_d$gene2)
STARfusion$cell_line <- as.character(STARfusion$cell_line)
arriba$cell_line <- as.character(arriba$cell_line)
fcirc$cell_line <- as.character(fcirc$cell_line)
infusion$cell_line <- as.character(infusion$cell_line)
fusioncatcher$cell_line <- as.character(fusioncatcher$cell_line)
jaffa_d$cell_line <- as.character(jaffa_d$cell_line)


# Create a list of table names
table_names <- c("STARfusion", "arriba", "fcirc", "infusion", "fusioncatcher", "jaffa_d")

# Combine all tables into a single dataframe with an additional column for table names
combined_df <- bind_rows(
  STARfusion %>% mutate(table_name = "STARfusion"),
  arriba %>% mutate(table_name = "arriba"),
  fcirc %>% mutate(table_name = "fcirc"),
  infusion %>% mutate(table_name = "infusion"),
  fusioncatcher %>% mutate(table_name = "fusioncatcher"),
  jaffa_d %>% mutate(table_name = "jaffa_d")
)

# Count the occurrences of each unique cell_line value for each table
counts <- combined_df %>%
  group_by(table_name, cell_line) %>%
  summarise(count = n()) %>%
  ungroup()

# Complete the dataset to include zeros
complete_counts <- counts %>%
  tidyr::complete(table_name, cell_line, fill=list(count=0)) %>%
  ungroup()

# Add a new column for log-transformed count
complete_counts$log_count <- log1p(complete_counts$count)


# Plot
p <- ggplot(complete_counts, aes(x = table_name, y = log_count, fill = cell_line)) +
  geom_bar(stat = "identity", position=position_dodge(width=0.9)) +
  theme_minimal() +
  labs(title = "Log Transformed Counts of Called Fusions by Caller", x = "Fusion Caller", y = "Log Transformed Count") +
  geom_text(aes(label=sprintf("%d", count), y = log_count + 0.1), position=position_dodge(width=0.9), color="black", size=3.5)

# Display the plot
print(p)


# Create a heatmap
ggplot(complete_counts, aes(x = table_name, y = cell_line)) +
  geom_tile(aes(fill = log_count), color = "white") +
  scale_fill_gradient(low = "white", high = "gold") + # You can change the color scheme
  theme_minimal() +
  labs(title = "Heatmap of Log Transformed Counts of Called Fusions by Caller",
       x = "Fusion Caller",
       y = "Cell Line",
       fill = "Log Transformed Count") +
  geom_text(aes(label = sprintf("%d", count)), vjust = 1) # To display count numbers
