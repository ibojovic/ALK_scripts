# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr) # for expand()

# Read your data
induced_data <- readRDS("C:/Users/itana/Desktop/rstudio-export/combined_df_public.rds")



# Specified cell lines and expected fusions

# Specified cell lines and expected fusions
cell_lines_induced <-c("TCGA-78-7163","TCGA-86-A4P8","TCGA-50-8460","TCGA-67-6215","TCGA-05-4397","TCGA-97-8176","TCGA-97-7941","TCGA-05-4433","TCGA-05-4382","TCGA-55-7284","TCGA-50-8457","TCGA-50-6594","TCGA-50-5072","TCGA-55-A490","TCGA-55-8094","TCGA-97-8175","TCGA-L9-A5IP","TCGA-67-6217","TCGA-99-7458","TCGA-50-6593","TCGA-93-A4JN","TCGA-99-8028","TCGA-38-A44F","TCGA-55-A494","TCGA-53-A4EZ","TCGA-99-8025","TCGA-55-8089","TCGA-44-7661","TCGA-97-8174","TCGA-49-AAQV","TCGA-97-8179","TCGA-L9-A50W","TCGA-MN-A4N1","TCGA-MN-A4N5","TCGA-97-7937","TCGA-99-AA5R","TCGA-73-A9RS","TCGA-55-A48Z","TCGA-55-7913","TCGA-86-A4P8","TCGA-55-7914","TCGA-55-8513")
expected_fusions_induced <- rep("EML4:ALK", length(cell_lines_induced))  


# Create a new column 'fusion' which is a combination of gene1 and gene2
induced_data$fusion <- paste(induced_data$gene1, induced_data$gene2, sep=":")

# Tag whether the fusion is expected or not
induced_data <- induced_data %>%
  group_by(cell_line, table_name, fusion) %>%
  mutate(is_duplicate = row_number() > 1) %>%
  ungroup()

induced_data$expected <- ifelse(
  induced_data$fusion %in% expected_fusions_induced & 
    induced_data$cell_line %in% cell_lines_induced & 
    !induced_data$is_duplicate, 
  "Expected", 
  "Not Expected"
)

# Create a count table
count_table <- induced_data %>%
  group_by(cell_line, table_name, expected) %>%
  summarise(n = n(), .groups = 'drop')

# Generate a full grid for cell lines in cell_lines_induced, table_name, and expected status
full_grid <- expand_grid(cell_line = cell_lines_induced,
                         table_name = unique(induced_data$table_name),
                         expected = "Not Expected")

# Find the "not called fusions" by anti-joining the full grid with the count_table
not_called_fusions <- anti_join(full_grid, count_table, by = c("cell_line", "table_name", "expected"))

# Add these "not called fusions" back into the data with n = 0
count_table <- bind_rows(count_table, mutate(not_called_fusions, n = 0))

# Calculate proportions
count_table <- count_table %>%
  group_by(cell_line, table_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Create point_data to include missing cell_lines
point_data <- induced_data %>%
  group_by(cell_line, table_name) %>%
  summarise(detected = any(expected == "Expected"), total = n(), .groups = 'drop')

# Calculate the average proportion for each table_name for "Expected" fusions
average_expected <- count_table %>%
  filter(expected == "Expected") %>%
  group_by(table_name) %>%
  summarize(avg_prop = mean(prop), .groups = 'drop')

# Create the bar plot
p <- ggplot(count_table, aes(x = cell_line, y = prop, fill = expected)) +
  geom_bar(data = subset(count_table, expected != "Not Expected"), stat = "identity", position = "stack") +
  geom_segment(data = subset(count_table, expected == "Not Expected"), aes(x = cell_line, xend = cell_line, y = 1, yend = 1), color = "grey", size = 1) +
  geom_text(data = point_data, aes(x = cell_line, y = 1, label = paste(total)), vjust = -0.5, inherit.aes = FALSE) +
  geom_point(data = point_data, aes(x = cell_line, y = 1, color = ifelse(detected, "Yes", "No"), size = log(total + 1)), inherit.aes = FALSE) +
  #geom_segment(data = average_expected, aes(x = "Average", xend = "Average", y = 0, yend = avg_prop), arrow = arrow(type = "closed", length = unit(0.2, "inches")), color = "black", size = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Expected" = "violetred2")) +
  scale_color_manual(values = c("Yes" = "green", "No" = "red")) +
  ylim(0, 1) +
  ylab("Proportion") +
  xlab("Cell Line") +
  ggtitle("Proportion of Expected vs Not Expected Fusions") +
  facet_grid(~ table_name) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Adjust the size here
    axis.text.y = element_text(size = 10)  # Adjust the size here
  )
# Show the plot
p













