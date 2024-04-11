# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr) # for expand()

# Read your data
induced_data <- readRDS("C:/Users/itana/Desktop/rstudio-export/combined_df_induced.rds")

# Specified cell lines and expected fusions
cell_lines_induced <- c("KIF5BE21d", "TFGC7_1d", "V1C1_1_d", "V3aA3_1d")
expected_fusions_induced <- c("KIF5B:ALK", "TFG:ALK", "EML4:ALK", "EML4:ALK")

# Create a new column 'fusion' which is a combination of gene1 and gene2
induced_data$fusion <- paste(induced_data$gene1, induced_data$gene2, sep=":")

# Create a new column 'expected' to tag whether the fusion is expected or not
induced_data$expected <- ifelse(induced_data$fusion %in% expected_fusions_induced & induced_data$cell_line %in% cell_lines_induced, "Expected", "Not Expected")

# Create a full grid of cell_line, table_name, and expected status
full_grid <- expand_grid(cell_line = c(unique(induced_data$cell_line), cell_lines_induced),
                         table_name = unique(induced_data$table_name),
                         expected = c("Expected", "Not Expected"))

# Find the "not called fusions" by anti-joining the full grid with the induced_data
not_called_fusions <- anti_join(full_grid,
                                induced_data %>% select(cell_line, table_name, expected),
                                by = c("cell_line", "table_name", "expected"))

# Add these not called fusions back into the data with n = 0
induced_data_summary <- induced_data %>%
  group_by(cell_line, expected, table_name) %>%
  summarize(n = n(), .groups = 'drop') %>%
  bind_rows(mutate(not_called_fusions, n = 0)) %>%
  group_by(cell_line, table_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(prop = ifelse(is.na(prop) & expected == "Not Expected", 1, prop))

# Create point_data to include missing cell_lines
original_point_data <- induced_data %>%
  group_by(cell_line, table_name) %>%
  summarise(detected = any(expected == "Expected"), total = n(), .groups = 'drop')

missing_point_data <- full_grid %>%
  anti_join(induced_data %>% select(cell_line, table_name), by = c("cell_line", "table_name")) %>%
  group_by(cell_line, table_name) %>%
  summarise(detected = FALSE, total = 0, .groups = 'drop')

point_data <- bind_rows(original_point_data, missing_point_data)

# Calculate the total number of unique cell lines
total_cell_lines <- length(unique(induced_data$cell_line))

# Calculate the average proportion for each table_name for "Expected" fusions
average_expected <- induced_data_summary %>%
  filter(expected == "Expected") %>%
  group_by(table_name) %>%
  summarize(avg_prop = sum(prop) / total_cell_lines, .groups = 'drop')

# Create the bar plot
p <- ggplot(induced_data_summary, aes(x = cell_line, y = prop, fill = expected)) +
  geom_bar(data = subset(induced_data_summary, expected != "Not Expected"), stat = "identity", position = "stack") +
  geom_segment(data = subset(induced_data_summary, expected == "Not Expected"), aes(x = as.numeric(cell_line) - 0.4, xend = as.numeric(cell_line) + 0.4, y = prop, yend = prop), color = "grey") +
  geom_text(data = point_data, aes(x = cell_line, y = 1, label = paste( total)), vjust = -0.5, inherit.aes = FALSE) +
  geom_point(data = point_data, aes(x = cell_line, y = 1, color = ifelse(detected, "Yes", "No"), size = log(total + 1)), inherit.aes = FALSE) +
  geom_segment(data = average_expected, aes(x = "Average", xend = "Average", y = 0, yend = avg_prop), arrow = arrow(type = "closed", length = unit(0.2, "inches")), color = "black", size = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Expected" = "lightseagreen")) +
  scale_color_manual(values = c("Yes" = "green", "No" = "red")) +
  ylim(0, 1) +
  ylab("Proportion") +
  xlab("Cell Line") +
  ggtitle("Proportion of Expected vs Not Expected Fusions") +
  facet_grid(~ table_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show the plot
p
