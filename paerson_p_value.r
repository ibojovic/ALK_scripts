all_data <- readRDS("C:/Users/itana/Desktop/rstudio-export/all_data_expression_proportion.rds")


# Calculate Pearson correlation and p-values for each tool
correlation_results_with_p_value <- all_data %>%
  group_by(source) %>%
  do(cor_test = cor.test(.$log2_TPM, .$num_fusions, method = "pearson")) %>%
  rowwise() %>%
  mutate(correlation = cor_test$estimate,
         p_value = cor_test$p.value) %>%
  select(-cor_test)

print(correlation_results_with_p_value)
