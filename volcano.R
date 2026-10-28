library(tidyverse)

input_file <- "DESeq2_results.tsv"
output_file <- "volcano.png"

results <- read_tsv(input_file, show_col_types = FALSE)
required_columns <- c("log2FoldChange", "padj")
if (!all(required_columns %in% names(results))) {
  stop("Input must contain log2FoldChange and padj columns.")
}
if (!all(vapply(results[required_columns], is.numeric, logical(1)))) {
  stop("log2FoldChange and padj must be numeric.")
}

results <- results %>%
  mutate(regulation = case_when(
    is.finite(log2FoldChange) & padj < 0.05 & log2FoldChange >= 2 ~ "Upregulated",
    is.finite(log2FoldChange) & padj < 0.05 & log2FoldChange <= 2 ~ "Downregulated",
    TRUE ~ "No change"
  ))

plot_data <- results %>%
  filter(is.finite(log2FoldChange), is.finite(padj), between(padj, 0, 1)) %>%
  mutate(neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin))) %>%
  arrange(regulation != "No change")
if (nrow(plot_data) == 0) {
  stop("No genes have finite log2FoldChange and valid padj values to plot.")
}

colors <- c(
  "Upregulated" = "#D55E00",
  "Downregulated" = "#0072B2",
  "No change" = "#999999"
)

volcano_plot <- ggplot(plot_data, aes(log2FoldChange, neg_log10_padj, color = regulation)) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#555555") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#555555") +
  scale_color_manual(values = colors, breaks = names(colors), drop = FALSE) +
  labs(title = "DESeq2 Volcano Plot", color = NULL,
       x = expression(log[2] ~ "fold change"),
       y = expression(-log[10] ~ "adjusted p-value")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

ggsave(output_file, plot = volcano_plot, width = 9, height = 7.5, dpi = 200)
message("Saved volcano plot to ", output_file)