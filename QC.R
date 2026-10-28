library(tidyverse)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. Install it with BiocManager::install('DESeq2').")
}

count_matrix_file <- "GSE243183_TC-32_RawCountFile_RSEM_genes.txt"
sample_sheet_file <- "samples.txt"
count_table <- read_tsv(count_matrix_file, show_col_types = FALSE)
samples <- read_tsv(sample_sheet_file, show_col_types = FALSE)
if (!"symbol" %in% names(count_table) ||
    !all(c("Sample", "Group") %in% names(samples))) {
  stop("Counts must contain symbol; metadata must contain Sample and Group.")
}
count_columns <- setdiff(names(count_table), "symbol")
if (!all(map_lgl(select(count_table, all_of(count_columns)), is.numeric))) {
  stop("All sample count columns must be numeric.")
}
counts <- count_table %>% select(all_of(count_columns)) %>% as.matrix()
if (any(!is.finite(counts)) || any(counts < 0) ||
    anyDuplicated(count_table$symbol) || anyNA(count_table$symbol)) {
  stop("Counts must be finite and nonnegative, with unique nonmissing gene IDs.")
}
rownames(counts) <- count_table$symbol
if (anyDuplicated(samples$Sample) ||
    !all(colnames(counts) %in% samples$Sample)) {
  stop("Metadata must uniquely identify every count-matrix sample.")
}
samples <- tibble(Sample = colnames(counts)) %>% left_join(samples, by = "Sample")
if (anyNA(samples$Group) || any(samples$Group == "")) {
  stop("Every sample must have a group.")
}
samples <- samples %>%
  mutate(Group = factor(Group),
         SampleLabel = str_remove(str_remove(Sample, "^sample[0-9]+_48hr_"), "_S[0-9]+$")) %>%
  group_by(Group) %>%
  mutate(Replicate = factor(row_number())) %>%
  ungroup()

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = samples %>% column_to_rownames("Sample"), design = ~ Group
)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
if (nrow(dds) < 2 || ncol(dds) < 3) {
  stop("QC requires at least two retained genes and three samples.")
}
dds <- DESeq2::estimateSizeFactors(dds)
transformed <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
expression_matrix <- SummarizedExperiment::assay(transformed)
stopifnot(all(is.finite(expression_matrix)))

palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
if (nlevels(samples$Group) > length(palette)) {
  stop("The QC palette supports up to six sample groups.")
}
group_colors <- setNames(palette[seq_len(nlevels(samples$Group))], levels(samples$Group))
expression_long <- expression_matrix %>%
  as_tibble(rownames = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(samples, by = "Sample") %>%
  mutate(SampleLabel = factor(SampleLabel, levels = samples$SampleLabel))
qc_theme <- theme_classic(base_size = 12) + theme(legend.position = "top")

gene_variances <- apply(expression_matrix, 1, var)
top_genes <- head(order(gene_variances, decreasing = TRUE), 500)
top_genes <- top_genes[gene_variances[top_genes] > 0]
if (length(top_genes) < 2) {
  stop("PCA requires at least two genes with nonzero variance.")
}
pca <- prcomp(t(expression_matrix[top_genes, , drop = FALSE]), scale. = FALSE)
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_data <- pca$x %>% as_tibble(rownames = "Sample") %>% left_join(samples, by = "Sample")
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleLabel), vjust = -1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  labs(title = "PCA: Most Variable Genes",
       x = sprintf("PC1 (%.1f%%)", variance_percent[1]),
       y = sprintf("PC2 (%.1f%%)", variance_percent[2])) +
  qc_theme

density_plot <- ggplot(expression_long, aes(Expression, color = Group,
                                          linetype = Replicate, group = Sample)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = group_colors) +
  labs(title = "Sample Expression Distributions",
       x = "Variance-Stabilized Expression", y = "Density") +
  qc_theme

boxplot_plot <- ggplot(expression_long, aes(SampleLabel, Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Sample Expression: Box and Whiskers",
       x = NULL, y = "Variance-Stabilized Expression") +
  qc_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_correlation <- cor(expression_matrix, method = "pearson")
sample_order <- hclust(as.dist(1 - sample_correlation))$order
ordered_samples <- colnames(expression_matrix)[sample_order]
sample_labels <- setNames(samples$SampleLabel, samples$Sample)
color_limits <- range(sample_correlation)
if (diff(color_limits) == 0) color_limits <- c(0, 1)
correlation_long <- sample_correlation %>%
  as_tibble(rownames = "SampleRow") %>%
  pivot_longer(-SampleRow, names_to = "SampleColumn", values_to = "Correlation") %>%
  mutate(SampleRow = factor(SampleRow, levels = ordered_samples),
         SampleColumn = factor(SampleColumn, levels = ordered_samples),
         LabelColor = if_else(Correlation > mean(color_limits), "black", "white"))
correlation_plot <- ggplot(correlation_long, aes(SampleColumn, SampleRow, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", Correlation), color = LabelColor), size = 3.5) +
  scale_color_identity() +
  scale_fill_viridis_c(limits = color_limits, name = "Pearson r") +
  scale_x_discrete(labels = sample_labels, expand = c(0, 0)) +
  scale_y_discrete(labels = sample_labels, expand = c(0, 0)) +
  coord_equal() +
  labs(title = "Sample Correlation (Pearson)", x = NULL, y = NULL) +
  qc_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

qc_plots <- list(
  "QC_PCA.png" = pca_plot,
  "QC_density.png" = density_plot,
  "QC_boxplot.png" = boxplot_plot,
  "QC_correlation_heatmap.png" = correlation_plot
)
iwalk(qc_plots, function(plot_object, filename) {
  ggsave(filename, plot = plot_object, width = 10, height = 8.5, dpi = 200)
  message("Saved ", filename)
})

message("QC complete: ", nrow(expression_matrix), " genes across ",
        ncol(expression_matrix), " samples; blind VST with DESeq2 size-factor normalization.")