if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. Install it with BiocManager::install('DESeq2').")
}

count_matrix_file <- "GSE243183_TC-32_RawCountFile_RSEM_genes.txt"
sample_sheet_file <- "samples.txt"
count_table <- read.delim(count_matrix_file, check.names = FALSE)
samples <- read.delim(sample_sheet_file, stringsAsFactors = FALSE)
if (!"symbol" %in% names(count_table) ||
    !all(c("Sample", "Group") %in% names(samples))) {
  stop("Counts must contain symbol; metadata must contain Sample and Group.")
}
count_columns <- setdiff(names(count_table), "symbol")
if (!all(vapply(count_table[count_columns], is.numeric, logical(1)))) {
  stop("All sample count columns must be numeric.")
}
counts <- as.matrix(count_table[count_columns])
if (any(!is.finite(counts)) || any(counts < 0) ||
    anyDuplicated(count_table$symbol) || anyNA(count_table$symbol)) {
  stop("Counts must be finite and nonnegative, with unique nonmissing gene IDs.")
}
rownames(counts) <- count_table$symbol
if (anyDuplicated(samples$Sample) ||
    !all(colnames(counts) %in% samples$Sample)) {
  stop("Metadata must uniquely identify every count-matrix sample.")
}
samples <- samples[match(colnames(counts), samples$Sample), , drop = FALSE]
if (anyNA(samples$Group) || any(samples$Group == "")) {
  stop("Every sample must have a group.")
}
rownames(samples) <- samples$Sample
samples$Group <- factor(samples$Group)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts), colData = samples, design = ~ Group
)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
if (nrow(dds) < 2 || ncol(dds) < 3) {
  stop("QC requires at least two retained genes and three samples.")
}
dds <- DESeq2::estimateSizeFactors(dds)
transformed <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
expression_matrix <- SummarizedExperiment::assay(transformed)
stopifnot(all(is.finite(expression_matrix)))

group_colors <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"),
  c(levels(samples$Group), rep("", 6 - nlevels(samples$Group)))
)
sample_colors <- unname(group_colors[as.character(samples$Group)])
sample_labels <- sub("_S[0-9]+$", "", sub("^sample[0-9]+_48hr_", "", samples$Sample))

save_plot <- function(filename, draw) {
  png(filename, width = 2000, height = 1700, res = 200)
  on.exit(dev.off())
  draw()
  message("Saved ", filename)
}
group_legend <- function() {
  legend("topright", legend = levels(samples$Group),
         col = group_colors[levels(samples$Group)], pch = 16, bty = "n")
}

gene_variances <- apply(expression_matrix, 1, var)
top_genes <- head(order(gene_variances, decreasing = TRUE), 500)
top_genes <- top_genes[gene_variances[top_genes] > 0]
if (length(top_genes) < 2) {
  stop("PCA requires at least two genes with nonzero variance.")
}
pca <- prcomp(t(expression_matrix[top_genes, , drop = FALSE]), scale. = FALSE)
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
save_plot("QC_PCA.png", function() {
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(pca$x[, 1], pca$x[, 2], col = sample_colors, pch = 19, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%%)", variance_percent[1]),
       ylab = sprintf("PC2 (%.1f%%)", variance_percent[2]),
       main = "PCA: Most Variable Genes", ylim = extendrange(pca$x[, 2], f = 0.2),
       xlim = extendrange(pca$x[, 1], f = 0.2))
  text(pca$x[, 1], pca$x[, 2], labels = sample_labels, pos = 3, cex = 0.8)
  group_legend()
})

densities <- lapply(seq_len(ncol(expression_matrix)), function(sample_index) {
  density(expression_matrix[, sample_index])
})
save_plot("QC_density.png", function() {
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(NA, xlim = range(vapply(densities, function(curve) range(curve$x), numeric(2))),
       ylim = c(0, max(vapply(densities, function(curve) max(curve$y), numeric(1))) * 1.15),
       xlab = "Variance-Stabilized Expression", ylab = "Density",
       main = "Sample Expression Distributions")
  for (sample_index in seq_along(densities)) {
    lines(densities[[sample_index]], col = sample_colors[sample_index],
          lty = (sample_index - 1) %% 3 + 1, lwd = 2)
  }
  legend("topright", legend = sample_labels, col = sample_colors,
         lty = (seq_along(densities) - 1) %% 3 + 1, lwd = 2, bty = "n", cex = 0.8)
})

save_plot("QC_boxplot.png", function() {
  par(mar = c(8, 5, 4, 2) + 0.1)
  boxplot(as.data.frame(expression_matrix), names = sample_labels,
          col = sample_colors, las = 2, cex.axis = 0.85,
          ylab = "Variance-Stabilized Expression", main = "Sample Expression: Box and Whiskers",
          outline = TRUE, pch = 16, cex = 0.25)
  group_legend()
})

sample_correlation <- cor(expression_matrix, method = "pearson")
sample_order <- hclust(as.dist(1 - sample_correlation))$order
ordered_correlation <- sample_correlation[sample_order, sample_order]
save_plot("QC_correlation_heatmap.png", function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  par(mar = c(8, 8, 4, 1) + 0.1)
  color_limits <- range(sample_correlation)
  if (diff(color_limits) == 0) color_limits <- c(0, 1)
  heat_colors <- hcl.colors(100, "Viridis")
  sample_positions <- seq_len(ncol(expression_matrix))
  image(sample_positions, sample_positions, ordered_correlation,
        col = heat_colors, zlim = color_limits, axes = FALSE,
        xlab = "", ylab = "", main = "Sample Correlation (Pearson)", asp = 1)
  axis(1, at = sample_positions, labels = sample_labels[sample_order], las = 2)
  axis(2, at = sample_positions, labels = sample_labels[sample_order], las = 2)
  positions <- expand.grid(horizontal = sample_positions, vertical = sample_positions)
  text(positions$horizontal, positions$vertical,
       labels = sprintf("%.3f", as.vector(ordered_correlation)),
       col = ifelse(as.vector(ordered_correlation) > mean(color_limits), "black", "white"),
       cex = 0.85)
  par(mar = c(8, 1, 4, 4) + 0.1)
  color_values <- seq(color_limits[1], color_limits[2], length.out = 100)
  image(1, color_values, matrix(color_values, nrow = 1),
        col = heat_colors, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 2)
})

message("QC complete: ", nrow(expression_matrix), " genes across ",
        ncol(expression_matrix), " samples; blind VST with DESeq2 size-factor normalization.")