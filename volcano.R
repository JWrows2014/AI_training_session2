input_file <- "DESeq2_results.tsv"
output_file <- "volcano.png"

results <- read.delim(input_file, check.names = FALSE)
required_columns <- c("log2FoldChange", "padj")
if (!all(required_columns %in% names(results))) {
  stop("Input must contain log2FoldChange and padj columns.")
}
if (!all(vapply(results[required_columns], is.numeric, logical(1)))) {
  stop("log2FoldChange and padj must be numeric.")
}

results$regulation <- "No change"
significant <- !is.na(results$padj) & results$padj < 0.05 &
  is.finite(results$log2FoldChange)
results$regulation[significant & results$log2FoldChange <= 2] <- "Downregulated"
results$regulation[significant & results$log2FoldChange >= 2] <- "Upregulated"

valid <- is.finite(results$log2FoldChange) & is.finite(results$padj) &
  results$padj >= 0 & results$padj <= 1
plot_data <- results[valid, ]
if (nrow(plot_data) == 0) {
  stop("No genes have finite log2FoldChange and valid padj values to plot.")
}
plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
plot_data <- plot_data[order(plot_data$regulation != "No change"), ]

colors <- c(
  "Upregulated" = "#D55E00",
  "Downregulated" = "#0072B2",
  "No change" = "#999999"
)

png(output_file, width = 1800, height = 1500, res = 200)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(
  plot_data$log2FoldChange,
  plot_data$neg_log10_padj,
  col = adjustcolor(colors[plot_data$regulation], alpha.f = 0.65),
  pch = 16,
  cex = 0.65,
  xlab = expression(log[2] ~ "fold change"),
  ylab = expression(-log[10] ~ "adjusted p-value"),
  main = "DESeq2 Volcano Plot"
)
abline(v = 2, h = -log10(0.05), lty = 2, col = "#555555")
legend(
  "topright", legend = names(colors), col = colors,
  pch = 16, bty = "n", cex = 0.85
)
invisible(dev.off())
message("Saved volcano plot to ", output_file)