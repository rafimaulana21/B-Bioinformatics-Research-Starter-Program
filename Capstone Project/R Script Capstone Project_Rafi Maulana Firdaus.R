################################################################
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
# Bioinformatic Analysis of pDC Response to Influenza A (GSE68849)
# Pipeline: Load -> Preprocess -> QC -> Paired Limma DEG -> Viz -> Enrichment
################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(enrichR)
})

# ==============================================================================
# 0. SETUP OUTPUT FOLDERS
# ==============================================================================
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
gset_list <- getGEO("GSE68849", GSEMatrix = TRUE, AnnotGPL = TRUE)

if (length(gset_list) > 1) {
  idx <- grep("GPL10558", attr(gset_list, "names"))
  if (length(idx) == 0) idx <- 1
} else {
  idx <- 1
}
gset <- gset_list[[idx]]

# Clean feature labels
fvarLabels(gset) <- make.names(fvarLabels(gset))

# ==============================================================================
# 2. PREPROCESSING (LOG2 CHECK + FILTER)
# ==============================================================================
ex <- exprs(gset)

# log2 transformation check (GEO2R style)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogC) {
  ex[ex <= 0] <- NA
  exprs(gset) <- log2(ex)
}

# Remove rows with NA after log
gset <- gset[complete.cases(exprs(gset)), ]
ex <- exprs(gset)  # IMPORTANT: refresh expression matrix after preprocessing

# ==============================================================================
# 3. DEFINE PHENOTYPE (Condition + Donor) FROM METADATA
# ==============================================================================
pd <- pData(gset)

# Choose a text field to detect condition
# Try commonly informative fields: title / characteristics_ch1 / source_name_ch1
meta_text <- rep("", nrow(pd))
if ("title" %in% colnames(pd)) meta_text <- paste(meta_text, pd$title)
if ("source_name_ch1" %in% colnames(pd)) meta_text <- paste(meta_text, pd$source_name_ch1)

# Add all "characteristics_ch1" columns if present
char_cols <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
if (length(char_cols) > 0) {
  meta_text <- paste(meta_text, apply(pd[, char_cols, drop = FALSE], 1, paste, collapse = " "))
}

meta_text <- tolower(meta_text)

# Detect Control vs Influenza
is_control <- grepl("control|no virus|mock", meta_text)
is_flu     <- grepl("influenza|flu|h1n1|a/pr/8|pr8|infect", meta_text)

# If some samples match both/none, fall back to simpler rule:
# prioritize influenza keywords, else control.
Condition <- ifelse(is_flu & !is_control, "Influenza",
                    ifelse(is_control & !is_flu, "Control",
                           ifelse(is_flu, "Influenza", "Control")))

Condition <- factor(Condition, levels = c("Control", "Influenza"))

# Donor inference:
# If donor info exists in metadata, use it. Otherwise assume paired ordering (1..5 each=2).
donor_col <- NULL
possible_donor_cols <- c("donor", "subject", "patient", "individual", "donor_id")
for (cname in colnames(pd)) {
  if (any(grepl(paste(possible_donor_cols, collapse="|"), cname, ignore.case = TRUE))) {
    donor_col <- cname
    break
  }
}

if (!is.null(donor_col)) {
  Donor <- factor(pd[[donor_col]])
} else {
  # Fallback (common for this dataset): 5 donors x 2 conditions
  # IMPORTANT: this assumes sample order is paired; if not, reorder using metadata.
  Donor <- factor(rep(1:5, each = 2))
}

gset$Condition <- Condition
gset$Donor <- Donor

# Save phenotype table for reproducibility
pheno_out <- data.frame(
  sample = rownames(pd),
  Condition = gset$Condition,
  Donor = gset$Donor,
  meta_text = meta_text
)
write.csv(pheno_out, "results/phenotype_table.csv", row.names = FALSE)

# ==============================================================================
# 4. QUALITY CONTROL (QC) VISUALIZATIONS
# ==============================================================================
palette_colors <- c(Control = "#2ca25f", Influenza = "#e34a33")

# A. Boxplot
png("figures/Fig1_Boxplot.png", width = 1400, height = 900, res = 180)
par(mar = c(8, 4, 3, 1))
boxplot(ex, boxwex = 0.6, notch = TRUE, outline = FALSE, las = 2,
        col = palette_colors[as.character(Condition)],
        main = "Boxplot: Distribusi Ekspresi Gen Antar Sampel",
        ylab = "Log2 Expression")
legend("topleft", legend = levels(Condition), fill = palette_colors, bty = "n")
dev.off()

# B. Density plot
png("figures/Fig2_DensityPlot.png", width = 1400, height = 900, res = 180)
par(mar = c(5, 4, 3, 1))
plotDensities(ex, group = Condition, col = palette_colors[levels(Condition)],
              main = "Density Plot: Sebaran Ekspresi per Kelompok", legend = "topright")
dev.off()

# ==============================================================================
# 5. DIFFERENTIAL EXPRESSION ANALYSIS (PAIRED LIMMA)
# ==============================================================================
design <- model.matrix(~ Donor + Condition)
fit <- lmFit(gset, design)
fit <- eBayes(fit)

tT_full <- topTable(fit, coef = "ConditionInfluenza", adjust = "fdr",
                    sort.by = "B", number = Inf)

# Keep safe columns if present
keep_cols <- c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")
tT <- tT_full[, intersect(keep_cols, colnames(tT_full)), drop = FALSE]

# Define significance
tT$Significance <- "Not Significant"
tT$Significance[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "Up-regulated"
tT$Significance[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "Down-regulated"

write.csv(tT, "results/DEG_full_table.csv", row.names = TRUE)

# ==============================================================================
# 6. VISUALIZATION OF RESULTS
# ==============================================================================
# A. PCA Plot
pca <- prcomp(t(ex), scale. = TRUE)
pca_data <- data.frame(pca$x, Condition = Condition, Donor = Donor)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Donor)) +
  geom_point(size = 4, alpha = 0.85) +
  theme_minimal() +
  labs(title = "PCA: Sample Clustering (Control vs Influenza)",
       x = "PC1", y = "PC2") +
  scale_color_manual(values = palette_colors)

ggsave("figures/Fig3_PCA.png", pca_plot, width = 7, height = 5, dpi = 300)

# B. Volcano Plot (use FDR on y-axis for consistency)
top_genes_volcano <- head(tT[order(tT$adj.P.Val), ], 15)

volcano_plot <- ggplot(tT, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.6) +
  scale_color_manual(values = c("Down-regulated" = "#2b8cbe",
                                "Not Significant" = "#bdbdbd",
                                "Up-regulated" = "#f03b20")) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", alpha = 0.5) +
  geom_text_repel(data = top_genes_volcano,
                  aes(label = Gene.symbol),
                  color = "black", size = 3.3, max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Volcano Plot: Influenza vs Control (Paired, FDR-based)",
       x = "log2 Fold Change", y = "-log10(FDR)")

ggsave("figures/Fig4_Volcano.png", volcano_plot, width = 7, height = 5, dpi = 300)

# C. Heatmap (Top 50 by FDR)
top_50_ids <- rownames(head(tT[order(tT$adj.P.Val), ], 50))
mat_top50 <- ex[top_50_ids, , drop = FALSE]

# Replace rownames with gene symbols if available
if ("Gene.symbol" %in% colnames(tT)) {
  gene_symbols <- tT[top_50_ids, "Gene.symbol"]
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- top_50_ids[is.na(gene_symbols) | gene_symbols == ""]
  # Ensure unique rownames for heatmap
  rownames(mat_top50) <- make.unique(gene_symbols)
}

annotation_col <- data.frame(Condition = Condition)
rownames(annotation_col) <- colnames(mat_top50)

ann_colors <- list(Condition = palette_colors)

png("figures/Fig5_Heatmap_Top50.png", width = 1600, height = 1200, res = 200)
pheatmap(mat_top50,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         show_colnames = FALSE,
         fontsize_row = 7,
         main = "Heatmap: Top 50 Differentially Expressed Genes")
dev.off()

# ==============================================================================
# 7. FUNCTIONAL ENRICHMENT (EnrichR) + DOTPLOTS
# ==============================================================================
# Use up-regulated genes (clean gene symbols)
if (!("Gene.symbol" %in% colnames(tT))) {
  stop("Kolom Gene.symbol tidak ditemukan. Pastikan AnnotGPL=TRUE dan platform annotation tersedia.")
}

sig_up <- subset(tT, Significance == "Up-regulated" & !is.na(Gene.symbol) & Gene.symbol != "")
genes_up <- sig_up$Gene.symbol
genes_up <- unlist(strsplit(genes_up, " /// |;|,"))
genes_up <- unique(trimws(genes_up))
genes_up <- genes_up[genes_up != ""]

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
enr <- enrichr(genes_up, dbs)
saveRDS(enr, "results/enrichr_upregulated.rds")

# --- GO DOTPLOT (Top 10) ---
go_res <- head(enr$GO_Biological_Process_2023, 10)
go_res$GeneCount <- as.numeric(sub("/.*", "", go_res$Overlap))
go_res$Term_Clean <- sub(" \\(GO:.*\\)", "", go_res$Term)

dotplot_go <- ggplot(go_res, aes(x = GeneCount, y = reorder(Term_Clean, GeneCount))) +
  geom_point(aes(size = GeneCount, color = Adjusted.P.value)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Top 10 GO Biological Process (Up-regulated)",
       x = "Number of Genes", y = "GO Term", color = "Adj. P-value", size = "Gene Count") +
  theme(axis.text.y = element_text(size = 10))

ggsave("figures/Fig6_GO_Dotplot.png", dotplot_go, width = 8, height = 5.5, dpi = 300)

# --- KEGG DOTPLOT (Top 10) ---
kegg_res <- head(enr$KEGG_2021_Human, 10)
kegg_res$GeneCount <- as.numeric(sub("/.*", "", kegg_res$Overlap))

dotplot_kegg <- ggplot(kegg_res, aes(x = GeneCount, y = reorder(Term, GeneCount))) +
  geom_point(aes(size = GeneCount, color = Adjusted.P.value)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Top 10 KEGG Pathways (Up-regulated)",
       x = "Number of Genes", y = "KEGG Pathway", color = "Adj. P-value", size = "Gene Count") +
  theme(axis.text.y = element_text(size = 10))

ggsave("figures/Fig7_KEGG_Dotplot.png", dotplot_kegg, width = 8, height = 5.5, dpi = 300)

# ==============================================================================
# 8. QUICK SUMMARY PRINTS
# ==============================================================================
cat("\n===== SUMMARY =====\n")
cat("Samples:", ncol(ex), "\n")
cat("Genes (after filtering):", nrow(ex), "\n")
cat("Condition counts:\n")
print(table(Condition))
cat("\nDEGs (FDR<0.05 & |logFC|>1):\n")
print(table(tT$Significance))
cat("\nTop 10 Up-regulated genes (by FDR):\n")
print(head(tT[tT$Significance == "Up-regulated", c("Gene.symbol", "logFC", "adj.P.Val")][order(tT[tT$Significance=="Up-regulated","adj.P.Val"]), ], 10))
cat("\nOutputs saved to: results/ and figures/\n")

# Cek Working Directory
getwd()
