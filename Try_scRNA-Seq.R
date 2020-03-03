library(scRNAseq)
data(allen)
print(allen)
sce <- as(allen, "SingleCellExperiment")
sce
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
library(scater)
# to avoid the error message of "counts" not in names by package scater
# rename the tophat_counts to counts
counts(sce) <- assay(sce, "tophat_counts")
sce
sce_cellMetrix <-perCellQCMetrics(
  sce,
  subsets = NULL,
  percent_top = c(50, 100, 200, 500),
  BPPARAM = SerialParam(),
)
sce

#'calculateQCMetrics' is deprecated.
#Use 'perCellQCMetrics' or 'perFeatureQCMetrics' instead. 
cell_density_plot <- plot(density(sce_cellMetrix$total), main = "Density - total_counts")
threshold <- 50000
abline(v=threshold)


sce_featureMetrix <- perFeatureQCMetrics(sce)
sce_featureMetrix
plot(density(sce_featureMetrix$detected), main = "Density - total_features")
threshold <- 20000
abline(v=threshold)

keep <- (sce_cellMetrix$total > 100000)
table(keep)

browseVignettes("scRNAseq")

# gene filtering
filter_genes <- apply(counts(sce), 1, function(x){length(x[x>=2]) >= 2})

# most variable genes
BiocManager::install("magrittr")
library(magrittr)
vars <- assay(sce) %>% log1p %>% rowVars 
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE) 
head(vars)

sce_sub <- sce[names(vars[1:50]),] 
sce_sub


logcounts <- log1p(assay(sce_sub)) 
pca <- prcomp(t(logcounts))
reducedDims(sce_sub) <- SimpleList(PCA = pca$x) 
sce_sub

colData(sce)

pca <- reducedDim(sce_sub, "PCA")[, 1:2]
col <- colData(sce)[, c("Primary.Type", "Core.Type")]
df <- data.frame(pca, col)
ggplot(df, aes(x = PC1, y = PC2, col = Primary.Type, shape = Core.Type)) +
  geom_point()


