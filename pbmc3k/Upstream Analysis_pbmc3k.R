if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Seurat")

# LOAD DATA

library(Seurat)
library(dplyr)

getwd()
# setwd( "/Users/pengl/Downloads/Bioinformatics-RNA-seq/Test_Seurat")
#load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Downloads/Bioinformatics-RNA-seq/Test_Seurat/filtered_gene_bc_matrices/hg19/")

# will load the appropriate libraries and import the raw data. 
#We will also be optimizing memory usage (important when dealing with large datasets) 
#using seurat’s sparse matricesExamine the memory savings between regular and sparse matrices

dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size

sparse.size <- object.size((x = pbmc.data))
sparse.size

dense.size/sparse.size

packageVersion(pkg = "Seurat")

# initialize the Seurat object using raw "non-normalized" data
# perform some initial prefiltering of the data at cell level and gene level by excluding any genes 
# that do not expressed in at least 3 cells, and excluding any cells that do not have a minimum of 200
# expressed genes in total. name our project "10X_PBMC"
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200, project = "10X_PBMC")


# perform some initial QC on our cells
# here we do a filtering on user-defined outlier cells
# we visualize gene and molecule counts, plot their relationship, and exclude cells 
# with a clear outlier number of genes detected as potential multiplets. 
# nFeature_RNA (number of genes) and nCount_RNA are automatically calculated for every object by Seurat
# % of Count (i.g., UMI) mapping to MT-genes is a common scRNA-Seq QC metric.


# calculate percentage of mitochondrial genes using PercentageFeatureSet function
# If your mitochondrial genes are named differently, 
# then you will need to adjust this pattern accordingly (e.g. “mt-“, “mt.”, or “MT_” etc.).
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# explor the SeuratObject
colnames(x = pbmc)  # actually is a tag sequence for differenctiate cells
Cells(pbmc) # also are tag sequences for differenctiate cells
rownames(x = pbmc)  # is gene names
ncol(x = pbmc)  #2700
nrow(x = pbmc)  # 13714
Cells(pbmc)
Idents(object = pbmc) # just label the datasource of every cell, here is 10X_PBMC
levels(pbmc) # also indicates the source is 10X_PBMC

# INITIAL CELL LEVEL QC 
# visulize QC matric then use them to filter outlier later
# filter cells that have unique feature counts over 2500 or less than 200
# filter cells that have >5% mitochondiral counts

# Visualize QC metrics as a violin plot, 
# nFeature_RNA is number of nuclear genes, 
# nCounts_RNA is sum of the number of nuclear genes and mitochondiral genes
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
quartz()
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
CombinePlots(plots = list(plot1,plot2))

# based on graphS as above, we can define a window of minimum 200 deteced genes per cell
# and a maximum 2500 detected gene per cell
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# NORMALIZING DATA
# Seurat defaultly use a global-scaling normaization called as "LogNormalize" 
# gene expr of a cell divided by total gene expr of all calls, multiple scale factor of 10,000
# then do log transformation

# pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- NormalizeData((object = pbmc))

# SELECT HIGHLY VARIABLE GENES
# Seurat only uses the highly variable genes (subsets of features that exhibit high cell-to-cell
# variation), default 2000 features to do downstream anaylsis, like PCA
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
quartz()
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


#SCALE USING ScaleData function
# linear transformation, e.g, shifts mean expression of gene across cells is 0
# and scales the variance of each gene expression across cells is 1.
# scRNA-Seq data likely contains "uninteresting" sources of variation
# including tenchnical noise, batch effects, even biologica sources of variation (cell cycle stage)
# We can regress out these variation before dimensionatlity reduction and clustering

all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)

# PERFORM LINEAR DIMENTISONAL REDUCTION USING PCA on prevously scaled data
# by default only the previously determined high variable features are used as input
# but we can change the features's choice
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# examine and visulize PCA resutls a few different ways
# print the first 5 PCAs and thier corresponding 5 representative genes in each PCA
print(pbmc[["pca"]], dims = 1:5, nfeatures=5)
quartz()
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
quartz()
DimPlot(pbmc, reduction = "pca")
# DimHeatmap allows for easy exploration of the primary sources of heterogeneity and be useful when
# trying to decide which PCs to include for further downstream analyses.
quartz()
DimHeatmap(pbmc, dims=1, cells=500, balanced=TRUE)
quartz()
DimHeatmap(pbmc, dims=1:15, cells=500, balanced=TRUE)

# DETERMINE THE "dimensionality" of the dataset (i.e., how many dimensions to choose)
# use JackStraw do resampling, permute 1% of the data, rerun PCA, repeat
# then choose "significant" PCs as those who have a strong enrichmen of low p-value features
# NOTE: This process can take a long time for big datasets, comment out for expediency. 
# a approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time is long
# Significant’ PCs will show a strong enrichment of genes with low p-values 
# (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant.

# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)

# a more ad hoc method for determing which PCs to use is to look at a plot of the sd
# ElbowPlot function to visulize the ranking of PCs based on the percentage of variance explained by each one
# the figure suggest the majority of true signal is captured in the first 10 PCs
# we might choose anything between PC 7-12 as a cutoff
quartz()
ElbowPlot(pbmc)


# CLUSTER THE CELLS
# Applies a graph-based clustering approach, embedding cells in a graph structure for example KNN
# First construct a KNN graph, refine the edge weights.
# then apply Louvain algorithm (default) or SLM to optimize
# increse of "granularity" leading to a greater number of clusters
# 0.4-1.2 is usually good for datasets of around 3k. For larger datasets, increase the resolution
# The clusters can be found using the Idents function

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)  # optimizer

head(Idents(pbmc), 5)


# RUN NON-LINEAR dimensional reduction using(UMAP/tSNE)
# As input to UMAP/tSNE, suggest using the same PCs as input to the clustering analysis
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
quartz()
DimPlot(pbmc, reduction = "umap")
quartz()
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")

# SAVE the oject at this pint so that it can be easily loaded without haveing to rerun steps performed before
# or easily shared with collaborators

saveRDS(pbmc, file = "./pbmc_tutorial.rds")

# FINDING differentially expressed features (cluster biomarker)
# seurat automatically find postive and negative markers (specified in ident.1) of a single cluster
# FindAllMarkers do for all cluster as 1: the left comparasion
# min.pct requires a feature to be detected at a minimum percentage in either of the two groups
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1= 1, min.pct = 0.25)
head(cluster1.markers, n=5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Seurat has four tests for differential expression which can be set with the test.use parameter: 
# ROC test (“roc”), t-test (“t”), LRT test based on zero-inflated data (“bimod”, default), 
# LRT test based on tobit-censoring models (“tobit”) 
# The ROC test returns the ‘classification power’ for any individual marker (ranging from 0 – random, to 1 – perfect).
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5)

# VISULIZE THE IDENTIFIED MARKER EXPRESSION
# VlnPlot (shows expression probability distributions across clusters), 
# FeaturePlot (visualizes feature expression on a tSNE or PCA plot) are most commonly used visualizations. 
# We also suggest exploring RidgePlot, CellScatter, and DotPlot as additional methods to view your dataset.
quartz()
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
quartz()
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

quartz()
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))


# plotting the top 20 markers (or all markers if less than 20) for each cluster.
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
quartz()
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


#Assigning cell type identity to clusters
#Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
levels(pbmc)   # 0-8 at the present
names(new.cluster.ids)   # is null at present
names(new.cluster.ids) <- levels(pbmc)    # after assign with levels, form a tuple/dictionary pair
pbmc <- RenameIdents(pbmc, new.cluster.ids)
quartz()
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "./pbmc3k_final.rds")

