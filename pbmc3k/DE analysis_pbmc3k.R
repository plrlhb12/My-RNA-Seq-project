library(Seurat)

# retrieve the saved file which have gone through the upstream data analysis
pbmc3k_final <- readRDS("~/Downloads/Bioinformatics-RNA-seq/Test_Seurat/pbmc3k_final.rds")
pbmc3k_final

# Seurat use non-paramteric Wilcoxon rank sum test to calculate DE.
# the two groups to be tested will be speified as ident.1 and ident.2 parameters
levels(pbmc)
#"Naive CD4 T","Memory CD4 T", "CD14+ Mono", "B" , "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet"   
# find DE betwween two groups which arespeified as ident.1 and ident.2 parameters
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
# pct.1: percentage of the cells where the feature is detected in the first group
# pct.2: percentage of the cells where the feature is detected in the second group
# p_val_adj: adjusted p value based on bonferroni correction using all features in the dataset
head(monocyte.de.markers)
nrow(monocyte.de.markers)   #1380



# find DE between one group between all the other cells, only search for positive markers
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = NULL, only.pos = TRUE)
nrow(monocyte.de.markers)    #151

# To increase the speed of marker discover, use prefilter features or cells
# filter out the features that are either expressed at simillar levels or very infrequently detected in either groups of cells

# detection percentage: the frequency of a feature in all the cells of the same cluster 

# pre-filter features that are detected at <50% frequency in either A or B clusters
# min.pct = 
nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.pct = 0.5)) # 432

# Pre-filter features that have less than a two-fold change between the average expression of cluster1 vs cluster2
# logfc.threshold = log(2)
nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))   #333

# pre-filter features whose detection percentages across the two clusters are similar (difference of pct lessl than 0.25)

nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.diff.pct = 0.25))    #389

#Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing, 
#but could also miss features that are prefiltered
# if cell number in custers are large, Subsample each group to a maximum of 200 cells Can be very useful for computationally-intensive DE tests

nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", max.cells.per.ident = 200))   #1380

#Perform DE analysis using alternative tests
# wilcox (default), bimod, roc, poisson, t, LR, MAST, DESeq2
# MAST: GLM-framework the treates cellular detection rate as a covariate
# DESeq2: base on a model using the negative binomial distribution

library(DESeq2)
library(MAST)

# Test for DE features using the MAST package
nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "MAST"))      #1380
# test for DE features using DESeq2, but DESeq2 can be computationally intensive for large datasets
#  We therefore suggest initially limiting the number of cells used for testing
nrow(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "DESeq2", max.cells.per.ident = 50))     #3414

