library(dplyr)
library(Seurat)
library(patchwork)

setwd("/media/gbim/bigdata/2025-40")

### Sample 2025A400 ###
dataset2025A400 <- Read10X(data.dir = "/media/gbim/bigdata/2025-40/2025A400/outs/raw_feature_bc_matrix")
sample2025A400 <- CreateSeuratObject(counts = dataset2025A400, project = "2025A400", min.cells = 3, min.features = 200)
sample2025A400

mt_genes_all <- c("ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
                  "COX1", "COX2", "COX3",
                  "ATP6", "ATP8", "CYTB")

mt_genes <- intersect(mt_genes_all, rownames(sample2025A400))
sample2025A400[["percent.mt"]] <- PercentageFeatureSet(sample2025A400, features = mt_genes)

#sample2025A400[["percent.mt"]] <- PercentageFeatureSet(sample2025A400, pattern = "^(?i)MT")
VlnPlot(sample2025A400, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(sample2025A400, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample2025A400, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


### Sample 2025A401 ###
dataset2025A401 <- Read10X(data.dir = "/media/gbim/bigdata/2025-40/2025A401/outs/raw_feature_bc_matrix")
sample2025A401 <- CreateSeuratObject(counts = dataset2025A401, project = "2025A401", min.cells = 3, min.features = 200)
sample2025A401

mt_genes_all <- c("ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
                  "COX1", "COX2", "COX3",
                  "ATP6", "ATP8", "CYTB")

mt_genes <- intersect(mt_genes_all, rownames(sample2025A401))
sample2025A401[["percent.mt"]] <- PercentageFeatureSet(sample2025A401, features = mt_genes)
VlnPlot(sample2025A401, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


### Sample 2025A402 ###
dataset2025A402 <- Read10X(data.dir = "/media/gbim/bigdata/2025-40/2025A402/outs/raw_feature_bc_matrix")
sample2025A402 <- CreateSeuratObject(counts = dataset2025A402, project = "2025A402", min.cells = 3, min.features = 200)
sample2025A402

mt_genes_all <- c("ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
                  "COX1", "COX2", "COX3",
                  "ATP6", "ATP8", "CYTB")

mt_genes <- intersect(mt_genes_all, rownames(sample2025A402))
sample2025A402[["percent.mt"]] <- PercentageFeatureSet(sample2025A402, features = mt_genes)
VlnPlot(sample2025A402, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Sample 2025A403 ###
dataset2025A403 <- Read10X(data.dir = "/media/gbim/bigdata/2025-40/2025A403/outs/raw_feature_bc_matrix")
sample2025A403 <- CreateSeuratObject(counts = dataset2025A403, project = "2025A403", min.cells = 3, min.features = 200)
sample2025A403

mt_genes_all <- c("ND1", "ND2", "ND3", "ND4L", "ND4", "ND5", "ND6",
                  "COX1", "COX2", "COX3",
                  "ATP6", "ATP8", "CYTB")

mt_genes <- intersect(mt_genes_all, rownames(sample2025A403))
sample2025A403[["percent.mt"]] <- PercentageFeatureSet(sample2025A403, features = mt_genes)
VlnPlot(sample2025A403, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
