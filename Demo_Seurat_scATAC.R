## Ref: https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
# library(ggplot2); library(cowplot)

## GitHub BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("BiocGenerics", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors"))

if(!require("Signac")) install.packages("Signac"); library(Signac)
if(!require("EnsDb.Hsapiens.v86")) BiocManager::install("EnsDb.Hsapiens.v86"); library(EnsDb.Hsapiens.v86)
# if(!require("biovizBase")) BiocManager::install("biovizBase"); library(biovizBase)

## GitHub
if(!require("remotes")) install.packages("remotes"); library(remotes)
if(!require("SeuratData")) remotes::install_github("satijalab/seurat-data"); library(SeuratData)

library(SeuratData)
InstallData("pbmcMultiome")



#### Load Data ####
# library(SeuratData)
# InstallData("pbmcMultiome")
# library(Seurat)
# pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
# pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

# load both modalities
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

pbmc.rna[["RNA"]] <- as(pbmc.rna[["RNA"]], Class = "Assay5")
# repeat QC steps performed in the WNN vignette
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

# Perform standard analysis of each modality independently RNA analysis
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# ATAC analysis add gene annotation information
# BiocManager::install("biovizBase")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Plot
p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2

plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))

#### Export ####
save.image("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing_Plot.RData")
save(pbmc.rna, pbmc.atac,
     file = paste0("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData"))
