## Ref: https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette


library(SeuratData)
InstallData("pbmcMultiome")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

library(SeuratData)
library(Seurat)
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
