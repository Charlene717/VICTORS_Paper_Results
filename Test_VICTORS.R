
##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

Rec_Time_Point.lt <- list()
Rec_Time_Spend.lt <- list()

Rec_Time_Point.lt[["Start_Time"]] <- Sys.time() # %>% as.character()

#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

#### Load Dataset ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231118_PBMC3K/Export_PBMC3K_MislabelB/20231119014028BQVNGA_Multi/20231119014028BQVNGA.RData")
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/pbmc3k/pbmc3k_UniformCellTypeName_Reference/Seurat_pbmc3k_Seed123_Ref1350_Ref_UniformCTName.RData")
seuratObject_Ref <- seuratObject

load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/pbmc3k/pbmc3k_UniformCellTypeName_Sample/Seurat_pbmc3k_Seed123_Ref1350_Samp_UniformCTName.RData")
seuratObject_Sample <- seuratObject

Rec_Time_Point.lt[["Load_data"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Load_data"]] <- Rec_Time_Point.lt[["Load_data"]] - Rec_Time_Point.lt[["Start_Time"]]

#### VICTORS ####
source("FUN_VICTORSPred.R")
seuratObject_Ref <- getFeatureSpace(seuratObject_Ref, "Actual_Cell_Type")   ## Get the feature space to train the classifiers
Rec_Time_Point.lt[["getFeatureSpace"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["getFeatureSpace"]] <- Rec_Time_Point.lt[["getFeatureSpace"]] - Rec_Time_Point.lt[["Load_data"]]

seuratObject_Ref <- trainModel(seuratObject_Ref) ## Train the model
Rec_Time_Point.lt[["trainModel"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["trainModel"]] <- Rec_Time_Point.lt[["trainModel"]] - Rec_Time_Point.lt[["getFeatureSpace"]]

seuratObject_Sample <- VICTORSPred(seuratObject_Sample, seuratObject_Ref) #, threshold = Set_scPredict_Thr)
Rec_Time_Point.lt[["VICTORS"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["VICTORS"]] <- Rec_Time_Point.lt[["VICTORS"]] - Rec_Time_Point.lt[["trainModel"]]


colnames(seuratObject_Sample@meta.data) <- gsub("_plus", "+", colnames(seuratObject_Sample@meta.data))
View(seuratObject_Sample@meta.data)



