# Seurat V5 Azimuth annotation
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://azimuth.hubmapconsortium.org/
#   Azimuth outputs Prediction scores, which range from 0 to 1 and reflect the confidence associated with each annotation. However, Azimuth does not provide a recommended threshold. Can I set it at 0.8 for our tests?

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)


#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna
