##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)


#### Load data ####
load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")



