##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)


if(!require("VICTOR")) devtools::install_github("Charlene717/VICTOR"); library(VICTOR)

if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)


#### Load Data ####
