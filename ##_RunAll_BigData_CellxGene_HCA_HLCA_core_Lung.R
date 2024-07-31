##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

#### Load Function ####
if(!require("devtools")) install.packages("devtools"); library(devtools)
if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
#open for Seurat Mult # trace("project_query", edit=TRUE) # layer = "data"
trace("project_query", edit=TRUE) # GetAssayData(new, "data") -> # GetAssayData(new, "RNA")
if(!require("VICTOR")) devtools::install_github("Charlene717/VICTOR"); library(VICTOR)

#### Set Parameter ####
GSE_Name = "HCA_HLCA_core"

#### Load Data Path ####
# Load a specific .RData file to get Actual_Cell_Type
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/CellxGene_HCA_HLCA_core_Lung/CellxGene_HCA_HLCA_core_Lung_Ref/CellxGene_HCA_HLCA_core_Lung_Ref_5000_Seed123.RData")
Actual_Cell_Type <- seuratObject@meta.data$Actual_Cell_Type %>% unique() %>% as.character()
Actual_Cell_Type <- c("None", Actual_Cell_Type)
rm(seuratObject)

# Define the folders for samples and references
Path_Sample_Folder <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/CellxGene_HCA_HLCA_core_Lung/CellxGene_HCA_HLCA_core_Lung_Query"
Path_Ref_Folder <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/CellxGene_HCA_HLCA_core_Lung/CellxGene_HCA_HLCA_core_Lung_Ref"

# Get the list of all .RData files from both folders
files_Sample <- list.files(Path_Sample_Folder, pattern = "\\.RData$", full.names = TRUE)
files_Ref <- list.files(Path_Ref_Folder, pattern = "\\.RData$", full.names = TRUE)

#### Run Main Code ####
# Helper function to extract matching part from file names
getMatchPart <- function(file_path, prefix) {
  sub(paste0(prefix, "_"), "", basename(file_path))
}

# Iterate through each Actual_Cell_Type
for (Cell_Type in Actual_Cell_Type) {
  # Setting Set_Ref_Delet according to Actual_Cell_Type
  Set_Ref_Delet <- Cell_Type

  # Iterate through each sample file
  for (sample_file in files_Sample) {
    # Extract Set_Sample from sample file name
    Set_Sample <- getMatchPart(sample_file, "CellxGene_HCA_HLCA_core_Lung_Query")

    # Iterate through each reference file
    for (ref_file in files_Ref) {
      # Extract Set_Reference from reference file name
      Set_Reference <- getMatchPart(ref_file, "CellxGene_HCA_HLCA_core_Lung_Ref_5000_Seed123_Delet")

      # Only proceed if the matching parts are identical
      if (Set_Sample == Set_Reference) {
        # Assign file paths
        Path_Sample <- sample_file
        Path_Ref <- ref_file

        # Source the R script
        try({
          source("##_RunAll_Main_Test.R")
        })

        # Clearing the variables except the ones needed for the next iteration
        rm(list = setdiff(ls(), c("GSE_Name", "ref_file", "sample_file", "Cell_Type",
                                  "Set_Sample", "Set_Reference",
                                  "files_Sample", "files_Ref", "Path_Sample_Folder", "Path_Ref_Folder",
                                  "getMatchPart", "Actual_Cell_Type", "Set_Ref_Delet")), envir = .GlobalEnv)
      }
    }
  }
}
