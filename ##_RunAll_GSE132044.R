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
trace("project_query", edit=TRUE) # layer = "data"
if(!require("VICTOR")) devtools::install_github("Charlene717/VICTOR"); library(VICTOR)

#### Set Parameter ####
GSE_Name = "GSE132044"

#### Load Data Path ####
# Load a specific .RData file to get Actual_Cell_Type
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Sample/GSE132044_10xV2_Sample_CellNum1663_Seed123_Processed.RData")
Actual_Cell_Type <- seuratObject@meta.data$Actual_Cell_Type %>% unique()
Actual_Cell_Type <- c("None", Actual_Cell_Type)
rm(seuratObject)

# Define the folders for samples and references
Path_Sample_Folder <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Sample_OriQC"
Path_Ref_Folder <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Ref_OriQC"

# Get the list of all .RData files from both folders
files_Sample <- list.files(Path_Sample_Folder, pattern = "\\.RData$", full.names = TRUE)
files_Ref <- list.files(Path_Ref_Folder, pattern = "\\.RData$", full.names = TRUE)



#### Run Main Code ####
# Helper function to extract Set_Sample and Set_Reference from file names
getSetName <- function(file_path) {
  parts <- unlist(strsplit(basename(file_path), "_", fixed = TRUE))
  return(paste(parts[1], parts[2], sep = "_"))
}

# Iterate through each Actual_Cell_Type
for (Cell_Type in Actual_Cell_Type) {
  # Setting Set_Ref_Delet according to Actual_Cell_Type
  Set_Ref_Delet <- Cell_Type

  # Iterate through each sample file
  for (sample_file in files_Sample) {
    # Extract Set_Sample from sample file name
    Set_Sample <- getSetName(sample_file)

    # Iterate through each reference file
    for (ref_file in files_Ref) {
      try({ detach("package:scPred", unload = TRUE) })
      # Extract Set_Reference from reference file name
      Set_Reference <- getSetName(ref_file)

      # Assign file paths
      Path_Sample <- sample_file
      Path_Ref <- ref_file

      # Source the R script
      try({
        # source("##_RunAll_CTAEvaluator_Main.R")
        source("##_RunAll_Main_Test.R")
      })

      # Clearing the variables except the ones needed for the next iteration
      rm(list = setdiff(ls(), c("ref_file", "sample_file", "Cell_Type",
                                "Set_Sample", "Set_Reference",
                                "files_Sample", "files_Ref", "Path_Sample_Folder", "Path_Ref_Folder",
                                "getSetName", "Actual_Cell_Type","Set_Ref_Delet")), envir = .GlobalEnv)
    }
  }

}
