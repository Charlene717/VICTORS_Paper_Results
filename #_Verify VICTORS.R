## Ref: https://satijalab.org/seurat/reference/addmodulescore
## Ref: https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/
## Paper: https://www.science.org/doi/10.1126/science.aad0501
### ChatGPT Record: https://chat.openai.com/share/58ed3d0f-f52f-4643-8ce7-8d2b4997829a

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

## Record time set
Rec_Time_Point.lt <- list()
Rec_Time_Spend.lt <- list()

Rec_Time_Point.lt[["Start_Time"]] <- Sys.time() # %>% as.character()

#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("RColorBrewer")) install.packages("RColorBrewer"); library(RColorBrewer)
if(!require("cowplot")) install.packages("cowplot"); library(cowplot) # library(ggpubr)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("magrittr")) install.packages("magrittr"); library(magrittr)
if(!require("stats")) install.packages("stats"); library(stats)

# if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
# remotes::install_github("satijalab/seurat-data"); library(SeuratData) # force = TRUE
# # InstallData("pbmc3k")
# # pbmc <- LoadData("pbmc3k")

Rec_Time_Point.lt[["Load_Packages"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Load_Packages"]] <- Rec_Time_Point.lt[["Load_Packages"]]-Rec_Time_Point.lt[["Start_Time"]]

##### Parameter, inport and export setting* #####
### Set inport ###
# ## seurat_pbmc3k DeletB
# Path_Sample <- "D:/Dropbox/##_GitHub/###_VUMC/Metrics2CellAnnot/Input_Dataset_Seurat_pbmc3k/Seurat_pbmc3k_Seed123_Ref1350_Samp_SingleR_DeletB.RData"
# Path_Ref <- "D:/Dropbox/##_GitHub/###_VUMC/Metrics2CellAnnot/Input_Dataset_Seurat_pbmc3k/Seurat_pbmc3k_Seed123_Ref1350_Ref.RData"
# Set_Sample <- "pbmc3k_Sam1350_SingleR_DeletB"
# Set_Reference <- "pbmc3k_Ref1350_Ref"  # "HumanPrimaryCellAtlasData", "Self","GSE132044"

source("Set_LoadData.R") ##* Comment out if Run All

### Parameter ###
source("Set_Parameter.R")

## Export parameter information to txt file
writeLines(readLines("Set_Parameter.R"),
           con = paste0(Name_ExportFolder,"/",Name_Export,"_Parameter_Settings_Record.txt"))

##### Load data* #####
## Load sample
load(Path_Sample) #; rm(Name_Export_o,Name_ExportFolder_o)

if(!is.null(seuratObject@meta.data$`Annotation`)){
  seuratObject@meta.data$`Cell_Type` <- seuratObject@meta.data$`Annotation`
}

# if(!is.null(seuratObject@meta.data$`Cell_Type`)){
#   seuratObject@meta.data$`Annotation` <- seuratObject@meta.data$`Cell_Type`
# }

seuratObject_Sample <- seuratObject; rm(seuratObject)
if(Set_Sam_Delet_Unknown){
  seuratObject_Sample <- subset(seuratObject_Sample, subset = Actual_Cell_Type != "Unknown")
  seuratObject_Sample <- subset(seuratObject_Sample, subset = Cell_Type != "Unknown")
}

## SetIdent for seuratObject_Sample
seuratObject_Sample <- seuratObject_Sample %>% SetIdent(value = "Annotation")
DimPlot(seuratObject_Sample, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Annotation")

## Load reference
load(Path_Ref) #; rm(Name_Export_o,Name_ExportFolder_o)
seuratObject_Ref <- seuratObject; rm(seuratObject)
if(Set_Ref_Delet_Unknown){ seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != "Unknown") }

## Set_Ref_Delet_CTMetric
if(Set_Ref_Delet_CTMetric){
  if(!is.null(seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Set_Ref_Delet"]])){
    seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Set_Ref_Delet"]])
    # seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != "Naive CD4 T")
  }
}

Rec_Time_Point.lt[["Load_data"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Load_data"]] <- Rec_Time_Point.lt[["Load_data"]] - Rec_Time_Point.lt[["Load_Packages"]]

## Record Para to SeuratObj
source("Record_SetParameterToSeuratObj.R")

## SetIdent for seuratObject_Ref
seuratObject_Ref <- seuratObject_Ref %>% SetIdent(value = "Actual_Cell_Type")
DimPlot(seuratObject_Ref, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type")


##### Data Preprocessing #####
## Seurat object Prepocessing
source("FUN_Seurat_Prepocessing.R")
if(Set_Run_Process_Sample){ seuratObject_Sample <- Seurat_Prepocessing(seuratObject_Sample, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }
if(Set_Run_Process_Ref){ seuratObject_Ref <- Seurat_Prepocessing(seuratObject_Ref, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }

Rec_Time_Point.lt[["Preprocessing"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Preprocessing"]] <- Rec_Time_Point.lt[["Preprocessing"]] - Rec_Time_Point.lt[["Load_data"]]

##### Cell Type Annotation #####
source("Run_Cell_Type_Annotation.R")

Rec_Time_Point.lt[["Cell_Type_Annotation"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Cell_Type_Annotation"]] <- Rec_Time_Point.lt[["Cell_Type_Annotation"]] - Rec_Time_Point.lt[["Preprocessing"]]


##### VICTORS #####
source("FUN_VICTORS.R")

## singleR
VICTORS.lt <- VICTORS(seuratObject_Sample, seuratObject_Ref,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "label_singleR_NoReject")
seuratObject_Sample <- VICTORS.lt$Sample
seuratObject_Ref <- VICTORS.lt$Reference

## scmap
VICTORS.lt <- VICTORS(seuratObject_Sample, seuratObject_Ref,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "label_scmap_NoReject")
seuratObject_Sample <- VICTORS.lt$Sample
seuratObject_Ref <- VICTORS.lt$Reference


## SCINA
VICTORS.lt <- VICTORS(seuratObject_Sample, seuratObject_Ref,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "label_SCINA_NoReject")
seuratObject_Sample <- VICTORS.lt$Sample
seuratObject_Ref <- VICTORS.lt$Reference


## scPred
VICTORS.lt <- VICTORS(seuratObject_Sample, seuratObject_Ref,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "label_scPred_NoReject")
seuratObject_Sample <- VICTORS.lt$Sample
seuratObject_Ref <- VICTORS.lt$Reference


##### Evaluate VICTORS results #####
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

VICTORS_DiagPara <- function(seurat_obj, stat_var, Diagnosis,
                                 annotation_col = "Annotation",
                                 actual_col = "Actual Cell Type") {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(
      !!Diagnosis := case_when(
        !!sym(annotation_col) == !!sym(actual_col) & !!sym(stat_var) == "T" ~ "TP",
        !!sym(annotation_col) != !!sym(actual_col) & !!sym(stat_var) == "F" ~ "TN",
        !!sym(annotation_col) == !!sym(actual_col) & !!sym(stat_var) == "F" ~ "FN",
        !!sym(annotation_col) != !!sym(actual_col) & !!sym(stat_var) == "T" ~ "FP",
        TRUE ~ NA_character_
      )
    )
  return(seurat_obj)
}

score_types <- c("label_singleR_NoReject","label_scmap_NoReject","label_SCINA_NoReject","label_scPred_NoReject")
for (score_type in score_types) {
  seuratObject_Sample <- VICTORS_DiagPara(seuratObject_Sample, paste0("Diag_VICTORS_", score_type,"_StatROC"),
                                        paste0("DiagPara_VICTORS_",score_type,"_ROC"),
                                        annotation_col = score_type,actual_col = "Actual_Cell_Type")

}


##### Visualization #####
# ## Plot cell type count
# try({ source("Run_CellCount.R") })
source("Plot_Histograms_AnnoDiagnosis.R")


##### Export #####
## Export Metadata
write.table(data.frame(ID=rownames(seuratObject_Sample@meta.data), seuratObject_Sample@meta.data),
            file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataSamp.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
write.table(data.frame(ID=rownames(seuratObject_Ref@meta.data), seuratObject_Ref@meta.data),
            file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataRef.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

## Export Time Record
Rec_Time_Spend.df <- data.frame(
  Item = names(Rec_Time_Spend.lt),
  Time = unlist(lapply(Rec_Time_Spend.lt, as.numeric))
)
write.table(Rec_Time_Spend.df,
            file=paste0(Name_ExportFolder,"/",Name_Export,"_TimeSpendRecord.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

## Record version and sessionInfo
info_output <- c("##_R Version Information:", capture.output(version), "",
                 "##_Session Information:", capture.output(sessionInfo()))

try({ writeLines(info_output, paste0(Name_ExportFolder,"/",Name_Export,"_VerSesInfo.txt")) })


## Export RData

# Remove Plot Object
plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])

save.image(paste0(Name_ExportFolder,"/",Name_Export,".RData"))

# # Save small RData
# save(seuratObject_Sample, seuratObject_Ref, ROC_Summarize.lt,
#      file = paste0(Name_ExportFolder, "/", Name_Export, "_S.RData"))


