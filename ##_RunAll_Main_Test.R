### Parameter ###
source("##_RunAll_Set_Parameter.R")

##### Load data #####
## Load sample
load(Path_Sample) #; rm(Name_Export_o,Name_ExportFolder_o)

if(is.null(seuratObject@meta.data$`Actual_Cell_Type`)){
  seuratObject@meta.data$`Actual_Cell_Type` <- seuratObject@meta.data$`Cell_Type`
}
seuratObject <- UpdateSeuratObject(seuratObject)
seuratObject_Sample <- seuratObject; rm(seuratObject)
if(Set_Sam_Delet_Unknown){
  seuratObject_Sample <- subset(seuratObject_Sample, subset = Actual_Cell_Type != "Unknown")
  # seuratObject_Sample <- subset(seuratObject_Sample, subset = Cell_Type != "Unknown")
}

## SetIdent for seuratObject_Sample
seuratObject_Sample <- seuratObject_Sample %>% SetIdent(value = "Annotation")
DimPlot(seuratObject_Sample, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Annotation")

## Load reference
load(Path_Ref) #; rm(Name_Export_o,Name_ExportFolder_o)
seuratObject <- UpdateSeuratObject(seuratObject)
seuratObject_Ref <- seuratObject; rm(seuratObject)
if(Set_Ref_Delet_Unknown){ seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != "Unknown") }


if(Set_Ref_Delet_CTMetric){
  seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != Set_Ref_Delet)
}

## SetIdent for seuratObject_Ref
seuratObject_Ref <- seuratObject_Ref %>% SetIdent(value = "Actual_Cell_Type")
DimPlot(seuratObject_Ref, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type")


#### Save information to seuratObject ####
## Save information to metadata
seuratObject_Sample@meta.data$Sample_Platform <- seuratObject_Sample@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@meta.data$Sample_DataID <- seuratObject_Sample@misc[["BasicInfo"]][["DataID"]]
seuratObject_Sample@meta.data$Ref_Platform <- seuratObject_Ref@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@meta.data$Ref_DataID <- seuratObject_Ref@misc[["BasicInfo"]][["DataID"]]

seuratObject_Sample@meta.data$Mislabel_CellType <- Set_Ref_Delet

## Save information to misc
seuratObject_Sample@misc[["CTAnnot"]][["Ref_Platform"]] <- seuratObject_Ref@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@misc[["CTAnnot"]][["Ref_DataID"]] <- seuratObject_Ref@misc[["BasicInfo"]][["DataID"]]
seuratObject_Sample@misc[["CTAnnot"]][["Ref_PMID"]] <- seuratObject_Ref@misc[["BasicInfo"]][["PMID"]]


## Create sub folder
Name_PlatForm <- paste0("Qry_", seuratObject_Sample@misc[["BasicInfo"]][["Platform"]],
                        "_Ref_",seuratObject_Ref@misc[["BasicInfo"]][["Platform"]])
Name_ExportFolder <- paste0(Name_ExportFolder, "/", Name_FileID,"_",Set_Ref_Delet_Mislabel_Name,"_",Name_PlatForm)
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder

## Export parameter information to txt file
writeLines(readLines("##_RunAll_Set_Parameter.R"),
           con = paste0(Name_ExportFolder,"/",Name_Export,"_Parameter_Settings_Record.txt"))


##### Data Preprocessing #####
## Seurat object Prepocessing
source("FUN_Seurat_Prepocessing.R")
if(Set_Run_Process_Sample){ seuratObject_Sample <- Seurat_Prepocessing(seuratObject_Sample, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }
if(Set_Run_Process_Ref){ seuratObject_Ref <- Seurat_Prepocessing(seuratObject_Ref, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }


##### Main: Cell Type Annotation & DiagnosticMetrics & VICTOR #####
source("#_RUN_CellTypeAnnot_ConfStat_VICTOR.R")


##### scReClassify #####

##### Export #####



