Set_RefAnnoCol <- "Actual_Cell_Type"

# Check if Set_Ref_Delet is NULL
if (!is.null(Set_Ref_Delet)) {
  # Corresponding to the cell type to be deleted
  cells_to_delete <- rownames(seuratObject_Ref@meta.data)[which(seuratObject_Ref@meta.data[["Actual_Cell_Type"]] %in% Set_Ref_Delet)]

  # Delete the corresponding cell from the seuratObject_Ref object
  seuratObject_RefM <- subset(seuratObject_Ref, cells = setdiff(Cells(seuratObject_Ref), cells_to_delete))
}


##### Run Cell Type Annotation ####
source("#_FUN_CellTypeAnnot.R")
source("Plot_CellAnnot_UMAP_Box.R")
source("FUN_Metrics_CellTypeAnnot.R")
source("Set_plot_color.R")

if(Set_Run_SingleR){
  ## Run SingleR
  # source("##_Run_singleR.R")
  SingleR.lt <- Run_singleR(seuratObject_Sample, seuratObject_RefM)
  seuratObject_Sample@meta.data$label_singleR_NoReject <- SingleR.lt$labels

  ## Annotation diagnostics
  seuratObject_Sample@meta.data[[paste0("label_singleR")]] <- SingleR.lt$pruned.labels
  seuratObject_Sample@meta.data$label_singleR <- ifelse(is.na(seuratObject_Sample@meta.data$label_singleR), "Unassign", seuratObject_Sample@meta.data$label_singleR)

  ## Classification and Metrics
  seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')
  seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                               "label_singleR_NoReject", "label_singleR")

  ## UMAP_Box
  try({
    plot_seurat_data(seuratObject_Sample, seuratObject_RefM,
                     label_column_name = "label_singleR",label_column_name2 = "label_singleR_NoReject",
                     DiagPara_column_name = "label_singleR_DiagPara",
                     export_folder = Name_ExportFolder, export_name = paste0(Name_Export,"_singleR"))
  })

}


if(Set_Run_scmap){
  seuratObject_Sample <- Run_scmap(seuratObject_Sample, seuratObject_RefM)

  ## Classification and Metrics
  seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scmap_NoReject')
  seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                               "label_scmap_NoReject", "label_scmap")

  ## UMAP_Box
  try({
    plot_seurat_data(seuratObject_Sample, seuratObject_RefM,
                     label_column_name = "label_scmap", label_column_name2 = "label_scmap_NoReject",
                     DiagPara_column_name = "label_scmap_DiagPara",
                     export_folder = Name_ExportFolder, export_name =  paste0(Name_Export,"_scmap"))
  })

}

if(Set_Run_SCINA){
  seuratObject_Sample <- Run_SCINA(seuratObject_Sample, seuratObject_RefM,
                                   ExportFolder = Name_ExportFolder, Export = Name_Export)

  ## Classification and Metrics
  seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_SCINA_NoReject')
  seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                               "label_SCINA_NoReject", "label_SCINA")
  ## UMAP_Box
  try({
    plot_seurat_data(seuratObject_Sample, seuratObject_RefM,
                     label_column_name = "label_SCINA", label_column_name2 = "label_SCINA_NoReject",
                     DiagPara_column_name = "label_SCINA_DiagPara",
                     export_folder = Name_ExportFolder, export_name =  paste0(Name_Export,"_SCINA"))

  })

}

if(Set_Run_scPred){
  if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
  seuratObject_RefM <- getFeatureSpace(seuratObject_RefM, "Actual_Cell_Type")   ## Get the feature space to train the classifiers
  seuratObject_RefM <- trainModel(seuratObject_RefM) ## Train the model

  seuratObject_Sample <- scPredict(seuratObject_Sample, seuratObject_RefM) #, threshold = Set_scPredict_Thr)
  # seuratObject_RefM <- scPredict(seuratObject_RefM, seuratObject_RefM)

  source("FUN_modify_colnames_scPred.R")
  seuratObject_Sample <- modify_colnames_scPred(seuratObject_Sample,"scPred")
  # seuratObject_RefM <- modify_colnames_scPred(seuratObject_RefM,"scPred")

  ## Classification and Metrics
  seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scPred_NoReject')
  seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                               "label_scPred_NoReject", "label_scPred")
  ## UMAP_Box
  try({plot_seurat_data(seuratObject_Sample, seuratObject_RefM,
                        label_column_name = "label_scPred", label_column_name2 = "label_scPred_NoReject",
                        DiagPara_column_name = "label_scPred_DiagPara",
                        export_folder = Name_ExportFolder, export_name =  paste0(Name_Export,"_scPred"))
  })
}


rm(seuratObject_RefM)


#### Save information to metadata ####
seuratObject_Sample@meta.data$Sample_Platform <- seuratObject_Sample@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@meta.data$Sample_DataID <- seuratObject_Sample@misc[["BasicInfo"]][["DataID"]]
seuratObject_Sample@meta.data$Ref_Platform <- seuratObject_Ref@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@meta.data$Ref_DataID <- seuratObject_Ref@misc[["BasicInfo"]][["DataID"]]

seuratObject_Sample@meta.data$Mislabel_CellType <- Set_Ref_Delet

seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Ref_Platform"]] <- seuratObject_Ref@misc[["BasicInfo"]][["Platform"]]
seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Ref_DataID"]] <- seuratObject_Ref@misc[["BasicInfo"]][["DataID"]]
seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Ref_PMID"]] <- seuratObject_Ref@misc[["BasicInfo"]][["PMID"]]
