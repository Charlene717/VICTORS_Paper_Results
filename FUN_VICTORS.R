# seuratObj_Samp <- seuratObject_Sample
# seuratObj_Samp$Annotation <- seuratObj_Samp$Cell_Type
# seuratObj_Ref <- seuratObject_Ref
# ActualCellTypeColumn = "Actual_Cell_Type"
# AnnotCellTypeColumn = "Cell_Type"

VICTORS <- function(seuratObj_Samp, seuratObj_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "Annotation",
                    Add_nROC = FALSE) {
  source("FUN_VICTORSPrep.R")

  #### Unified column name ####
  seuratObj_Ref[["Actual_Cell_Type"]] <- seuratObj_Ref[[ActualCellTypeColumn]]
  # seuratObj_Samp[["Annotation"]] <- seuratObj_Samp[[AnnotCellTypeColumn]]

  #### VICTORSPrep score ####
  if(is.null(seuratObj_Ref@misc[["VICTORS"]])) {
    seuratObj_Ref <- getFeatureSpace(seuratObj_Ref, "Actual_Cell_Type")  ## Get the feature space to train the classifiers
  }

  if(length(seuratObj_Ref@misc[["VICTORS"]]@train) == 0) {
    seuratObj_Ref <- trainModel(seuratObj_Ref,model = "glmnet") # Train the model
  }

  ## Ref
  if(!"VICTORS_max" %in% colnames(seuratObj_Ref@meta.data)) {
    seuratObj_Ref <- VICTORSPrep(seuratObj_Ref, seuratObj_Ref) #, threshold = Set_scPredict_Thr)
  }

  ## Sample
  if(!"VICTORS_max" %in% colnames(seuratObj_Samp@meta.data)) {
    seuratObj_Samp <- VICTORSPrep(seuratObj_Samp, seuratObj_Ref) #, threshold = Set_scPredict_Thr)
  }

  #### Modify column names ####
  modify_metadata_cols_VICTORS <- function(seuratObj) {
    colnames(seuratObj@meta.data) <- sapply(colnames(seuratObj@meta.data), function(name) {
      if(grepl("^VICTORS_", name) && !name %in% c("VICTORS_max", "VICTORS_prediction", "VICTORS_no_rejection")) {
        paste0(gsub("_", " ", gsub("_plus", "+", gsub("^VICTORS_", "", name))), " VICTORSScore")
      } else { name }
    })
    return(seuratObj) # Return the modified Seurat object
  }

  seuratObj_Ref = modify_metadata_cols_VICTORS(seuratObj_Ref)
  seuratObj_Samp = modify_metadata_cols_VICTORS(seuratObj_Samp)

  #### Ref: Optimization threshold ####
  source("FUN_ROC.R")
  ROC_Ref.lt <-   roc_analysis(seuratObj_Ref, Set_ACT = "Actual_Cell_Type",
                               Set_Anno = "Actual_Cell_Type", Set_Col = " VICTORSScore",
                               DefaultThr = 0.5)
  ROC_Ref.df <- generate_roc_data(ROC_Ref.lt, title_prefix = paste("Reference"))
  seuratObj_Ref@misc[["ROC"]][["ROC_Ref.lt"]] <- ROC_Ref.lt
  seuratObj_Ref@misc[["ROC"]][["ROC_Ref.df"]] <- ROC_Ref.df


  #### VICTORS Diagnosis ####
  score_type <- "VICTORS"
  metadata <- seuratObj_Samp@meta.data
  metadata_col <- paste0(AnnotCellTypeColumn, "_", score_type, "Score")

  ## Calculate the score and update metadata
  metadata[[metadata_col]] <- apply(metadata, 1, function(row) {
    score_column <- paste0(row[[AnnotCellTypeColumn]], " ", score_type, "Score")
    if (score_column %in% colnames(metadata)) row[[score_column]] else NA
  })

  metadata[[metadata_col]] <- metadata[[metadata_col]] %>% as.numeric()

  if(Add_nROC){
    ## State
    ScoreStat_Thr <- 0.5
    metadata[[paste0("Diag_", score_type,"_", AnnotCellTypeColumn, "_Stat")]] <- ifelse(metadata[[metadata_col]] < ScoreStat_Thr, "F", "T")
  }

  ## State ROC
  thresholds <- sapply(metadata[[AnnotCellTypeColumn]], function(cell_type) {
    if (cell_type %in% ROC_Ref.df$CellType) {
      ROC_Ref.df[ROC_Ref.df$CellType == cell_type, "OptimalThreshold"]
    } else {
      NA
    }
  })

  metadata[[paste0("Diag_", score_type,"_", AnnotCellTypeColumn,  "_StatROC")]] <- ifelse(metadata[[metadata_col]] < thresholds, "F", "T")

  seuratObj_Samp@meta.data <- metadata

  #### Export ####
  Output.lt <- list(Sample = seuratObj_Samp, Reference = seuratObj_Ref)
  return(Output.lt)
}


# seuratObject_Sample2 <- seuratObject_Sample
# seuratObject_Ref2 <- seuratObject_Ref
#
# VICTORS.lt <- VICTORS(seuratObject_Sample2,seuratObject_Ref2,
#                       ActualCellTypeColumn = "Actual_Cell_Type",
#                       AnnotCellTypeColumn = "Cell_Type")



