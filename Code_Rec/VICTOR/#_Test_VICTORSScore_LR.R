## Load Packages
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("glmnet")) install.packages("glmnet"); library(glmnet)

#### Rank Score ####
Set_RankScore = FALSE # Set_RankScore = TRUE
if(Set_RankScore == TRUE){
  # 從seuratObject選出符合條件的欄位
  score_columns_ref <- grep("VICTORSScore$", colnames(seuratObject_Ref@meta.data), value=TRUE)
  score_columns_sample <- grep("VICTORSScore$", colnames(seuratObject_Sample@meta.data), value=TRUE)

  # Function to rank values in a row and return the rankings
  rank_row <- function(row) { rank(row) }

  # Apply the ranking for seuratObject_Ref
  ranked_data_ref <- t(apply(seuratObject_Ref@meta.data[, score_columns_ref], 1, rank_row))
  colnames(seuratObject_Ref@meta.data) <- gsub("VICTORSScore$", "VICTORSScore.woRank", colnames(seuratObject_Ref@meta.data))
  seuratObject_Ref@meta.data <- cbind(seuratObject_Ref@meta.data, ranked_data_ref)

  # Apply the ranking for seuratObject_Sample
  ranked_data_sample <- t(apply(seuratObject_Sample@meta.data[, score_columns_sample], 1, rank_row))
  colnames(seuratObject_Sample@meta.data) <- gsub("VICTORSScore$", "VICTORSScore.woRank", colnames(seuratObject_Sample@meta.data))
  seuratObject_Sample@meta.data <- cbind(seuratObject_Sample@meta.data, ranked_data_sample)
}


#### Logistic Regression ####
# 從seuratObject_Ref@meta.data選出符合條件的欄位
score_columns_ref <- grep("VICTORSScore$", colnames(seuratObject_Ref@meta.data), value=TRUE)

# 針對每一種細胞類型建立羅吉斯回歸模型
cell_types <- unique(seuratObject_Ref@meta.data$"Actual_Cell_Type")
models <- list()

for (cell in cell_types) {
  y <- ifelse(seuratObject_Ref@meta.data$"Actual_Cell_Type" == cell, 1, 0)
  X <- as.matrix(seuratObject_Ref@meta.data[, score_columns_ref])

  model <- glmnet(X, y, family="binomial", alpha = 0.5) # alpha=1 for lasso

  models[[cell]] <- model
}

# summary(model)

# 使用來自所有的seuratObject_Ref模型預測seuratObject_Sample
predict_seurat_LR_All <- function(seuratObject, models, cell_types, score_columns) {
  results_df <- data.frame(matrix(nrow=nrow(seuratObject@meta.data), ncol=length(cell_types)))
  colnames(results_df) <- paste0(cell_types, " VICTORSScore")

  for (cell in cell_types) {
    X_data <- as.matrix(seuratObject@meta.data[, score_columns])

    prob <- predict(models[[cell]], newx=X_data, type="response", s=0.01)
    results_df[[paste0(cell, " VICTORSScore")]] <- prob[,1]

  }


  # Rename meta.data
  colnames(seuratObject@meta.data) <- gsub("VICTORSScore$", "VICTORSScore.1", colnames(seuratObject@meta.data))

  # Update metadata
  seuratObject@meta.data <- cbind(seuratObject@meta.data, results_df)

  return(seuratObject)
}

seuratObject_Sample = predict_seurat_LR_All(seuratObject_Sample, models, cell_types, score_columns_ref)
seuratObject_Ref = predict_seurat_LR_All(seuratObject_Ref, models, cell_types, score_columns_ref)

## Summarize LRScore
# 根據特定列名稱提取對應的score的函數
extract_value_from_column <- function(row, metadata, set_score_column = " score",
                                      cell_type_column = "Annotation", Set_PvStat = TRUE) {

  cell_type <- row[[cell_type_column]]
  score_column <- paste0(cell_type, set_score_column)
  score_value <- ifelse(score_column %in% colnames(metadata), row[[score_column]], NA)

  return(list(score_value ))
}

## Sample
metadata <- seuratObject_Sample@meta.data
values_list_LR <- lapply(1:nrow(metadata), function(i) extract_value_from_column(metadata[i,], metadata, set_score_column = " LRScore"))

# metadata$ScoreLR <- sapply(values_list_LR, `[[`, 1)
# metadata$LRStat <- ifelse(metadata$ScoreLR < 0.5, "F", "T")
# seuratObject_Sample@meta.data <- metadata

## Reference
metadata_Ref <- seuratObject_Ref@meta.data
values_list_LR_Ref <- lapply(1:nrow(metadata_Ref), function(i) extract_value_from_column(metadata_Ref[i,], metadata_Ref, set_score_column = " LRScore", cell_type_column = "Actual_Cell_Type"))

# metadata_Ref$ScoreLR <- sapply(values_list_LR_Ref, `[[`, 1)
# metadata_Ref$LRStat <- ifelse(metadata_Ref$ScoreLR < 0.5, "F", "T")
# seuratObject_Ref@meta.data <- metadata_Ref

# head(seuratObject_Sample@meta.data)


# #### ROC ####
# # try({
# #   if (Set_LGBM_ThrROC) {
# #     # source("##_Metrics2CellAnnot_MarkerScore_LR_ROC.R")
# #     source("FUN_ROC_Process_and_Store.R")
# #     if (!exists("ROC_Summarize.lt")) { ROC_Summarize.lt <- list() }
# #     ROC_Summarize.lt <- process_and_store_ROC(ROC_Summarize.lt,seuratObject_Sample,seuratObject_Ref,
# #                                               "LRScore", " LRScore", Name_ExportFolder, Name_Export)
# #   }
# # })
#
#
# if (Set_LR_ThrROC) {
#   metadata <- seuratObject_Sample@meta.data
#   # ref_data <- seuratObject_Sample@misc[["LR_ROC_Sum_df"]][seuratObject_Sample@misc[["LR_ROC_Sum_df"]]$Source == "Reference", ]
#   ref_data <- subset(ROC_Summarize.lt[["LRScore"]][["ROC_Set_Sum.df"]], Source == "Reference")
#
#   # 使用sapply遍歷metadata$Annotation，根據每個Cell Type獲取對應的OptimalThreshold
#   thresholds <- sapply(metadata$Annotation, function(cell_type) {
#     if (cell_type %in% ref_data$CellType) {
#       return(ref_data[ref_data$CellType == cell_type, "OptimalThreshold"])
#     } else {
#       return(NA)  # 若metadata中的CellType不存在於ref_data中，返回NA
#     }
#   })
#
#   # 根據新的thresholds設定metadata$LRStatROC的值
#   metadata$LRStatROC <- ifelse(metadata$ScoreLR < thresholds, "F", "T")
#
#   # Update metadata
#   seuratObject_Sample@meta.data <- metadata
# }
#
# #### Classification ####
# source("FUN_Classification.R")
#
# # Classification
# seuratObject_Sample <- AnnotMetric_4Class(seuratObject_Sample, "LRStat", "ClassificationLR")
# if (Set_LR_ThrROC) {seuratObject_Sample <- AnnotMetric_4Class(seuratObject_Sample, "LRStatROC", "ClassificationLRROC")}
#
#
# ##################################################################################
# ##### Visualization ####
# source("PlotFun_Histogram.R")
# source("PlotFun_CombineFigs.R")
# source("Set_pbmc3k_plot.R")
#
# #### Histograms ####
# metadata <- seuratObject_Sample@meta.data
#
# variables_LR <- list(
#   LR = list(Stat_var = "LRStat", Classification = "ClassificationLR", Note_Title = Note_TitleACT, color_vector = color_class),
#   ROCLR = list(Stat_var = "LRStatROC", Classification = "ClassificationLRROC", Note_Title = Note_TitleACT, color_vector = color_class)
# )
#
# # Generate plots
# plots_Stat_LR <- map(variables_LR, function(vars) {
#   if (vars$Stat_var == "LRStatROC" && !Set_LR_ThrROC) return(NULL)
#
#   map(c("Actual Cell Type","seurat clusters", "Annotation"), function(v) {
#     plot_histograms_for_vars(metadata, v, vars$Stat_var, vars$Note_Title)
#   })
# })
#
# plots_4Class_LR <- map(variables_LR, function(vars) {
#   if (vars$Stat_var == "LRStatROC" && !Set_LR_ThrROC) return(NULL)
#
#   map(c("Actual Cell Type","seurat clusters", "Annotation"), function(v) {
#     plot_histograms_for_vars(metadata, v, vars$Classification, vars$Note_Title, vars$color_vector)
#   })
# })
#
#
# #### UMAP ####
# if (exists("color_class")) {
#   Plot_UMAP_ClassLR <- DimPlot(seuratObject_Sample, group.by = "ClassificationLR", label = TRUE) +
#     scale_color_manual(values = color_class) # + theme(legend.position = "bottom")
# } else {
#   Plot_UMAP_ClassLR <- DimPlot(seuratObject_Sample, group.by = "ClassificationLR", label = TRUE)
# }
#
# if (Set_LR_ThrROC) {
#   if (exists("color_class")) {
#     Plot_UMAP_ClassLRROC <- DimPlot(seuratObject_Sample, group.by = "ClassificationLRROC", label = TRUE) +
#       scale_color_manual(values = color_class) # + theme(legend.position = "bottom")
#   } else {
#     Plot_UMAP_ClassLRROC <- DimPlot(seuratObject_Sample, group.by = "ClassificationLRROC", label = TRUE)
#   }
# }
#
# seuratObject_Sample_Temp <- seuratObject_Sample[,!is.na(seuratObject_Sample@meta.data$ScoreLR) ]
#
# Plot_UMAP_ScoreLR <- FeaturePlot(seuratObject_Sample_Temp, features = "ScoreLR", pt.size = 0.5, min.cutoff = 'q05', max.cutoff = 'q95')
#
# Plot_UMAP_LRStat <- DimPlot(seuratObject_Sample_Temp, group.by = "LRStat", label = TRUE)
#
# if (Set_LR_ThrROC) {Plot_UMAP_LRStatROC <- DimPlot(seuratObject_Sample_Temp, group.by = "LRStatROC", label = TRUE)}
#
# rm(seuratObject_Sample_Temp)
#
# #### Export pdf ####
# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_LRStat.pdf"),
#     width = 17, height = 17)
#
# ## LRROC
# if (Set_LR_ThrROC) {
#   combine_Figs(plots_4Class_LR$ROCLR)
#   combine_Figs(plots_Stat_LR$ROCLR)
# }
#
# if (Set_LR_ThrROC) {
#   try({
#     grid.arrange(
#       Plot_UMAP_ClassLRROC,
#       Plot_UMAP_ScoreLR,
#       Plot_UMAP_LRStatROC,
#       ncol=2)
#   })
# }
#
# ## LR
# combine_Figs(plots_4Class_LR$LR)
# combine_Figs(plots_Stat_LR$LR)
#
# try({
#   grid.arrange(
#     Plot_UMAP_ClassLR,
#     Plot_UMAP_ScoreLR,
#     Plot_UMAP_LRStat,
#     ncol=2)
# })
#
# dev.off()
#
#
# ## -[] Export Dataframe
# ## -[] Score Box plots
#
# #################################################################################
# #### Old Version ####
# # # 檢查ClassificationLR中的類別
# # table(seuratObject_Sample@meta.data$ClassificationLR)
# #
# # # 檢查color_class
# # if (exists("color_class")) {
# #   print(color_class)
# #
# #   # 如果上面的table的結果中的類別不完全符合color_class中的顏色，則手動調整color_class
# #   # 例如:
# #   # color_class <- c("T" = "red", "F" = "blue")
# #
# #   Plot_LR_UMAP_ClassPv <- DimPlot(seuratObject_Sample, group.by = "ClassificationLR", label = TRUE) +
# #     scale_color_manual(values = color_class)
# # } else {
# #   Plot_LR_UMAP_ClassPv <- DimPlot(seuratObject_Sample, group.by = "ClassificationLR", label = TRUE)
# # }
# #
# #
# # color_class <- list("Class 1" = "#b58b2a",
# #                     "Class 2" = "#e8bc56",
# #                     "Class 3" = "#368a5b",
# #                     "Class 4" = "#73bd94")
# # Plot_LR_UMAP_ClassPv <- DimPlot(seuratObject_Sample, group.by = "ClassificationLR", label = TRUE) +
# #   scale_color_manual(values = color_class)
