if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)


# 定義函數以計算各細胞類型的Accuracy並更新Seurat物件
FUN_CTAnnot_Accuracy <- function(seuratObject, actualCellTypeField, labelField) {
  if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
  if (!inherits(seuratObject, "Seurat")) { stop("The input is not a valid Seurat object.") }

  # 從meta.data中提取實際細胞類型和註釋結果
  meta_data <- seuratObject@meta.data
  actual_types <- meta_data[[actualCellTypeField]]
  predicted_types <- meta_data[[labelField]]

  # 計算總體和各細胞類型的Accuracy
  accuracy_data <- meta_data %>%
    transmute(
      Actual = actual_types,
      Predicted = predicted_types,
      Is_Correct = Actual == Predicted
    )

  overall_accuracy <- mean(accuracy_data$Is_Correct)

  accuracy_per_type <- accuracy_data %>%
    group_by(Actual) %>%
    summarise(
      Per_Type_Accuracy = mean(Is_Correct),
      .groups = 'drop'
    )

  # 創建或更新misc槽位中的對應列表
  accuracy_list <- list(
    Per_Cell_Type = accuracy_per_type$Per_Type_Accuracy,
    Overall = overall_accuracy
  )

  labelField_M <- paste0(labelField, "_Accuracy")
  seuratObject@misc[[labelField_M]] <- accuracy_list

  return(seuratObject)
}

# ## Test Function
# seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoRejection')

################################################################################
FUN_Confusion_Matrix <- function(seuratObject, actualCellTypeField, originalLabelField, filteredLabelField) {
  if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
  if (!inherits(seuratObject, "Seurat")) { stop("The input is not a valid Seurat object.")}

  meta_data <- seuratObject@meta.data

  # 生成混淆矩陣狀態
  ConfusionStatus <- mapply(function(actual, original, filtered) {
    if (actual == original && filtered == original) {
      "TP"  # True Positive
    } else if (actual == original && filtered != original) {
      "FN"  # False Negative
    } else if (actual != original && (filtered == "Unassign" || filtered == actual)) {
      "TN"  # True Negative
    } else if (actual != original && filtered == original && original != "Unassign") {
      "FP"  # False Positive
    } else {
      "Other"
    }
  }, meta_data[[actualCellTypeField]], meta_data[[originalLabelField]], meta_data[[filteredLabelField]], SIMPLIFY = FALSE)

  # 更新meta.data
  meta_data$ConfusionStatus <- unlist(ConfusionStatus)
  seuratObject@meta.data[[paste0(filteredLabelField, "_ConfStat")]] <- meta_data$ConfusionStatus

  # 計算每類型細胞的性能指標
  per_cell_type_metrics <- meta_data %>%
    group_by(Actual_Cell_Type = meta_data[[actualCellTypeField]]) %>%
    summarise(
      TP = sum(ConfusionStatus == "TP"),
      FP = sum(ConfusionStatus == "FP"),
      FN = sum(ConfusionStatus == "FN"),
      TN = sum(ConfusionStatus == "TN")
    ) %>%
    mutate(
      Accuracy = (TP + TN) / (TP + FP + FN + TN),
      Sensitivity = TP / (TP + FN),
      Specificity = TN / (TN + FP),
      Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
      F1_Score = ifelse((Precision + Sensitivity) > 0, 2 * (Precision * Sensitivity) / (Precision + Sensitivity), 0)
    )

  # 計算總體性能指標
  overall_metrics <- summarise(per_cell_type_metrics,
                               Overall_TP = sum(TP),
                               Overall_FP = sum(FP),
                               Overall_FN = sum(FN),
                               Overall_TN = sum(TN)
  ) %>%
    mutate(
      Overall_Accuracy = (Overall_TP + Overall_TN) / (Overall_TP + Overall_FP + Overall_FN + Overall_TN),
      Overall_Sensitivity = Overall_TP / (Overall_TP + Overall_FN),
      Overall_Specificity = Overall_TN / (Overall_TN + Overall_FP),
      Overall_Precision = ifelse((Overall_TP + Overall_FP) > 0, Overall_TP / (Overall_TP + Overall_FP), 0),
      Overall_F1_Score = ifelse((Overall_Precision + Overall_Sensitivity) > 0, 2 * (Overall_Precision * Overall_Sensitivity) / (Overall_Precision + Overall_Sensitivity), 0)
    )

  # 儲存性能指標到 misc 槽
  seuratObject@misc[[paste0(filteredLabelField, "_Metrics")]] <- list(
    Per_Cell_Type = per_cell_type_metrics,
    Overall = overall_metrics
  )

  return(seuratObject)
}

# ## Test Function
# seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
#                                              "label_singleR_NoReject", "label_singleR")


################################################################################
FUN_Confusion_Matrix_DiagTools <- function(seurat_obj, stat_var, diagnosis_col,
                                           annotation_col = "Annotation",
                                           actual_col = "Actual_Cell_Type") {

  # 添加診斷列
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(
      !!diagnosis_col := case_when(
        !!sym(annotation_col) == !!sym(actual_col) & !!sym(stat_var) == "T" ~ "TP",
        !!sym(annotation_col) != !!sym(actual_col) & !!sym(stat_var) == "F" ~ "TN",
        !!sym(annotation_col) == !!sym(actual_col) & !!sym(stat_var) == "F" ~ "FN",
        !!sym(annotation_col) != !!sym(actual_col) & !!sym(stat_var) == "T" ~ "FP",
        TRUE ~ NA_character_
      )
    )

  # 使用基礎R函數來計算混淆矩陣
  confusion_matrix <- table(seurat_obj@meta.data[[diagnosis_col]])

  # 確認混淆矩陣包含所有需要的列
  all_labels <- c("TP", "TN", "FP", "FN")
  for (label in all_labels) {
    if (!(label %in% names(confusion_matrix))) {
      confusion_matrix[label] <- 0
    }
  }
  confusion_matrix <- as.data.frame(as.matrix(confusion_matrix))

  # 獲取性能指標
  TP <- confusion_matrix[confusion_matrix$diagnosis_vec == "TP", "Freq"]
  TN <- confusion_matrix[confusion_matrix$diagnosis_vec == "TN", "Freq"]
  FP <- confusion_matrix[confusion_matrix$diagnosis_vec == "FP", "Freq"]
  FN <- confusion_matrix[confusion_matrix$diagnosis_vec == "FN", "Freq"]

  overall_metrics <- list(
    Accuracy = (TP + TN) / (TP + TN + FP + FN),
    Precision = TP / (TP + FP),
    Recall = TP / (TP + FN),
    F1_Score = 2 * (TP / (TP + FP)) * (TP / (TP + FN)) / ((TP / (TP + FP)) + (TP / (TP + FN)))
  )

  # 計算每個細胞類型的性能指標
  per_cell_type_metrics <- seurat_obj@meta.data %>%
    group_by(!!sym(actual_col)) %>%
    summarise(
      Accuracy = sum(!!sym(diagnosis_col) == "TP" | !!sym(diagnosis_col) == "TN") / n(),
      Precision = sum(!!sym(diagnosis_col) == "TP") / sum(!!sym(diagnosis_col) == "TP" | !!sym(diagnosis_col) == "FP"),
      Recall = sum(!!sym(diagnosis_col) == "TP") / sum(!!sym(diagnosis_col) == "TP" | !!sym(diagnosis_col) == "FN"),
      F1_Score = 2 * Precision * Recall / (Precision + Recall)
    )

  # 儲存性能指標到 misc 槽
  seurat_obj@misc[[paste0(diagnosis_col, "_Metrics")]] <- list(
    Per_Cell_Type = per_cell_type_metrics,
    Overall = overall_metrics
  )

  return(seurat_obj)
}

# ## Test Function
# score_type <- "VICTOR"
# seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_",score_type,"_Stat"),
#                                                       paste0("DiagPara_", score_type),
#                                                       annotation_col = "Annotation")
