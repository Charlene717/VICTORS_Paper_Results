if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)


# 定義函數以計算各細胞類型的Accuracy並更新Seurat物件
FUN_Accuracy <- function(seuratObject, actualCellTypeField, labelField) {
  # 確認seuratObject是Seurat物件
  if (!inherits(seuratObject, "Seurat")) {
    stop("The input is not a valid Seurat object.")
  }

  # 從meta.data中提取實際細胞類型和註釋結果
  actual_types <- seuratObject@meta.data[[actualCellTypeField]]
  predicted_types <- seuratObject@meta.data[[labelField]]

  # 計算總體Accuracy
  overall_accuracy <- sum(actual_types == predicted_types) / length(actual_types)

  # 計算每種細胞類型的Accuracy
  accuracy_per_type <- actual_types %>%
    as.factor() %>%
    levels() %>%
    setNames(nm = .) %>%
    lapply(function(cell_type) {
      actual <- actual_types[actual_types == cell_type]
      predicted <- predicted_types[actual_types == cell_type]
      sum(actual == predicted) / length(actual)
    }) %>%
    unlist() %>%
    as.data.frame()
  colnames(accuracy_per_type) <- "Accuracy"

  # 創建或更新misc槽位中的對應列表
  if (!labelField %in% names(seuratObject@misc)) {
    seuratObject@misc[[labelField]] <- list()
  }

  # 將總體Accuracy和每種細胞類型的Accuracy存入misc槽
  seuratObject@misc[[labelField]][["Accuracy_Per_Cell_Type"]] <- accuracy_per_type
  seuratObject@misc[[labelField]][["Overall_Accuracy"]] <- overall_accuracy

  # 返回更新後的Seurat物件
  return(seuratObject)
}

# ## Test Function
# seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoRejection')

################################################################################
FUN_DiagnosticMetrics <- function(seuratObject, actualCellTypeField,
                                  originalLabelField, filteredLabelField) {

  # 确认seuratObject是Seurat对象
  if (!inherits(seuratObject, "Seurat")) {
    stop("The input is not a valid Seurat object.")
  }

  # 获取实际和预测的细胞类型
  actual_types <- seuratObject@meta.data[[actualCellTypeField]]
  original_labels <- seuratObject@meta.data[[originalLabelField]]
  filtered_labels <- seuratObject@meta.data[[filteredLabelField]]

  # # 计算TP, FP, FN, TN
  # DiagPara <- ifelse(actual_types == original_labels & actual_types == filtered_labels, "TP",
  #                             ifelse(actual_types == original_labels & actual_types != filtered_labels, "FN",
  #                                    ifelse(actual_types != original_labels & filtered_labels == "Unassigned", "TN", "FP")))

  ## 计算诊断度量（TP, FP, FN, TN）
  # DiagPara <- mapply(function(actual, original, filtered) {
  #   if (actual == original && actual == filtered) {
  #     "TP"  # True Positive
  #   } else if (actual == original && filtered != actual) {
  #     "FN"  # False Negative
  #   } else if (actual != original && filtered == "Unassign") {
  #     "TN"  # True Negative
  #   } else if (actual != original && filtered != "Unassign") {
  #     "FP"  # False Positive
  #   }
  # }, actual_types, original_labels, filtered_labels)

  DiagPara <- mapply(function(actual, original, filtered) {
    if (actual == original && filtered == original) {
      "TP"  # True Positive
    } else if (actual == original && filtered != original) {
      "FN"  # False Negative
    } else if (actual != original && (filtered == "Unassign" || filtered == actual)) {
      "TN"  # True Negative
    } else if (actual != original && filtered == original) {
      "FP"  # False Positive
    } else {"Other"}
  }, actual_types, original_labels, filtered_labels)

  # 新增DiagPara列到metadata
  seuratObject@meta.data$DiagPara <- DiagPara
  seuratObject@meta.data[[paste0(filteredLabelField,"_DiagPara")]] <- DiagPara

  # 按照类型分组，计算各个指标
  group_by_type <- seuratObject@meta.data %>%
    group_by(Actual_Cell_Type = actual_types) %>%
    summarise(
      TP = sum(DiagPara == "TP"),
      FP = sum(DiagPara == "FP"),
      FN = sum(DiagPara == "FN"),
      TN = sum(DiagPara == "TN"),
      .groups = 'drop'
    )

  # 计算性能指标
  group_by_type <- group_by_type %>%
    mutate(
      Accuracy = (TP + TN) / (TP + FP + FN + TN),
      Sensitivity = TP / (TP + FN),
      Specificity = TN / (TN + FP),
      Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
      F1_Score = ifelse((Precision + Sensitivity) > 0, 2 * (Precision * Sensitivity) / (Precision + Sensitivity), 0)
    )

  # 存储各个指标到misc槽
  seuratObject@misc[[paste0(filteredLabelField,"_Diag")]] <- list()
  seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Per_Cell_Type"]] <- group_by_type

  # 计算总体指标
  overall_metrics <- group_by_type %>%
    summarise(
      TP = sum(TP),
      FP = sum(FP),
      FN = sum(FN),
      TN = sum(TN)
    ) %>%
    mutate(
      Overall_Accuracy = (TP + TN) / (TP + FP + FN + TN),
      Overall_Sensitivity = TP / (TP + FN),
      Overall_Specificity = TN / (TN + FP),
      Overall_Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
      Overall_F1_Score = ifelse((Overall_Precision + Overall_Sensitivity) > 0, 2 * (Overall_Precision * Overall_Sensitivity) / (Overall_Precision + Overall_Sensitivity), 0)
    )

  # 将总体指标也存入misc槽
  seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Overall"]] <- overall_metrics


  # 定义一个函数来将dataframe中的NA和NaN替换为NaN
  replace_na_with_NaN_df <- function(df) {
    # 逐个元素检查是否为NA或NaN，并替换
    df[] <- lapply(df, function(x) ifelse(is.na(x) | is.nan(x), NaN, x))
    return(df)
  }

  seuratObject@meta.data$DiagPara <- NULL

  seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Per_Cell_Type"]] <- replace_na_with_NaN_df(seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Per_Cell_Type"]])
  seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Overall"]] <- replace_na_with_NaN_df(seuratObject@misc[[paste0(filteredLabelField,"_Diag")]][["Overall"]])


  return(seuratObject)
}

# ## Test Function
# seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
#                                              "label_singleR_NoReject", "label_singleR")


