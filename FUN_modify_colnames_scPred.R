modify_colnames_scPred <- function(seuratObject, scPredType) {
  # 获取当前的列名
  current_colnames <- colnames(seuratObject@meta.data)

  # 创建一个新的列名列表
  new_colnames <- character(length(current_colnames))

  # 遍历当前列名
  for (i in 1:length(current_colnames)) {
    colname <- current_colnames[i]
    if (startsWith(colname, "scpred_") || startsWith(colname, "scPred_")) {
      if (colname %in% c("scpred_max", "scPred_max",
                         "scpred_label","scPred_label",
                         "scpred_no_rejection", "scPred_no_rejection",
                         "scpred_prediction", "scPred_prediction")) {
        # 保留不修改的列名
        new_colnames[i] <- colname
      } else {
        # 根据规则修改列名
        new_colnames[i] <- gsub("scpred_", "", colname)
        new_colnames[i] <- gsub("scPred_", "", new_colnames[i])
        new_colnames[i] <- gsub("\\.", " ", new_colnames[i])
        new_colnames[i] <- gsub("_plus", "+", new_colnames[i])

        new_colnames[i] <- paste0(new_colnames[i], paste0(" ",scPredType,"Score"))

      }
    } else {
      # 非以"scpred_"或"scPred_"开头的列名保持不变
      new_colnames[i] <- colname
    }
  }

  # 更新Seurat对象的列名
  colnames(seuratObject@meta.data) <- new_colnames

  # Dynamically generating the column name
  label_scPred_col <- paste0("label_", scPredType)

  colnames(seuratObject@meta.data) <- gsub("sc[pP]red_prediction", label_scPred_col, colnames(seuratObject@meta.data))
  colnames(seuratObject@meta.data) <- gsub("sc[pP]red_no_rejection", paste0(label_scPred_col,"_NoReject"), colnames(seuratObject@meta.data))
  colnames(seuratObject@meta.data) <- gsub("sc[pP]red_max",  paste0(label_scPred_col,"_Score"), colnames(seuratObject@meta.data))

  # Dynamically modifying the column in the metadata of the seuratObject
  seuratObject@meta.data <- seuratObject@meta.data %>%
    mutate(!!label_scPred_col := if_else(!!sym(label_scPred_col) == "unassigned", "Unassign", !!sym(label_scPred_col)))

  return(seuratObject)
}

# ## Test Function
# seuratObject_Sample <- modify_colnames_scPred(seuratObject_Sample, "scPredScore")
# seuratObject_Ref <- modify_colnames_scPred(seuratObject_Ref, "scPredScore")
#
# seuratObject_Sample <- modify_colnames_scPred(seuratObject_Sample, "scPredLRglmScore")
# seuratObject_Sample <- modify_colnames_scPred(seuratObject_Sample, "scPredLRglmnetScore")


