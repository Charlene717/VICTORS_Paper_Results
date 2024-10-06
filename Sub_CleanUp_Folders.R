## 設定路徑
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231107/Export_GSE132044_Test"
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231113_GSE132044_OriQC/#_AllMix_Export_GSE132044_Sum"
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240111_IntCT_scRNAseqPanc_Muraro/#_Export_20240111_Muraro_Baron_MislabelNone"
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240205_GSE132044_OriQC_2LR_Rank_padj001/#_AllMix_Export_20240205_GSE132044_OriQC_2LR_Rank_padj001"
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240205_GSE132044_OriQC_2LR_Rank_padj005/#_AllMix_Export_20240205_GSE132044_OriQC_2LR_Rank_padj005"
# path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231209_GSE132044_OriQC_DeBug/#_AllMix_Export_GSE132044_OriQC_DeBug_Sum"
# path <-"D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240205_GSE132044_OriQC_2LR_padj005/#_AllMix_Export_20240205_GSE132044_OriQC_2LR_padj005"
# path <-"D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results_HLCA_core/Export_HLCA_core_20240808"

path <-"D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell"


# 使用list.files()找到所有的.RData檔案
rdata_files <- list.files(path = path, pattern = "\\.RData$", recursive = TRUE, full.names = TRUE)

# 使用file.remove()刪除找到的.RData檔案
result <- file.remove(rdata_files)

# 檢查是否所有的檔案都被成功刪除
if(all(result)) {
  cat("All .RData files have been deleted successfully!\n")
} else {
  cat("Some .RData files were not deleted. Please check manually.\n")
}


# #### Clean More ####
# # 使用list.files()找到所有的.RData檔案
# rdata_files <- list.files(path = path, pattern = "\\.RData$", recursive = TRUE, full.names = TRUE)
# rdata_files2 <- list.files(path = path, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
#
#
# # 使用file.remove()刪除找到的.RData檔案
# result <- file.remove(rdata_files)
# result2 <- file.remove(rdata_files2)
#
# # 檢查是否所有的檔案都被成功刪除
# if(all(result) && all(result2)) {
#   cat("All .RData files have been deleted successfully!\n")
# } else {
#   cat("Some .RData files were not deleted. Please check manually.\n")
# }


# #################################################################################
# ## 你想要留下的是結尾為_S.RData的檔案，其他的.RData檔案要刪除。我們可以透過正則表達式來達成這個需求。
# # 設定路徑
# path <- "D:/Dropbox/##_GitHub/###_VUMC/Metrics2CellAnnot/Export_GSE132044_RefDeletB cell"
#
# # 使用list.files()找到所有的.RData檔案
# all_rdata_files <- list.files(path = path, pattern = "\\.RData$", recursive = TRUE, full.names = TRUE)
#
# # 篩選出不是以_S.RData結尾的檔案
# rdata_files_to_remove <- all_rdata_files[!grepl("_S\\.RData$", all_rdata_files)]
#
# # 使用file.remove()刪除找到的.RData檔案
# result <- file.remove(rdata_files_to_remove)
#
# # 檢查是否所有的檔案都被成功刪除
# if(all(result)) {
#   cat("All undesired .RData files have been deleted successfully!\n")
# } else {
#   cat("Some .RData files were not deleted. Please check manually.\n")
# }
