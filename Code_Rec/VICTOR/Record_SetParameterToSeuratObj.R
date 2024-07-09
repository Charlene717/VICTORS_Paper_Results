###############################################################################
# 創建一個名為Para_CellTypeAnnotMetrics的list來儲存所有參數
Para_CellTypeAnnotMetrics <- list()

# 將所有參數加入到Para_CellTypeAnnotMetrics中
Para_CellTypeAnnotMetrics$Set_Sample <- Set_Sample
Para_CellTypeAnnotMetrics$Set_Reference <- Set_Reference

if(!is.null(seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Set_Ref_Delet"]])){
  Para_CellTypeAnnotMetrics$CTAnnot_RefDelet <- seuratObject_Sample@misc[["Para_CellTypeAnnot"]][["Set_Ref_Delet"]]
}

try({ Para_CellTypeAnnotMetrics$GSE_Name <- GSE_Name })

Para_CellTypeAnnotMetrics$Name_ExportFolder <- Name_ExportFolder
Para_CellTypeAnnotMetrics$Name_Sup <- Name_Sup
Para_CellTypeAnnotMetrics$Name_Test <- Name_Test
Para_CellTypeAnnotMetrics$Name_CP <- Name_CP
Para_CellTypeAnnotMetrics$Name_time_wo_micro <- Name_time_wo_micro
Para_CellTypeAnnotMetrics$Name_time_day <- Name_time_day
try({ Para_CellTypeAnnotMetrics$Name_FileID <- Name_FileID })

Para_CellTypeAnnotMetrics$Name_Dataset <- Name_Dataset

Para_CellTypeAnnotMetrics$Set_PvStat_Thr <- Set_PvStat_Thr
Para_CellTypeAnnotMetrics$Set_BHPvStat_Thr <- Set_BHPvStat_Thr
Para_CellTypeAnnotMetrics$Set_PvStat_ThrROC <- Set_PvStat_ThrROC
Para_CellTypeAnnotMetrics$Set_LR_ThrROC <- Set_LR_ThrROC
Para_CellTypeAnnotMetrics$Set_LR_Fun <- Set_LR_Fun

try({ Para_CellTypeAnnotMetrics$Set_SVM_ThrROC <- Set_SVM_ThrROC })
try({ Para_CellTypeAnnotMetrics$Set_SVM_Kernel <- Set_SVM_Kernel })

try({ Para_CellTypeAnnotMetrics$Set_LGBM_ThrROC <- Set_LGBM_ThrROC })

try({ Para_CellTypeAnnotMetrics$Set_Score_Rank <- Set_Score_Rank })
try({ Para_CellTypeAnnotMetrics$Set_Ref_Delet_CTMetric <- Set_Ref_Delet_CTMetric })
try({ Para_CellTypeAnnotMetrics$Set_Sam_Delet_Unknown <- Set_Sam_Delet_Unknown })
try({ Para_CellTypeAnnotMetrics$Set_Ref_Delet_Unknown <- Set_Ref_Delet_Unknown })

Para_CellTypeAnnotMetrics$Set_abs_logFC <- Set_abs_logFC
Para_CellTypeAnnotMetrics$Set_fltMK_allCT <- Set_fltMK_allCT
Para_CellTypeAnnotMetrics$Set_MarkerSource <- Set_MarkerSource
Para_CellTypeAnnotMetrics$transform_method <- transform_method
Para_CellTypeAnnotMetrics$MK_logFC_Thr <- MK_logFC_Thr
Para_CellTypeAnnotMetrics$MK_pvalue_Thr <- MK_pvalue_Thr
Para_CellTypeAnnotMetrics$MK_padj_Thr <- MK_padj_Thr
Para_CellTypeAnnotMetrics$MK_Num <- MK_Num
Para_CellTypeAnnotMetrics$Set_testMK_Method <- Set_testMK_Method
Para_CellTypeAnnotMetrics$Set_Score <- Set_Score
Para_CellTypeAnnotMetrics$Set_Reference_aggregate <- Set_Reference_aggregate
Para_CellTypeAnnotMetrics$Set_Reference_FltMarker <- Set_Reference_FltMarker
Para_CellTypeAnnotMetrics$Set_MarkerbyFindVar <- Set_MarkerbyFindVar
Para_CellTypeAnnotMetrics$Set_nfeatures_Simi <- Set_nfeatures_Simi
Para_CellTypeAnnotMetrics$Set_Smooth <- Set_Smooth
Para_CellTypeAnnotMetrics$Set_Smooth_NumNBHD <- Set_Smooth_NumNBHD
try({Para_CellTypeAnnotMetrics$Set_Scale_Meta <- Set_Scale_Meta})
Para_CellTypeAnnotMetrics$Set_Run_Seurat_Integration <- Set_Run_Seurat_Integration
Para_CellTypeAnnotMetrics$Set_Integ_FeaNum <- Set_Integ_FeaNum
Para_CellTypeAnnotMetrics$Set_MdSC_CtrlNum <- Set_MdSC_CtrlNum
Para_CellTypeAnnotMetrics$Set_Refer_Sam <- Set_Refer_Sam
Para_CellTypeAnnotMetrics$Set_Seed <- Set_Seed
Para_CellTypeAnnotMetrics$Set_RefNum <- Set_RefNum
Para_CellTypeAnnotMetrics$Num_PCA <- Num_PCA
Para_CellTypeAnnotMetrics$Set_nfeatures <- Set_nfeatures

Para_CellTypeAnnotMetrics$Set_Run_Process_Sample <- Set_Run_Process_Sample
Para_CellTypeAnnotMetrics$Set_Run_Process_Ref <- Set_Run_Process_Ref

Para_CellTypeAnnotMetrics$Set_Run_scPred <- Set_Run_scPred
Para_CellTypeAnnotMetrics$Set_Run_SVGLRglm <- Set_Run_SVGLRglm
Para_CellTypeAnnotMetrics$Set_Run_SVGLRglmnet <- Set_Run_SVGLRglmnet
Para_CellTypeAnnotMetrics$Set_Run_EnsembleE2 <- Set_Run_EnsembleE2


# 將Para_CellTypeAnnotMetrics這個list儲存到seuratObject_Sample的misc之中
seuratObject_Sample@misc$Para_CellTypeAnnotMetrics <- Para_CellTypeAnnotMetrics
seuratObject_Ref@misc$Para_CellTypeAnnotMetrics <- Para_CellTypeAnnotMetrics

Para_CellTypeAnnotMetrics.df <- Para_CellTypeAnnotMetrics %>% as.data.frame() %>% t() %>% as.data.frame()
Para_CellTypeAnnotMetrics.df <- data.frame(Item = row.names(Para_CellTypeAnnotMetrics.df),
                                           Content = Para_CellTypeAnnotMetrics.df[,1])

write.table(Para_CellTypeAnnotMetrics.df,
            file=paste0(Name_ExportFolder,"/",Name_Export,"_Parameter.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

rm(Para_CellTypeAnnotMetrics,Para_CellTypeAnnotMetrics.df)
