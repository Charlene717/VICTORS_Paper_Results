### Parameter ###
# GSE_Name = "GSE132044"
GSE_Name = "GSE136831"

## Set Annotation Check ##
Set_Annot <- c("label_singleR_NoReject", "label_scmap_NoReject", "label_SCINA_NoReject",
               "label_scPred_NoReject_Annot") #, "label_SVGLRglmnet_NoReject")
# Set_Annot <- "label_scPred_NoReject" # "label_scPred" # Set_Annot <- "Annotation"

## Set score type ##
Set_Score <- c("scPred", "SVGLRglmnet", "EnsembleE2")
# Set_Score <- "SVGLRglmnet" # "EnsembleE2Score"
#                                 # "scPredScore" # "scPredLRglmScore" # "scPredLRglmnetScore"
#                                 ## "ModuleScore","AUCellScore","UCellScore, "JASMINE"
#                                 ## "MKExpScore", "cosine", "spearman", "chatterjee", "Enrichment"



## Set Run ##
Set_Run_Process_Sample <- FALSE; Set_Run_Process_Ref <- FALSE
Set_Run_Marker_Genes <- FALSE

Set_Run_scPred <- TRUE
Set_Run_SVGLRglm <- FALSE
Set_Run_SVGLRglmnet <- TRUE
Set_Run_EnsembleE2 <- TRUE

## Set Prepocessing ##
Num_PCA <- 50
Set_nfeatures <- 2000 # nrow(seuratObject@assays[["RNA"]]), 2000

## Set Sam & Ref Delet ##
Set_Ref_Delet_CTMetric <- TRUE
Set_Ref_Delet_Unknown <- TRUE
Set_Sam_Delet_Unknown <- TRUE

## Settings for sampling the references ##
Set_Refer_Sam <- FALSE; Set_Seed <- 123; Set_RefNum <- 1350

## Settings for Integration ##
Set_Run_Seurat_Integration <- FALSE # TRUE, FALSE
Set_Integ_FeaNum <- 2000
if(Set_Run_Seurat_Integration){
  Set_MdSC_CtrlNum <- 80
}else{Set_MdSC_CtrlNum <- 100}

## Cell Type Annotation ##
Set_RefAnnoCol <- "Actual_Cell_Type"
Set_Run_SingleR <- TRUE
Set_Run_scPred <- TRUE
Set_Run_scmap <- TRUE
Set_Run_SCINA <- TRUE

# Set_Ref_Delet <- "B cell" # Set_Ref_Delet <- "None" # Set_Ref_Delet <- c("Naive CD4 T") # Set_Ref_Delet <- c("B","NK")
# # Set_Ref_Delet_Mislabel_Name <- paste0("Mislabel",Set_Ref_Delet)
Set_Ref_Delet_Mislabel_Name <- paste0("Mislabel",Set_Ref_Delet)

## Settings for filtering marker genes ##
Set_abs_logFC <- FALSE # TRUE, FALSE
Set_fltMK_allCT = FALSE # TRUE, FALSE
Set_MarkerSource <- "" # "Literature" ""

transform_method <- "default" #  "default", "log10plus1",  "log1p"
MK_logFC_Thr = 0.25; MK_pvalue_Thr = 0.01; MK_padj_Thr = 0.01; MK_Num = 20
Set_testMK_Method = "wilcox" # "negbinom" # "wilcox"

Set_Reference_aggregate <- "median" # median, mean, mode_stat_nonzero
Set_Reference_FltMarker <- TRUE # TRUE, FALSE
if(Set_Reference_FltMarker){
  Set_MarkerbyFindVar <- "NaN"
  Set_nfeatures_Simi <- "NaN"
}else{
  Set_MarkerbyFindVar <-  TRUE
  Set_abs_logFC <- "NaN"
  Set_nfeatures_Simi <- 200
}


## Settings for smoothing scores ##
Set_Smooth <- TRUE
if(Set_Smooth){
  Set_Smooth_NumNBHD = 5
}else{
  Set_Smooth_NumNBHD = ""
}

## Settings for Scaling ##
Set_Scale_Meta <- NULL # NULL # "zscore" "iqr_median"


## Set score rank ##
Set_Score_Rank <- FALSE

## Set LR ##
# Set_LR_Rank <- FALSE
Set_LR_Fun <- "glm"

## Set SVM ##
Set_SVM_Kernel = "linear" # "linear", "radial"
Set_SVM_ThrROC = TRUE

## Set LGBM ##
# Set_LGBM_Rank <- FALSE
Set_LGBM_ThrROC <- TRUE

## Threshold for Discrimination Results ##
Set_PvStat_Thr <- 0.01; Set_BHPvStat_Thr <- 0.01
Set_PvStat_ThrROC <- TRUE
Set_LR_ThrROC <- TRUE
Set_ScoreStat_Thr <- 0.55







## Set export ###
## Set export name
Name_Sup <- ""
Name_Test <- "" # Name_Test <- "_DeletB" # Name_Test <- "_RmRefInSam"
Name_CP <- "MSINB"
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

#Bug set.seed(as.integer(Sys.time()))
# Save Set seed
previous_seed <- .Random.seed
set.seed(NULL)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 6), collapse = ""))
# Recover Set seed
.Random.seed <- previous_seed

Name_time_day <- gsub("-", "_",Sys.Date())
Name_Dataset <- GSE_Name

if(Set_Refer_Sam){
  Name_Ref <- paste0("_RefSam", Set_Refer_Sam, "_RefSeed",Set_Seed, "_RefNum",Set_RefNum)
}else{
  Name_Ref <- ""
}

Name_Note <- paste0(Name_Sup,
                    # Set_Score,"_Smooth",Set_Smooth, Set_Smooth_NumNBHD,
                    # "_PvThr", gsub("\\.", "",Set_PvStat_Thr),
                    # "_Integ",Set_Run_Seurat_Integration, "_FltMKbyAllCT", Set_fltMK_allCT,  "_abslogFC",Set_abs_logFC,
                    # "_MKSrc",Set_MarkerSource, "_MKlogFC",gsub("\\.", "",MK_logFC_Thr), "_MKp",gsub("\\.", "",MK_pvalue_Thr),"_MKpadj",gsub("\\.", "",MK_padj_Thr), "_MKNum",MK_Num, "_TFM",transform_method,
                    # "_Fea",Set_nfeatures,"_nPC",Num_PCA,
                    "_Samp",Set_Sample,"_Ref", Set_Reference, Name_Ref, Name_Test,"_",Set_Ref_Delet_Mislabel_Name)

Fig_Title <- paste0(Name_FileID, " Samp:",Set_Sample," Ref:",Set_Reference, Name_Ref,
                    " Smooth:",Set_Smooth," (Num", Set_Smooth_NumNBHD,")",
                    # "_MKSrc",Set_MarkerSource,
                    # " PvStatThr:",Set_PvStat_Thr,
                    # " Integ:",Set_Run_Seurat_Integration, " FltMKbyAllCT:", Set_fltMK_allCT, " abslogFC:",Set_abs_logFC,
                    " MKlogFC:",MK_logFC_Thr, " MKp:",MK_pvalue_Thr," MKBHp:",MK_padj_Thr, " MKNum:",MK_Num, # " TFM:",transform_method,
                    # "_Fea",Set_nfeatures,"_nPC",Num_PCA,
                    " ",Name_Test)


# Name_Export <- paste0(Name_CP,"_",Name_FileID,
#                       # Name_time_day,"_",
#                       "_",Name_Note) # as.numeric(Sys.time(), units = "secs") %>% trunc())

Name_Export <- paste0(Name_FileID)

## Set ExportFolder ##
# Name_ExportFolder <- "Export"
Note_ExportFolder <- "" # Note_ExportFolder <- paste0("_",Set_Fin_Anno)
Name_ExportFolder <- paste0("Export","_",Name_Dataset, Note_ExportFolder,"_",Set_Ref_Delet_Mislabel_Name)


# Name_ExportFolder <- paste0("Export_",Name_Dataset,"_Ref",seuratObject@misc[["CellTypeAnnot_Para"]][["Set_Ref_Delet_Name"]])
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
if (length(Set_Score) == 1) {
  Name_ExportFolder <- paste0(Name_ExportFolder, "/", Name_FileID, "_", Set_Score)
} else {
  Name_ExportFolder <- paste0(Name_ExportFolder, "/", Name_FileID, "_Multi")
}
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
