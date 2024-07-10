### Parameter ###
# GSE_Name = "GSE132044"

## Set Run ##
Set_Run_Process_Sample <- FALSE; Set_Run_Process_Ref <- FALSE
# Set_Run_Marker_Genes <- FALSE


## Set Prepocessing ##
Num_PCA <- 50
Set_nfeatures <- 2000 # nrow(seuratObject@assays[["RNA"]]), 2000

## Set Sam & Ref Delet ##
Set_Ref_Delet_CTMetric <- TRUE
Set_Ref_Delet_Unknown <- TRUE
Set_Sam_Delet_Unknown <- TRUE

## Settings for sampling the references ##
Set_Refer_Sam <- FALSE; Set_Seed <- 123; Set_RefNum <- 1350


## Cell Type Annotation ##
Set_RefAnnoCol <- "Actual_Cell_Type"


# Set_Ref_Delet <- "B cell" # Set_Ref_Delet <- "None" # Set_Ref_Delet <- c("Naive CD4 T") # Set_Ref_Delet <- c("B","NK")
# # Set_Ref_Delet_Mislabel_Name <- paste0("Mislabel",Set_Ref_Delet)
Set_Ref_Delet_Mislabel_Name <- paste0("Mislabel",Set_Ref_Delet)

# ## Settings for filtering marker genes ##
# Set_abs_logFC <- FALSE # TRUE, FALSE
# Set_fltMK_allCT = FALSE # TRUE, FALSE
# Set_MarkerSource <- "" # "Literature" ""
#
# transform_method <- "default" #  "default", "log10plus1",  "log1p"
# MK_logFC_Thr = 0.25; MK_pvalue_Thr = 0.01; MK_padj_Thr = 0.01; MK_Num = 20
# Set_testMK_Method = "wilcox" # "negbinom" # "wilcox"
#
# Set_Reference_aggregate <- "median" # median, mean, mode_stat_nonzero
# Set_Reference_FltMarker <- TRUE # TRUE, FALSE
# if(Set_Reference_FltMarker){
#   Set_MarkerbyFindVar <- "NaN"
#   Set_nfeatures_Simi <- "NaN"
# }else{
#   Set_MarkerbyFindVar <-  TRUE
#   Set_abs_logFC <- "NaN"
#   Set_nfeatures_Simi <- 200
# }


################################################################################

#### Set export #####
## Set export name
Name_CP <- "MSINB"
Name_Sup <- ""
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
                    # "_MKSrc",Set_MarkerSource, "_MKlogFC",gsub("\\.", "",MK_logFC_Thr), "_MKp",gsub("\\.", "",MK_pvalue_Thr),"_MKpadj",gsub("\\.", "",MK_padj_Thr), "_MKNum",MK_Num, "_TFM",transform_method,
                    # "_Fea",Set_nfeatures,"_nPC",Num_PCA,
                    "_Samp",Set_Sample,"_Ref", Set_Reference, Name_Ref, "_",Set_Ref_Delet_Mislabel_Name,
                    Name_Sup)

Fig_Title <- paste0(Name_FileID, " Samp:",Set_Sample," Ref:",Set_Reference, Name_Ref,
                    # " MKlogFC:",MK_logFC_Thr, " MKp:",MK_pvalue_Thr," MKBHp:",MK_padj_Thr, " MKNum:",MK_Num, # " TFM:",transform_method,
                    # "_Fea",Set_nfeatures,"_nPC",Num_PCA,
                    Name_Sup)


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


# Name_ExportFolder <- paste0(Name_ExportFolder, "/", Name_FileID,"_",Set_Ref_Delet_Mislabel_Name)
# if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
#
# # Name_PlatForm <- paste0("Qry_", seuratObject_Sample@misc[["BasicInfo"]][["Platform"]],
# #                         "_Ref_",seuratObject_Ref@misc[["BasicInfo"]][["Platform"]])
# # Name_ExportFolder <- paste0(Name_ExportFolder, "/", Name_FileID,"_",Set_Ref_Delet_Mislabel_Name,"_",Name_PlatForm)
# # if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
