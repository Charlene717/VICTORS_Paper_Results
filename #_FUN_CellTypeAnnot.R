##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("devtools")) install.packages("devtools"); library(devtools)
# if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}; library(SeuratData)
if (!requireNamespace("SeuratData", quietly = TRUE)) {devtools::install_github('satijalab/seurat-data')}; library(SeuratData)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
if(!require("scater")) BiocManager::install("scater") ; library(scater)
if(!require("celldex")) BiocManager::install("celldex"); library(celldex)
if(!require("scuttle")) BiocManager::install("scuttle"); library(scuttle)


# devtools::install_github('dviraran/SingleR'); library(SingleR)
# # this might take long, though mostly because of the installation of Seurat.


# # Ensure the default assay is set to RNA
# DefaultAssay(Query_Seurat) <- "RNA"
# DefaultAssay(Reference_Seurat) <- "RNA"



#### singleR ####
Run_singleR <- function(Query_Seurat, Reference_Seurat,
                        Set_RefAnnoCol = "Actual_Cell_Type", seurat_version = "V5", ...) {
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("SingleR", quietly = TRUE)) BiocManager::install("SingleR"); library(SingleR)
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("scater", quietly = TRUE)) BiocManager::install("scater"); library(scater)


  # Preprocess reference dataset
  ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  ref_sce$label <- ref_sce[[Set_RefAnnoCol]]
  ref_sce <- ref_sce[, !is.na(ref_sce$label)]

  # Log-normalize reference dataset if not already done
  if(seurat_version == "V5"){
    if (is.null(Reference_Seurat@assays[["RNA"]]@data)) { ref_sce <- logNormCounts(ref_sce) }
  }else if(seurat_version == "V5M"){
    if (is.null(Reference_Seurat@assays[["RNA"]]@layers[["data"]])) { ref_sce <- logNormCounts(ref_sce) }
  }else{
    if (is.null(Reference_Seurat@assays[["RNA"]]@data)) { ref_sce <- logNormCounts(ref_sce) }
  }


  # Preprocess query dataset
  query_sce <- as.SingleCellExperiment(Query_Seurat)
  try(query_sce <- query_sce[, colSums(counts(query_sce)) > 0])

  # Log-normalize query dataset if not already done
  if(seurat_version == "V5"){
    # if (is.null(Query_Seurat@assays[["RNA"]]@counts) || is.null(Query_Seurat@assays[["RNA"]]@data)) { query_sce <- logNormCounts(query_sce) }
    if (is.null(Query_Seurat@assays[["RNA"]]@data)) { query_sce <- scater::logNormCounts(query_sce) }
  }else if(seurat_version == "V5M"){
    # try( if (is.null(Query_Seurat@assays[["RNA"]]@layers[["counts"]]) || is.null(Query_Seurat@assays[["RNA"]]@layers[["data"]])) { query_sce <- logNormCounts(query_sce) } )
    if (is.null(Query_Seurat@assays[["RNA"]]@layers[["data"]])) { query_sce <- scater::logNormCounts(query_sce) }
  }else{
    # if (is.null(Query_Seurat@assays[["RNA"]]@counts) || is.null(Query_Seurat@assays[["RNA"]]@data)) { query_sce <- logNormCounts(query_sce) }
    if (is.null(Query_Seurat@assays[["RNA"]]@data)) { query_sce <- scater::logNormCounts(query_sce) }

  }

  # Run SingleR
  SingleR.lt <- SingleR(test = query_sce, ref = ref_sce, assay.type.test = 1,
                        labels = ref_sce$label, de.method = "classic")

  Query_Seurat@meta.data$label_singleR_NoReject <- SingleR.lt$labels

  ## Annotation diagnostics
  Query_Seurat@meta.data[[paste0("label_singleR")]] <- SingleR.lt$pruned.labels
  Query_Seurat@meta.data$label_singleR <- ifelse(is.na(Query_Seurat@meta.data$label_singleR), "Unassign", Query_Seurat@meta.data$label_singleR)

  Query_Seurat@misc$CTAnnot$singleR_Scores <- SingleR.lt@listData[["scores"]]
  Query_Seurat@misc$CTAnnot$singleR_Delta <- SingleR.lt@listData[["delta.next"]]

  return(Query_Seurat)
}


#### scmap ####
Run_scmap <- function(Query_Seurat, Reference_Seurat,
                      Set_RefAnnoCol = "Actual_Cell_Type",
                      Set_Threshold = 0.7, Set_NumFeatures = 500, ...) {
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("scmap", quietly = TRUE)) BiocManager::install("scmap"); library(scmap)
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("Seurat", quietly = TRUE)) install.packages("Seurat");library(Seurat)

  # Convert Seurat object to SingleCellExperiment
  ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
  colData(ref_sce)$cluster <- Reference_Seurat@meta.data[[Set_RefAnnoCol]]

  # Perform feature selection and index clusters
  ref_sce <- selectFeatures(ref_sce, n_features = Set_NumFeatures)
  ref_sce <- indexCluster(ref_sce, cluster_col = "cluster")
  index_list <- list(ref_sce = metadata(ref_sce)$scmap_cluster_index)

  # Convert unlabelled Seurat object to SingleCellExperiment
  query_sce <- as.SingleCellExperiment(Query_Seurat)
  rowData(query_sce)$feature_symbol <- rownames(query_sce)

  # Perform projection
  projection_result <- scmapCluster(projection = query_sce, index_list = index_list, threshold = Set_Threshold)
  projection_result_All <- scmapCluster(projection = query_sce, index_list = index_list, threshold = 0)

  # Add predicted cell types to original Seurat object
  Query_Seurat$label_scmap_NoReject <- projection_result_All[["scmap_cluster_labs"]] %>% as.character()
  Query_Seurat$label_scmap <- projection_result[["scmap_cluster_labs"]] %>% as.character()
  Query_Seurat$label_scmap_Scores <- projection_result[["scmap_cluster_siml"]] %>% as.numeric()

  Query_Seurat@misc$CTAnnot$scmap_Scores <- projection_result[["scmap_cluster_siml"]] %>% as.numeric()


  Query_Seurat@meta.data <- Query_Seurat@meta.data %>%
    mutate(label_scmap = if_else(label_scmap == "unassigned", "Unassign", label_scmap),
           label_scmap_NoReject = if_else(label_scmap_NoReject == "unassigned", "Unassign", label_scmap_NoReject))

  return(Query_Seurat)
}



#### SCINA ####
Run_SCINA <- function(Query_Seurat, Reference_Seurat,
                      Set_RefAnnoCol = "Actual_Cell_Type",
                      ExportFolder = getwd(), Export = "",seurat_version = "V5", ...) {
  # Load necessary package
  if (!requireNamespace("SCINA", quietly = TRUE)) install.packages("SCINA"); library(SCINA)

  # Extract representation matrix from Seurat objects
  if(seurat_version == "V5M"){
    exp <- Seurat::GetAssayData(Query_Seurat, assay = "RNA", slot = "data")
  }else{
    exp <- Seurat::GetAssayData(Query_Seurat, assay = "RNA", slot = "counts")
  }


  # Set cell type annotation column as active ident
  Reference_Seurat <- Seurat::SetIdent(Reference_Seurat, value = Set_RefAnnoCol)

  # Find marker genes for each cell type
  markers <- Seurat::FindAllMarkers(Reference_Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers <- markers %>% dplyr::filter(p_val_adj < 0.01)
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(50, wt = avg_log2FC) %>%
    ungroup()

  # Create a list containing marker genes for each cell type
  cell_types <- unique(Reference_Seurat@meta.data[[Set_RefAnnoCol]])
  signatures <- lapply(cell_types, function(ct) {
    top_markers[top_markers$cluster == ct,]$gene
  })
  names(signatures) <- cell_types

  try({
    # Cell type annotation with SCINA
    results_All <- SCINA(exp, signatures, max_iter = 100, convergence_n = 10,
                         convergence_rate = 0.99, sensitivity_cutoff = 0.9,
                         rm_overlap = TRUE, allow_unknown = FALSE,
                         log_file = paste0(ExportFolder, "/", Export, "_SCINA_All.log"))
    results <- SCINA(exp, signatures, max_iter = 100, convergence_n = 10,
                     convergence_rate = 0.99, sensitivity_cutoff = 0.9,
                     rm_overlap = TRUE, allow_unknown = TRUE,
                     log_file = paste0(ExportFolder, "/", Export, "_SCINA.log"))

    # Add predicted cell types to original Seurat object
    Query_Seurat$label_SCINA_NoReject <- results_All$cell_labels
    Query_Seurat$label_SCINA <- results$cell_labels
    Query_Seurat@meta.data <- Query_Seurat@meta.data %>%
      mutate(label_SCINA = if_else(label_SCINA == "unknown", "Unassign", label_SCINA))
  })

  return(Query_Seurat)
}

#### scPred ####
# if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
# # # trace("project_query", edit=TRUE) # new_data <- GetAssayData(new, "data")[shared_features, ] # new_data <- GetAssayData(new, layer = "data")[shared_features, ]
# # Extract the original code of project_query function
# project_query_code <- deparse(body(project_query))
#
# # Find and replace the target line
# modified_project_query_code <- gsub(
#   'GetAssayData\\(new, "data"\\)',
#   'GetAssayData(new, layer = "data")',
#   project_query_code
# )
#
# # Reassemble the modified code into a function body
# new_project_query_body <- paste(modified_project_query_code, collapse = "\n")
#
# # cat(new_project_query_body) # Print the new function body for inspection
#
# # Redefine the project_query function with the modified function body
# # project_query <- eval(parse(text = paste("function(new, reference, max.iter.harmony = 20, recompute_alignment = TRUE, seed = 66, ...) {", new_project_query_body, "}")))
# #NG# assign("project_query", new_project_query_body, envir = .GlobalEnv)
# #NG# body(project_query) <- new_project_query_body
# # eval(parse(text = paste("project_query <- function(new, reference, max.iter.harmony = 20, recompute_alignment = TRUE, seed = 66, ...) {", new_project_query_body, "}")))
# eval(parse(text = paste0("project_query <- function(new, reference, max.iter.harmony = 20, recompute_alignment = TRUE, \n  seed = 66, ...) \n", new_project_query_body, "")))
#
# # trace("project_query", edit=TRUE)


# trace("project_query", edit=TRUE) # new_data <- GetAssayData(new, "data")[shared_features, ] # new_data <- GetAssayData(new, layer = "data")[shared_features, ]
Run_scPred <- function(Query_Seurat, Reference_Seurat,
                       Set_RefAnnoCol = "Actual_Cell_Type", ...) {
  # # Load necessary packages
  # if(!require("devtools")) install.packages("devtools"); library(devtools)
  # if(!require("harmony")) devtools::install_github("immunogenomics/harmony"); library(harmony)
  # if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)

  # Ensure the default assay is set to RNA
  DefaultAssay(Query_Seurat) <- "RNA"
  DefaultAssay(Reference_Seurat) <- "RNA"

  Reference_Seurat <- getFeatureSpace(Reference_Seurat, Set_RefAnnoCol)   ## Get the feature space to train the classifiers
  Reference_Seurat <- trainModel(Reference_Seurat) ## Train the model

  # trace("project_query", edit=TRUE) # new_data <- GetAssayData(new, "data")[shared_features, ] # new_data <- GetAssayData(new, layer = "data")[shared_features, ]
  Query_Seurat <- scPredict(Query_Seurat, Reference_Seurat) #, threshold = Set_scPredict_Thr)
  # Reference_Seurat <- scPredict(Reference_Seurat, Reference_Seurat)

  modify_colnames_scPred <- function(seuratObject, scPredType) {
    colnames(seuratObject@meta.data) <- sapply(colnames(seuratObject@meta.data), function(colname) {
      if (grepl("^sc[pP]red_", colname) && !grepl("_(max|label|no_rejection|prediction)$", colname)) {
        modified_name <- gsub("sc[pP]red_", "", colname) %>% gsub("\\.", " ", .) %>% gsub("_plus", "+", .)
        paste0(modified_name, " ", scPredType, "Score")
      } else {colname}
    })

    label_scPred_col <- paste0("label_", scPredType)
    colnames(seuratObject@meta.data) <- gsub("sc[pP]red_prediction", label_scPred_col, colnames(seuratObject@meta.data))
    colnames(seuratObject@meta.data) <- gsub("sc[pP]red_no_rejection", paste0(label_scPred_col, "_NoReject"), colnames(seuratObject@meta.data))
    colnames(seuratObject@meta.data) <- gsub("sc[pP]red_max", paste0(label_scPred_col, "_Score"), colnames(seuratObject@meta.data))

    seuratObject@meta.data[[label_scPred_col]] <- ifelse( seuratObject@meta.data[[label_scPred_col]] == "unassigned", "Unassign", seuratObject@meta.data[[label_scPred_col]] )

    return(seuratObject)
  }


  Query_Seurat <- modify_colnames_scPred(Query_Seurat,"scPred")
  # Reference_Seurat <- modify_colnames_scPred(Reference_Seurat,"scPred")

  ## Save the scPred cell type scores in misc slot
  score_columns <- grep(" scPredScore$", colnames(Query_Seurat@meta.data), value = TRUE)
  Query_Seurat@misc$CTAnnot$scPred_Scores <- Query_Seurat@meta.data[, score_columns, drop = FALSE]
  ## Remove score_columns from meta.data
  Query_Seurat@meta.data <- Query_Seurat@meta.data[, !colnames(Query_Seurat@meta.data) %in% score_columns]


  return(Query_Seurat)
}





#### CHETAH ####
Run_CHETAH <- function(Query_Seurat, Reference_Seurat,
                       Set_RefAnnoCol = "Actual_Cell_Type", ...) {
  # PMID31226206 CHETAH # https://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("CHETAH", quietly = TRUE)) BiocManager::install("CHETAH"); library(CHETAH)

  # Convert Seurat objects to SingleCellExperiment
  ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  query_sce <- as.SingleCellExperiment(Query_Seurat)

  # Prepare Reference Data
  ref_sce$celltypes <- Reference_Seurat@meta.data[[Set_RefAnnoCol]]

  # Run CHETAH classifier
  query_sce_All <- CHETAHclassifier(input = query_sce, ref_cells = ref_sce, thresh = 0)
  query_sce <- CHETAHclassifier(input = query_sce, ref_cells = ref_sce)

  # Extract cell types
  celltypes_All <- query_sce_All$celltype_CHETAH
  celltypes <- query_sce$celltype_CHETAH

  # Rename unassigned cell
  celltypes_All <- ifelse(grepl("^Node", celltypes_All), "Unassign", celltypes_All)
  celltypes_All <- ifelse(celltypes_All == "Unassigned", "Unassign", celltypes_All)
  celltypes <- ifelse(grepl("^Node", celltypes), "Unassign", celltypes)
  celltypes <- ifelse(celltypes == "Unassigned", "Unassign", celltypes)

  # Update Seurat Object
  Query_Seurat$label_CHETAH_NoReject <- celltypes_All
  Query_Seurat$label_CHETAH <- celltypes

  return(Query_Seurat)
}


#### scClassify ####
Run_scClassify <- function(Query_Seurat, Reference_Seurat,
                           Set_RefAnnoCol = "Actual_Cell_Type", seurat_version = "V5", ...) {
  # PMID32567229 scClassify  # https://sydneybiox.github.io/scClassify/articles/pretrainedModel.html
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
  if (!require("scClassify", quietly = TRUE)) devtools::install_github("SydneyBioX/scClassify"); library(scClassify)


  # Convert Seurat objects to SingleCellExperiment
  ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  query_sce <- as.SingleCellExperiment(Query_Seurat)

  # Prepare Reference Data
  ref_sce$celltypes <- Reference_Seurat@meta.data[[Set_RefAnnoCol]]


  # Run scClassify
  if(seurat_version == "V5M"){
    exprsMat_Qry <- as.matrix(query_sce@assays@data@listData[["logcounts"]])
    exprsMat_Ref <- as.matrix(counts(ref_sce))

  }else{
    exprsMat_Qry <- as.matrix(counts(query_sce))
    exprsMat_Ref <- as.matrix(counts(ref_sce))
  }

  result_All <- scClassify( exprsMat_train = exprsMat_Ref, cellTypes_train = ref_sce$celltypes,
                            prob_threshold = 0, cor_threshold_static =0 , pSig = 0.1, # pSig = 0.05,
                            exprsMat_test = exprsMat_Qry)
  result_All.df <- as.data.frame(result_All[["testRes"]][["test"]][["pearson_WKNN_limma"]][["predLabelMat"]])

  result <- scClassify( exprsMat_train = exprsMat_Ref, cellTypes_train = ref_sce$celltypes,
                        exprsMat_test = exprsMat_Qry)
  result.df <- as.data.frame(result[["testRes"]][["test"]][["pearson_WKNN_limma"]][["predLabelMat"]])

  # Add the predicted cell types to the Seurat object
  query_sce$label_scClassify_NoReject <- ifelse(result_All.df[,ncol(result_All.df)] == "unassigned", "Unassign", result_All.df[,ncol(result_All.df)])
  query_sce$label_scClassify <- ifelse(result.df[,ncol(result.df)] == "unassigned", "Unassign", result.df[,ncol(result.df)])

  # Update Seurat object
  Query_Seurat$label_scClassify_NoReject <- query_sce$label_scClassify_NoReject
  Query_Seurat$label_scClassify <- query_sce$label_scClassify

  return(Query_Seurat)
}



#### Seurat ####
Run_Seurat_Annot <- function(Query_Seurat, Reference_Seurat,
                             Set_RefAnnoCol = "Actual_Cell_Type",
                             Set_NumPC = 50, ...) {
  # Mapping and annotating query datasets # https://satijalab.org/seurat/articles/integration_mapping
  # Load necessary packages
  if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat"); library(Seurat)
  if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse"); library(tidyverse)
  if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret"); library(caret)

  # Find transfer anchors and transfer data
  anchors <- FindTransferAnchors(reference = Reference_Seurat, query = Query_Seurat, dims = 1:Set_NumPC)
  predictions <- TransferData(anchorset = anchors, refdata = Reference_Seurat@meta.data[[Set_RefAnnoCol]], dims = 1:Set_NumPC)
  Query_Seurat <- AddMetaData(Query_Seurat, metadata = predictions)

  # Rename the predicted.id column to label_Seurat_NoReject
  colnames(Query_Seurat@meta.data)[which(colnames(Query_Seurat@meta.data) == "predicted.id")] <- "label_Seurat_NoReject"

  # # Calculate the maximum prediction score for each cell
  # Query_Seurat$mapping.score <- apply(Query_Seurat@meta.data[, grep("prediction.score", colnames(Query_Seurat@meta.data))], 1, max)
  Query_Seurat$mapping.score <- Query_Seurat$prediction.score.max

  # Add label_Seurat column based on label_Seurat_NoReject but mark cells with mapping.score < 0.8 as Unassign
  Query_Seurat$label_Seurat <- ifelse(Query_Seurat$mapping.score < 0.8, "Unassign", Query_Seurat$label_Seurat_NoReject)
  colnames(Query_Seurat@meta.data)[which(colnames(Query_Seurat@meta.data) == "mapping.score")] <- "label_Seurat_Scores"


  # Save the Seurat cell type scores in misc slot
  Query_Seurat@misc$CTAnnot$Seurat_Scores <- Query_Seurat@meta.data[, grep("prediction.score", colnames(Query_Seurat@meta.data))]
  Query_Seurat@misc$CTAnnot$Seurat_Scores <- Query_Seurat@misc$CTAnnot$Seurat_Scores[, !grepl("prediction.score.max", colnames(Query_Seurat@misc$CTAnnot$Seurat_Scores))]
  colnames(Query_Seurat@misc$CTAnnot$Seurat_Scores) <- gsub("^prediction\\.score\\.", "", colnames(Query_Seurat@misc$CTAnnot$Seurat_Scores))
  colnames(Query_Seurat@misc$CTAnnot$Seurat_Scores) <- gsub("\\.", " ", colnames(Query_Seurat@misc$CTAnnot$Seurat_Scores))

  # Remove prediction.score columns except prediction.score.max
  score_columns <- grep("^prediction.score", colnames(Query_Seurat@meta.data), value = TRUE)
  Query_Seurat@meta.data <- Query_Seurat@meta.data[, !colnames(Query_Seurat@meta.data) %in% score_columns]

  return(Query_Seurat)
}

# ## Test Function
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")
#
# seuratObject_Sample <- pbmc.rna
# seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
#
# seuratObject_Sample <- Run_Seurat_Annot(seuratObject_Sample, seuratObject_Ref)
#
# p1 <- DimPlot(seuratObject_Sample, group.by = "label_Seurat", label = FALSE, label.size = 3) # + NoLegend()
# p2 <- DimPlot(seuratObject_Sample, group.by = "label_Seurat_NoReject")
# # p3 <- DimPlot(seuratObject_Sample, group.by = "seurat_annotations")
# p1 + p2
#
# # Visualize mapping quality score
# FeaturePlot(seuratObject_Sample, features = "mapping.score")


# #### Seurat_Annot_Multimodal ####
#
# Run_Seurat_Annot_Multimodal <- function(seuratObject_Sample, seuratObject_Ref){
#   # Seurat v4 Reference Mapping
#   # https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
#
# }
#






################################################################################
#### scReClassify ####
Fun_scReClassify <- function(Query_Seurat, Reference_Seurat,
                             Set_classifier = "svm", Set_percent = 1,
                             Set_L = 10, ...) {
  # PMID31874628 scReClassify # https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html

  # Load necessary packages
  if (!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)
  if (!require("DT")) install.packages("DT"); library(DT)
  if (!require("mclust")) install.packages("mclust"); library(mclust)

  # Convert Seurat objects to SingleCellExperiment
  # ref_sce <- as.SingleCellExperiment(Reference_Seurat)
  sample_sce <- as.SingleCellExperiment(Query_Seurat)

  # Standardize the data and create 'logNorm' assay
  # Log-normalize dataset if not already done
  if (is.null(Reference_Seurat@assays[["RNA"]]@layers[["data"]])) { ref_sce <- logNormCounts(ref_sce)}
  if (is.null(Query_Seurat@assays[["RNA"]]@layers[["data"]])) { sample_sce <- logNormCounts(sample_sce)}


  # # Create 'logNorm' layer in assays
  # SummarizedExperiment::assay(ref_sce, "logNorm") <- SummarizedExperiment::assay(ref_sce, "logcounts")
  SummarizedExperiment::assay(sample_sce, "logNorm") <- SummarizedExperiment::assay(sample_sce, "logcounts")

  # Remove constant rows
  remove_constant_rows <- function(mat) { mat[rowSums(mat != 0) > 0, ]}

  # logNorm_ref <- remove_constant_rows(assay(ref_sce, "logNorm"))
  logNorm_sample <- remove_constant_rows(assay(sample_sce, "logNorm"))


  # ## Dimension reduction
  # pca_ref <- stats::prcomp(t(logNorm_ref), center = TRUE, scale. = TRUE) # reducedDim(ref_sce, "matPCs") <- matPCs(ref_sce, assay = "logNorm", 0.7)
  pca_sample <- stats::prcomp(t(logNorm_sample), center = TRUE, scale. = TRUE) # reducedDim(sample_sce, "matPCs") <- matPCs(sample_sce, assay = "logNorm", 0.7)

  # reducedDim(ref_sce, "matPCs") <- pca_ref$x
  reducedDim(sample_sce, "matPCs") <- pca_sample$x

  # Cell types
  cellTypes <- sample_sce[["cellTypes"]]


  # Check for NA values and remove them from ref_sce
  valid_cells <- !is.na(cellTypes) & complete.cases(reducedDim(sample_sce, "matPCs"))
  cellTypes <- cellTypes[valid_cells]
  sample_sce <- sample_sce[, valid_cells]

  # Check for balanced classes and adjust if necessary in sample_sce
  cellType_counts <- table(cellTypes)
  LimNum <- 5
  if (any(cellType_counts < LimNum)) {
    warning(paste0("Some classes in reference have fewer than ",LimNum," samples. This may cause issues with the classifier."))
    min_class <- names(cellType_counts[cellType_counts < LimNum])
    for (cls in min_class) {
      sample_sce <- sample_sce[, cellTypes != cls]
      cellTypes <- cellTypes[cellTypes != cls]
    }
  }


  # Run scReClassify
  set.seed(123)
  cellTypes.reclassify <- multiAdaSampling(sample_sce, cellTypes, reducedDimName = "matPCs",
                                           classifier = Set_classifier, percent = Set_percent, L = Set_L)


  # Save the new annotations to meta.data
  Query_Seurat$label_scReClassify <- cellTypes.reclassify$final


  return(Query_Seurat)
}

# # Example usage
# Query_Seurat <- seuratObject_Sample
# Query_Seurat@meta.data$cellTypes <- Query_Seurat@meta.data$seurat_annotations
# # Query_Seurat@meta.data$cellTypes <- gsub("_", "  ", Query_Seurat@meta.data$cellTypes)
#
# Reference_Seurat <- seuratObject_Ref
# seuratObject_Sample_T <- Fun_scReClassify(seuratObject_Sample, seuratObject_Ref)
#
