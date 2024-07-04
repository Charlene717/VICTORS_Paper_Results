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

#### singleR ####
Run_singleR <- function(Query_Seurat, Reference_Seurat, Set_RefAnnoCol = "Actual_Cell_Type") {
  # Load necessary packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("SingleR", quietly = TRUE)) BiocManager::install("SingleR"); library(SingleR)
  if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if (!require("scater", quietly = TRUE)) BiocManager::install("scater"); library(scater)


  # Preprocess reference dataset
  CTFeatures <- as.SingleCellExperiment(Reference_Seurat)
  CTFeatures$label <- CTFeatures[[Set_RefAnnoCol]]
  CTFeatures <- CTFeatures[, !is.na(CTFeatures$label)]

  # Log-normalize reference dataset if not already done
  if (is.null(Reference_Seurat@assays[["RNA"]]@layers[["data"]])) {
    CTFeatures <- logNormCounts(CTFeatures)
  }

  # Preprocess query dataset
  scRNA_Tar <- as.SingleCellExperiment(Query_Seurat)
  scRNA_Tar <- scRNA_Tar[, colSums(counts(scRNA_Tar)) > 0]

  # Log-normalize query dataset if not already done
  if (is.null(Query_Seurat@assays[["RNA"]]@layers[["counts"]]) || is.null(Query_Seurat@assays[["RNA"]]@layers[["data"]])) {
    scRNA_Tar <- logNormCounts(scRNA_Tar)
  }

  # Run SingleR
  SingleR.lt <- SingleR(test = scRNA_Tar, ref = CTFeatures, assay.type.test = 1,
                        labels = CTFeatures$label, de.method = "classic")

  Query_Seurat@meta.data$label_singleR_NoReject <- SingleR.lt$labels

  ## Annotation diagnostics
  Query_Seurat@meta.data[[paste0("label_singleR")]] <- SingleR.lt$pruned.labels
  Query_Seurat@meta.data$label_singleR <- ifelse(is.na(Query_Seurat@meta.data$label_singleR), "Unassign", Query_Seurat@meta.data$label_singleR)

  Query_Seurat@misc$CTAnnot$singleR_scores <-SingleR.lt@listData[["scores"]]
  Query_Seurat@meta.data$singleR_delta <-SingleR.lt@listData[["delta.next"]]

  return(Query_Seurat)
}


#### scmap ####
Run_scmap <- function(seuratObject_Sample, seuratObject_Ref) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if(!require("scmap")) BiocManager::install("scmap"); library(scmap)

  if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

  ## Assume seuratObject_Ref is your Seurat object

  ## Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(seuratObject_Ref)

  ## Add gene names to rowData
  rowData(sce)$feature_symbol <- rownames(sce)

  ## Get cell type metadata
  cell_types <- seuratObject_Ref@meta.data$Actual_Cell_Type

  ## Perform feature selection
  sce <- selectFeatures(sce, n_features = 500)

  ## Add cell type labels to sce object
  colData(sce)$cluster <- cell_types

  ## Index clusters
  sce <- indexCluster(sce, cluster_col = "cluster")

  ## Create a list with index
  index_list <- list(sce = metadata(sce)$scmap_cluster_index)


  ## Assume seuratObject_Sample is your unlabelled Seurat object
  ## Convert unlabelled Seurat object to SingleCellExperiment
  sce_unlabelled <- as.SingleCellExperiment(seuratObject_Sample)

  ## Add gene names to rowData
  rowData(sce_unlabelled)$feature_symbol <- rownames(sce_unlabelled)

  ## Perform projection
  projection_result <- scmapCluster(projection = sce_unlabelled, index_list = index_list)
  # sum(projection_result[["scmap_cluster_labs"]] == projection_result[["combined_labs"]])

  projection_result_All <- scmapCluster(projection = sce_unlabelled, index_list = index_list,
                                        threshold = 0)
  # sum(projection_result[["scmap_cluster_siml"]] == projection_result_All[["scmap_cluster_siml"]])


  # Add predicted cell types to original Seurat object
  seuratObject_Sample$label_scmap <- projection_result[["scmap_cluster_labs"]] %>% as.character()
  seuratObject_Sample$label_scmap_Score <- projection_result[["scmap_cluster_siml"]] %>% as.numeric()
  seuratObject_Sample$label_scmap_NoReject <- projection_result_All[["scmap_cluster_labs"]] %>% as.character()

  # DimPlot(seuratObject_Sample, reduction = "umap", group.by ="label_scmap" ,label = TRUE, pt.size = 0.5)# + NoLegend()
  # table(seuratObject_Sample$label_scmap)

  seuratObject_Sample@meta.data <- seuratObject_Sample@meta.data %>%
    mutate(label_scmap = if_else(label_scmap == "unassigned", "Unassign", label_scmap))
  seuratObject_Sample@meta.data <- seuratObject_Sample@meta.data %>%
    mutate(label_scmap_NoReject = if_else(label_scmap_NoReject == "unassigned", "Unassign", label_scmap_NoReject))
  return(seuratObject_Sample)
}



#### SCINA ####
Run_SCINA <- function(seuratObject_Sample, seuratObject_Ref, ExportFolder = getwd(), Export = "") {
  if(!require("SCINA")) install.packages("SCINA"); library(SCINA)

  ## Extract representation matrix from Seurat objects
  exp <- Seurat::GetAssayData(seuratObject_Sample, assay = "RNA", slot = "counts")

  # ## Set "Actual_Cell_Type" as active ident
  # seuratObject_Sample <- SetIdent(seuratObject_Sample, value = "Actual_Cell_Type")

  ## Find marker genes for each cell type
  seuratObject_Ref <- SetIdent(seuratObject_Ref, value = "Actual_Cell_Type")
  markers <- Seurat::FindAllMarkers(seuratObject_Ref, only.pos = TRUE,
                                    min.pct = 0.25, logfc.threshold = 0.25)

  markers <- markers %>% dplyr::filter(p_val_adj < 0.01)
  markers$avg_log2FC <- as.numeric(markers$avg_log2FC)
  top_markers <- markers %>%
    group_by(cluster) %>%
    top_n(50, wt = avg_log2FC) %>%
    ungroup()
  top_markers$cluster <- as.character(top_markers$cluster)

  ## Create a list containing marker genes for each cell type
  cell_types <- unique(seuratObject_Ref@meta.data$Actual_Cell_Type)
  signatures <- lapply(cell_types, function(ct) {
    top_markers[top_markers$cluster == ct,]$gene
  })
  names(signatures) <- cell_types

  try({
    ## Cell type annotation with SCINA
    results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10,
                    convergence_rate = 0.99, sensitivity_cutoff = 0.9,
                    rm_overlap = TRUE, allow_unknown = TRUE,
                    log_file = paste0(ExportFolder,"/",Export,"_SCINA.log"))

    results_All = SCINA(exp, signatures, max_iter = 100, convergence_n = 10,
                        convergence_rate = 0.99, sensitivity_cutoff = 0.9,
                        rm_overlap = TRUE, allow_unknown = FALSE,
                        log_file = paste0(ExportFolder,"/",Export,"_SCINA_All.log"))

    # ## View Results
    # View(results$cell_labels)
    # View(results$probabilities)
    #
    # ## Draw heatmap (may cause crash)
    # plotheat.SCINA(exp, results, signatures)

    ## Save label_SCINA to metadata in seurat object
    seuratObject_Sample@meta.data$label_SCINA <- results$cell_labels
    seuratObject_Sample@meta.data$label_SCINA_NoReject <- results_All$cell_labels
    seuratObject_Sample@meta.data <- seuratObject_Sample@meta.data %>%
      mutate(label_SCINA = if_else(label_SCINA == "unknown", "Unassign", label_SCINA))

  })


  return(seuratObject_Sample)
}


#### CHETAH ####
Run_CHETAH <- function(seuratObject_Sample, seuratObject_Ref,...){
  # PMID31226206 CHETAH
  # https://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if(!require("CHETAH")) BiocManager::install("CHETAH"); library(CHETAH)

  # Convert Seurat objects to SingleCellExperiment
  ref_sce <- as.SingleCellExperiment(seuratObject_Ref)
  sample_sce <- as.SingleCellExperiment(seuratObject_Sample)


  # Prepare Reference Data
  ref_sce$celltypes <- seuratObject_Ref@meta.data[["Actual_Cell_Type"]]


  # Run CHETAH classifier
  sample_sce <- CHETAHclassifier(input = sample_sce, ref_cells = ref_sce)
  sample_sce_All <- CHETAHclassifier(input = sample_sce, ref_cells = ref_sce, thresh = 0)

  # Plot classification
  PlotCHETAH(sample_sce)

  # Extract cell types
  celltypes <- sample_sce$celltype_CHETAH
  celltypes_All <- sample_sce_All$celltype_CHETAH


  # Rename unassigned cell
  celltypes <- ifelse(grepl("^Node", celltypes), "Unassign", celltypes)
  celltypes <- ifelse(celltypes == "Unassigned", "Unassign", celltypes)

  celltypes_All <- ifelse(grepl("^Node", celltypes_All), "Unassign", celltypes_All)
  celltypes_All <- ifelse(celltypes_All == "Unassigned", "Unassign", celltypes_All)

  # Update Seurat Object
  seuratObject_Sample$label_CHETAH <- celltypes
  seuratObject_Sample$label_CHETAH_NoReject <- celltypes_All

  return(seuratObject_Sample)
}

# ## Test function
# seuratObject_Sample <- Run_CHETAH(seuratObject_Sample, seuratObject_Ref)
#
# plot_CHETAH <- DimPlot(seuratObject_Sample,group.by = "label_CHETAH", reduction = "umap")
# plot_seurat <- DimPlot(seuratObject_Sample,group.by = "seurat_annotations", reduction = "umap")
# plot_CHETAH_All <- DimPlot(seuratObject_Sample,group.by = "label_CHETAH_NoReject", reduction = "umap")
#
# plot_seurat + plot_CHETAH + plot_CHETAH_All



#### scClassify ####
Run_scClassify <- function(seuratObject_Sample, seuratObject_Ref){
  # PMID32567229 scClassify
  # https://sydneybiox.github.io/scClassify/articles/pretrainedModel.html
  # https://www.embopress.org/doi/full/10.15252/msb.20199389

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # BiocManager::install(c("S4Vectors", "hopach", "limma"))
  # BiocManager::install(c("Rhdf5lib", "rhdf5filters", "rhdf5", "sparseMatrixStats", "graph", "HDF5Array", "DelayedMatrixStats", "GSEABase", "Cepo"))
  if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)

  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if(!require("scClassify")) devtools::install_github("SydneyBioX/scClassify"); library(scClassify)



  # Convert Seurat objects to SingleCellExperiment
  ref_sce <- as.SingleCellExperiment(seuratObject_Ref)
  sample_sce <- as.SingleCellExperiment(seuratObject_Sample)


  # Prepare Reference Data
  ref_sce$celltypes <- seuratObject_Ref@meta.data[["Actual_Cell_Type"]]


  # Run scClassify
  result <- scClassify(
    exprsMat_train = as.matrix(counts(ref_sce)),
    cellTypes_train = ref_sce$celltypes,
    exprsMat_test = as.matrix(counts(sample_sce))
  )

  result.df <- result[["testRes"]][["test"]][["pearson_WKNN_limma"]][["predLabelMat"]] %>% as.data.frame()

  result_All <- scClassify(
    exprsMat_train = as.matrix(counts(ref_sce)),
    cellTypes_train = ref_sce$celltypes,
    prob_threshold = 0,
    exprsMat_test = as.matrix(counts(sample_sce))
  )

  result_All.df <- result_All[["testRes"]][["test"]][["pearson_WKNN_limma"]][["predLabelMat"]] %>% as.data.frame()


  # Add the predicted cell types to the Seurat object
  sample_sce$label_scClassify <- result.df[,ncol(result.df)]
  sample_sce$label_scClassify_NoReject <- result_All.df[,ncol(result_All.df)]

  # Replace "unassigned" with "Unassign"
  sample_sce$label_scClassify <- ifelse(sample_sce$label_scClassify == "unassigned", "Unassign", sample_sce$label_scClassify)
  sample_sce$label_scClassify_NoReject <- ifelse(sample_sce$label_scClassify_NoReject == "unassigned", "Unassign", sample_sce$label_scClassify_NoReject)

  # Update Seurat object
  seuratObject_Sample$label_scClassify <- sample_sce$label_scClassify
  seuratObject_Sample$label_scClassify_NoReject <- sample_sce$label_scClassify_NoReject

  return(seuratObject_Sample)
}


# ## Test function
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")
#
# seuratObject_Sample <- pbmc.rna
# seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
#
# seuratObject_Sample <- Run_scClassify(seuratObject_Sample, seuratObject_Ref)
#
# plot_scClassify <- DimPlot(seuratObject_Sample,group.by = "label_scClassify", reduction = "umap")
# plot_seurat <- DimPlot(seuratObject_Sample,group.by = "seurat_annotations", reduction = "umap")
# plot_scClassify_All <- DimPlot(seuratObject_Sample,group.by = "label_scClassify_NoReject", reduction = "umap")
#
# plot_seurat + plot_scClassify + plot_scClassify_All



#### Seurat ####
Run_Seurat_Annot <- function(Query_Seurat, Reference_Seurat){
  # Mapping and annotating query datasets
  # https://satijalab.org/seurat/articles/integration_mapping

  #### Load Packages ####
  if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
  if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if(!require("caret")) install.packages("caret"); library(caret)

  ## Find transfer anchors and transfer data
  # Find transfer anchors
  anchors <- FindTransferAnchors(reference = Reference_Seurat, query = Query_Seurat, dims = 1:30)

  # Transfer data
  predictions <- TransferData(anchorset = anchors, refdata = Reference_Seurat$Actual_Cell_Type, dims = 1:30)

  # Add predicted results to the query dataset
  Query_Seurat <- AddMetaData(Query_Seurat, metadata = predictions)

  # Rename the predicted.id column to label_Seurat_NoReject
  colnames(Query_Seurat@meta.data)[which(colnames(Query_Seurat@meta.data) == "predicted.id")] <- "label_Seurat_NoReject"

  ## Mapping QC
  # Calculate the maximum prediction score for each cell
  Query_Seurat$mapping.score <- apply(Query_Seurat@meta.data[, grep("prediction.score", colnames(Query_Seurat@meta.data))], 1, max)

  # Add label_Seurat column based on label_Seurat_NoReject but mark cells with mapping.score < 0.8 as Unassign
  Query_Seurat$label_Seurat <- ifelse(Query_Seurat$mapping.score < 0.8, "Unassign", Query_Seurat$label_Seurat_NoReject)

  return(Query_Seurat)
}


# ## Test Function
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
