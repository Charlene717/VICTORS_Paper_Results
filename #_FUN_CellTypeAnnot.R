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
Run_singleR <- function(seuratObject_Sample, seuratObject_Ref) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if(!require("SingleR")) BiocManager::install("SingleR"); library(SingleR)
  if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
  if(!require("scater")) BiocManager::install("scater") ; library(scater)

  ## Prepossessing
  CTFeatures <- as.SingleCellExperiment(seuratObject_Ref)
  CTFeatures$label <- CTFeatures@colData@listData[[Set_RefAnnoCol]]
  # One should normally do cell-based quality control at this point, but for brevity's sake, we will just remove the unlabelled libraries here.
  CTFeatures <- CTFeatures[,!is.na(CTFeatures$label)]

  ## SingleR() expects reference scRNAdatasets.chr to be normalized and log-transformed.
  # library(scuttle)

  # Check if log normalization has been done
  if (is.null(seuratObject_Ref@assays$RNA@counts) || is.null(seuratObject_Ref@assays$RNA@data)) {
    CTFeatures <- logNormCounts(CTFeatures)
  }

  ## Set Target SeuObj
  ## Prepossessing
  scRNA_Tar <- as.SingleCellExperiment(seuratObject_Sample)
  scRNA_Tar <- scRNA_Tar[,colSums(counts(scRNA_Tar)) > 0] # Remove libraries with no counts.
  # Check if log normalization has been done
  if (is.null(seuratObject_Sample@assays$RNA@counts) || is.null(seuratObject_Sample@assays$RNA@data)) {
    scRNA_Tar <- logNormCounts(scRNA_Tar)
  }


  ## Run SingleR
  library(SingleR)
  SingleR.lt <- SingleR(test = scRNA_Tar, ref = CTFeatures, assay.type.test=1,
                        labels = CTFeatures$label , de.method= "classic") # de.method = c("classic", "wilcox", "t")

  return(SingleR.lt)

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
  celltypes <- ifelse(grepl("^Node", celltypes), "Unassigned", celltypes)
  celltypes_All <- ifelse(grepl("^Node", celltypes_All), "Unassigned", celltypes_All)


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




