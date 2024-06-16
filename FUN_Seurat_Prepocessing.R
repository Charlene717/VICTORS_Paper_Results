Seurat_Prepocessing <- function(seurat_obj, Num_PCA = 50 ,Set_nfeatures = 2000 ) {

  seurat_obj <- seurat_obj  %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = Set_nfeatures) %>%
    ScaleData() %>%
    RunPCA(npcs = Num_PCA) %>%
    FindNeighbors(dims = 1:Num_PCA) %>%
    FindClusters() %>% # resolution = 0.5
    # RunTSNE(dims = 1:Num_PCA) %>%
    RunUMAP(dims = 1:Num_PCA)

}

# # Set_nfeatures <- 2000
# # Num_PCA <- 50
# set.seed(Set_Seed)
# seuratObject<- seuratObject %>%
#   NormalizeData() %>%
#   FindVariableFeatures(nfeatures = Set_nfeatures) %>%
#   ScaleData() %>%
#   RunPCA(npcs = Num_PCA) %>%
#   FindNeighbors(dims = 1:Num_PCA) %>%
#   FindClusters() %>% # resolution = 0.5
#   # RunTSNE(dims = 1:Num_PCA) %>%
#   RunUMAP(dims = 1:Num_PCA)
