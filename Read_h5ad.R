##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

library(zellkonverter)

library(reticulate)
use_condaenv("zellkonverterAnnDataEnv")  # 使用適當的 Conda 環境
anndata <- import("anndata")

# 嘗試讀取 .h5ad 檔案
adata <- readH5AD("C:/Users/user/Downloads/2aa90e63-9a6d-444d-8343-8fc2a9921797.h5ad")
