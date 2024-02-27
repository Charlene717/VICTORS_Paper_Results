##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
library(data.table)



#### Load data ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_GSE136831_COPD/#_AllMix_Export_GSE136831_Sum/20240110193729_GSE136831_Sum_SVDLRglmnet_ROC.RData")
# load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231209_GSE132044_OriQC_DeBug/#_AllMix_Export_GSE132044_OriQC_DeBug_Sum/20240101035050_GSE132044_Sum_SVDLRglmnet_ROC.RData")

folder_path <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240111_IntCT_scRNAseqPanc_Muraro/#_AllMix_Export_20240111_IntCT_scRNAseqPanc_Baron_Muraro_Sum/"
load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20240111_IntCT_scRNAseqPanc_Muraro/#_AllMix_Export_20240111_IntCT_scRNAseqPanc_Baron_Muraro_Sum/20240111113332_Muraro_Baron_Sum_SVDLRglmnet_ROC.RData")


# folder_path <-"D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_GSE135893_GSE128033_PropSame/#_AllMix_Export_GSE135893_GSE128033_Sum/"
# load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_GSE135893_GSE128033_PropSame/#_AllMix_Export_GSE135893_GSE128033_Sum/20240115025955_GSE135893_GSE128033_Sum_SVDLRglmnet_ROC.RData")


# all_data$Actual_Cell_Type %>% unique()

all_data_ori <- all_data

if(grepl("GSE135893", folder_path)) {

}else{
  all_data <- all_data %>%
    filter(!(Sample_Platform == Ref_Platform) &
             !(grepl("^10x", Sample_Platform) & grepl("^10x", Ref_Platform)))
}



#### Set Parameter ####
# Set_Obs_CellType <- "CD4+ T cell" # "CD4+ T cell" # "Cytotoxic T cell"
Set_Obs_CellType <- "Alpha cell" # "Beta cell" # "Alpha cell" # "Acinar cell"

# Set_Obs_CellType <- "Epithelial cell" # "Myeloid cell" # "Natural killer T cell"

Set_Ref_State <- "with" # FALSE  # TRUE # "with" #"lack" #"Comp

if(Set_Ref_State == "with" ){
  Set_Title_End <- paste0(" (Ref with ",Set_Obs_CellType,")")
}else if (Set_Ref_State == "Comp"){
  Set_Title_End <- paste0(" (Ref with all cell type",")")
}else{
  Set_Title_End <- paste0(" (Ref lack ",Set_Obs_CellType,")")
}

if(Set_Ref_State == "with" ){
  Set_Title_End2 <- paste0(" when reference with",Set_Obs_CellType)
}else if (Set_Ref_State == "Comp"){
  Set_Title_End <- paste0(" when reference with all cell type")
}else{
  Set_Title_End2 <- paste0(" when reference lack ",Set_Obs_CellType)
}

# Name_Note <- paste0(Name_Note,"_",Set_Obs_CellType,Set_Title_End)
Name_Note <- paste0(Set_Obs_CellType,Set_Title_End)

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_ExportFolder <- folder_path # Name_ExportFolder <- paste0("Export_IntegAll")
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
# Name_ExportFolder <- paste0(Name_ExportFolder,"/",Name_time_wo_micro)
# if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
Name_Export <- paste0(Name_time_wo_micro,"_",Name_Note)



# 加載所需的套件
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("tidyr")) install.packages("tidyr"); library(tidyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)


# 處理數據集，避免重複的數據並選擇相關列
selected_data <- all_data %>%
  select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType, Actual_Cell_Type, contains("DiagPara"), contains("label")) %>%
  distinct()

# Replace 'unknown' with 'Unassign' in the label_SCINA_NoReject column
selected_data <- selected_data %>%
  mutate(label_SCINA_NoReject = ifelse(label_SCINA_NoReject == "unknown", "Unassign", label_SCINA_NoReject))

# 篩選出 Actual_Cell_Type 為 Set_Obs_CellType 的數據
if(Set_Ref_State == "with" ){
  alpha_cell_data <- selected_data %>%
    filter(Mislabel_CellType != Set_Obs_CellType) %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
}else if (Set_Ref_State == "Comp"){
    filter(Mislabel_CellType == "None") %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
}else{
  alpha_cell_data <- selected_data %>%
    filter(Mislabel_CellType == Set_Obs_CellType) %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
}


# Reshape the data for visualization
alpha_cell_long <- alpha_cell_data %>%
  select(Actual_Cell_Type, label_singleR_NoReject, label_singleR, label_scmap_NoReject, label_scmap, label_SCINA_NoReject, label_SCINA, label_scPred_NoReject, label_scPred) %>%
  pivot_longer(cols = -Actual_Cell_Type, names_to = "Labeling_Method", values_to = "Predicted_Cell_Type")

# Modify labeling methods and order them
alpha_cell_long$Labeling_Method <- gsub("label_", "", alpha_cell_long$Labeling_Method)
alpha_cell_long$Labeling_Method <- factor(alpha_cell_long$Labeling_Method, levels = c("singleR_NoReject", "singleR", "scmap_NoReject", "scmap", "SCINA_NoReject", "SCINA", "scPred_NoReject", "scPred"))


# 創建堆疊條形圖
source("Set_plot_color.R")
ggplot(alpha_cell_long, aes(x = Labeling_Method, fill = Predicted_Cell_Type)) +
  geom_bar() +
  theme_minimal() +
  labs(title = paste0("Labeling across Different Methods on ",Set_Obs_CellType,Set_Title_End), x = "Labeling Method", y = "Count", fill = "Predicted Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the stacked bar chart with the defined color scheme
Plot_bar <- ggplot(alpha_cell_long, aes(x = Labeling_Method, fill = Predicted_Cell_Type)) +
  geom_bar() +
  scale_fill_manual(values = unlist(color_CellType)) +
  theme_minimal() +
  labs(title = paste0("Labeling across Different Methods on ",Set_Obs_CellType,Set_Title_End), x = "Labeling Method", y = "Count", fill = "Predicted Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x axis labels size
        axis.text.y = element_text(size = 16), # Increase y axis labels size
        axis.title = element_text(size = 18), # Increase axis title size
        plot.title = element_text(size = 18), # Increase plot title size
        legend.title = element_text(size = 14), # Increase legend title size
        legend.text = element_text(size = 12)) # Increase legend text size-> Plot_bar

Plot_bar

#################################################################################

#### Heatamp ####
# 加載所需的套件
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("tidyr")) install.packages("tidyr"); library(tidyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

# Calculate count for each labeling method and predicted cell type
cell_count <- alpha_cell_long %>%
  group_by(Labeling_Method, Predicted_Cell_Type) %>%
  summarise(Count = n())

# Create the heatmap with increased font sizes
Plot_Heatmap <- ggplot(cell_count, aes(x = Labeling_Method, y = Predicted_Cell_Type, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "#c3d8eb", high = "steelblue") +
  # scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal(base_size = 14) + # Increase base font size
  # labs(title = "Cell Type Prediction Frequency by Labeling Method", x = "Labeling Method", y = "Predicted Cell Type", fill = "Frequency") +
  labs(title = paste0("Cell Type Prediction Frequency on ",Set_Obs_CellType," by Labeling Method",Set_Title_End), x = "Labeling Method", y = "Predicted Cell Type", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x axis labels size
        axis.text.y = element_text(size = 12), # Increase y axis labels size
        axis.title = element_text(size = 14), # Increase axis title size
        plot.title = element_text(size = 16), # Increase plot title size
        legend.title = element_text(size = 14), # Increase legend title size
        legend.text = element_text(size = 12)) # Increase legend text size
Plot_Heatmap




########################################################################################
#### Bubble Plot ####
# Load necessary libraries
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

# Define the color scheme for classes
color_Class <- list(
  "TP" = "#b58b2a",
  "TN" = "#e8bc56",
  "FN" = "#368a5b",
  "FP" = "#73bd94",
  "Other" = "#919191"
)

# Assuming selected_data is already loaded into your R environment

# # Filter data for Set_Obs_CellType
# alpha_cell_data <- selected_data %>%
#   filter(Mislabel_CellType == Set_Obs_CellType) %>%
#   filter(Actual_Cell_Type == Set_Obs_CellType)


# Prepare and reshape the data for each method
prepare_data <- function(data, method, diag_methods, Set_CellType, Set_CellType_Reverse) {
  if(Set_CellType_Reverse == "with"){
    data %>%
      filter(Actual_Cell_Type == Set_CellType) %>%
      filter(Mislabel_CellType != Set_CellType) %>%
      select(Predicted_Cell_Type = {{method}}, all_of(diag_methods)) %>%
      pivot_longer(cols = -Predicted_Cell_Type, names_to = "Diag_Method", values_to = "Class") %>%
      mutate(Method_Type = as.character(substitute(method)))
    }else if(Set_CellType_Reverse == "Comp"){
      data %>%
        filter(Actual_Cell_Type == Set_CellType) %>%
        filter(Mislabel_CellType == "None") %>%
        select(Predicted_Cell_Type = {{method}}, all_of(diag_methods)) %>%
        pivot_longer(cols = -Predicted_Cell_Type, names_to = "Diag_Method", values_to = "Class") %>%
        mutate(Method_Type = as.character(substitute(method)))
    }else{
      data %>%
        filter(Actual_Cell_Type == Set_CellType) %>%
        filter(Mislabel_CellType == Set_CellType) %>%
        select(Predicted_Cell_Type = {{method}}, all_of(diag_methods)) %>%
        pivot_longer(cols = -Predicted_Cell_Type, names_to = "Diag_Method", values_to = "Class") %>%
        mutate(Method_Type = as.character(substitute(method)))
    }

  }

singleR_data <- prepare_data(selected_data, label_singleR_NoReject, c("label_singleR_DiagPara", "DiagPara_label_singleR_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
scmap_data <- prepare_data(selected_data, label_scmap_NoReject, c("label_scmap_DiagPara", "DiagPara_label_scmap_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
SCINA_data <- prepare_data(selected_data, label_SCINA_NoReject, c("label_SCINA_DiagPara", "DiagPara_label_SCINA_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
scPred_data <- prepare_data(selected_data, label_scPred_NoReject_Annot, c("label_scPred_DiagPara_Annot", "DiagPara_label_scPred_NoReject_Annot_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)

# Combine all method data
combined_data <- bind_rows(singleR_data, scmap_data, SCINA_data, scPred_data) %>%
  group_by(Method_Type, Diag_Method, Predicted_Cell_Type, Class) %>%
  summarise(Count = n(), .groups = 'drop')

# Rename Diag_Method names and reorder them
combined_data <- combined_data %>%
  mutate(Diag_Method = case_when(
    Diag_Method == "DiagPara_label_SCINA_NoReject_SVGLRglmnet_ROC" ~ "SCINA_VICTORS",
    Diag_Method == "DiagPara_label_scPred_NoReject_Annot_SVGLRglmnet_ROC" ~ "scPred_VICTORS",
    Diag_Method == "DiagPara_label_scmap_NoReject_SVGLRglmnet_ROC" ~ "scmap_VICTORS",
    Diag_Method == "DiagPara_label_singleR_NoReject_SVGLRglmnet_ROC" ~ "singleR_VICTORS",
    Diag_Method == "label_SCINA_DiagPara" ~ "SCINA",
    Diag_Method == "label_scPred_DiagPara_Annot" ~ "scPred",
    Diag_Method == "label_scmap_DiagPara" ~ "scmap",
    Diag_Method == "label_singleR_DiagPara" ~ "singleR",
    TRUE ~ Diag_Method
  )) %>%
  mutate(Diag_Method = factor(Diag_Method, levels = c("singleR", "singleR_VICTORS", "scmap", "scmap_VICTORS", "SCINA", "SCINA_VICTORS", "scPred", "scPred_VICTORS")))

# Create a Bubble Chart
ggplot(combined_data, aes(x = Diag_Method, y = Predicted_Cell_Type, size = Count, color = Class)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +
  scale_color_manual(values = unlist(color_Class)) +
  theme_minimal(base_size = 14) +
  # labs(title = paste0(Set_Title_End, ": Distribution of TP, FP, FN, TN in Predicted Cell Types Across Methods"),
  # labs(title = paste0("Distribution of Metrics in Predicted Cell Types for ",Set_Obs_CellType," Across Methods",Set_Title_End),
  labs(title = paste0("Annotation on ",Set_Obs_CellType, Set_Title_End2),

       x = "Diagnostic Method", y = "Predicted Cell Type", size = "Count", color = "Class") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        aspect.ratio = 1,
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) -> Plot_Bubble
Plot_Bubble

filtered_dataframe_FP <- combined_data %>% filter(Class == 'FP')



## Filter data for specific cell types
# combined_data2 <- combined_data %>%
#   filter(Predicted_Cell_Type %in% c("Gamma cell", "Epsilon cell", "Delta cell", "Beta cell")) %>%
#   filter(Class %in% c("TN", "FP", "TP", "FN"))

combined_data2 <- combined_data %>%
  filter(!Predicted_Cell_Type %in% c("Unassign"))

ggplot(combined_data2, aes(x = Diag_Method, y = Predicted_Cell_Type, size = Count, color = Class)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +
  scale_color_manual(values = unlist(color_Class)) +
  theme_minimal(base_size = 14) +
  # labs(title = paste0(Set_Title_End, ": Distribution of TP, FP, FN, TN in Predicted Cell Types Across Methods"),
  # labs(title = paste0("Distribution of Metrics in Predicted Cell Types for ",Set_Obs_CellType," Across Methods",Set_Title_End),
  labs(title = paste0("Annotation on ",Set_Obs_CellType, Set_Title_End2),
       x = "Diagnostic Method", y = "Predicted Cell Type", size = "Count", color = "Class") +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    aspect.ratio = 1,
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)) -> Plot_Bubble2
Plot_Bubble2



pdf(paste0(Name_ExportFolder, "/", Name_Export,"_Check.pdf"),
    width = 8, height = 5) #  width = 17, height = 17)

print(Plot_Bubble2)
print(Plot_Bubble)
print(Plot_Heatmap)
print(Plot_bar)

dev.off()

