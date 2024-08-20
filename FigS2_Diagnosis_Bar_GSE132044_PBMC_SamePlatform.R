##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("tidyr")) install.packages("tidyr"); library(tidyr)


if(!require(cowplot)) install.packages("cowplot"); library(cowplot)


#### Set Parameter ####
# Set_Obs_CellType <- "Alpha cell" # "Beta cell" # "Alpha cell" # "Acinar cell"

Set_Ref_State <- "lack" # "with" #"lack" #"Comp"
Set_RmUnassign <- FALSE

Set_TopN <- 1
# 设置参数以指定要比较的方法
Set_Method <- "singleR" # "singleR", "scmap", "SCINA", "scPred", "CHETAH", "scClassify" ,"Seurat"等
Victors_Method <- paste(Set_Method, "VICTOR", sep = "_")

Set_Dataset <- "scRNAseqPanc"
# "GSE132044_SamePlatform" # "GSE132044_CrossPlatform" # "scRNAseqPanc" # "HLCA_core"

#### Load data ####

# load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/Export_IntegrateAll_20240715161840IVL_GSE132044/20240715161840IVL_GSE132044_IntegrateAll.RData")

if (Set_Dataset %in% c("GSE132044_SamePlatform", "GSE132044_CrossPlatform")) {
  folder_path <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/Export_IntegrateAll_20240715161840IVL_GSE132044/"
  load(paste0(folder_path, "20240715161840IVL_GSE132044_IntegrateAll.RData"))
}else if(Set_Dataset == "scRNAseqPanc" ){
  folder_path <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/Export_IntegrateAll_20240715165652CTZ_scRNAseqPanc/"
  load(paste0(folder_path, "20240715165652CTZ_scRNAseqPanc_IntegrateAll.RData"))
}else if(Set_Dataset == "HLCA_core" ){
  folder_path <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_IntegrateAll_20240817172921OSN_HLCA_core/"
  load(paste0(folder_path, "IntegrateAll_20240817172921OSN_HLCA_core.RData"))

}


if (Set_Dataset %in% c("GSE132044_SamePlatform", "HLCA_core")){
  Set_Platform <- "same"
}else if(Set_Dataset %in% c("scRNAseqPanc","GSE132044_CrossPlatform")){
  Set_Platform <- "cross"
}else{
  Set_Platform <- "other"
}



all_data <- combined_data
all_data_ori <- all_data

if (Set_Platform == "same") {
  Summary.df <- data_same_platform
} else if (Set_Platform == "cross") {
  Summary.df <- data_cross_platform
}


## Remove Nonw
all_data <- all_data %>%  filter(Mislabel_CellType != "None")
Summary.df <- Summary.df %>%  filter(Mislabel_CellType != "None")



#### Data prepocessing ####
# 處理數據集，避免重複的數據並選擇相關列
selected_data <- all_data %>%
  select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType, Actual_Cell_Type, contains("DiagPara"), contains("label")) %>%
  distinct()


## 設定 Same_platform 和 Cross_platform
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
# 定义 platform_group 函数
platform_group <- function(platform) {
  if (platform %in% c("10xV2", "10xV2A", "10xV2B")) {
    return("10xV2_Group")
  } else {
    return(platform)
  }
}

# 根据 Sample_Platform 和 Ref_Platform 分组
selected_data <- selected_data %>%
  mutate(Sample_Platform_Group = sapply(Sample_Platform, platform_group),
         Ref_Platform_Group = sapply(Ref_Platform, platform_group))

# Set_Platform 逻辑
if (Set_Platform == "same") {
  # 筛选 Same_platform
  selected_data <- selected_data %>%
    filter(Sample_Platform_Group == Ref_Platform_Group)
} else if (Set_Platform == "cross") {
  # 筛选 Cross_platform
  selected_data <- selected_data %>%
    filter(Sample_Platform_Group != Ref_Platform_Group)
}else{
  selected_data <- selected_data
}


library(dplyr)
library(tidyr)

# 第一步：找到差值最大的Mislabel_CellType
max_diff_cell_type <- Summary.df %>%
  filter(Metric == "Accuracy") %>%
  filter(Method %in% c(Set_Method, Victors_Method)) %>%
  group_by(Mislabel_CellType, Method) %>%
  summarise(Median_Value = median(Value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Method, values_from = Median_Value) %>%
  mutate(Difference = !!sym(Victors_Method) - !!sym(Set_Method)) %>%
  filter(Difference == max(Difference)) %>%
  pull(Mislabel_CellType)

Set_Obs_CellType <- max_diff_cell_type %>% as.character()

# 第二步：根据条件筛选TargetCell_data
TargetCell_data <- selected_data %>%
  filter(Mislabel_CellType == Set_Obs_CellType)

# # 第三步：找到对应Mislabel_CellType且Method为Set_Method的Value最小的FileID
# min_value_fileID <- Summary.df %>%
#   filter(Metric == "Accuracy", Mislabel_CellType %in% Set_Obs_CellType, # Actual_Cell_Type %in% Set_Obs_CellType,
#          Method == Set_Method) %>%
#   arrange(Value) %>%
#   slice(1) %>%
#   pull(FileID)

################################################################################
library(dplyr)
library(tidyr)  # for pivot_wider

# 创建方法的映射，用于转换 Method 名称
method_mapping <- c(
  "label_singleR_ConfStat" = "singleR",
  "label_scmap_ConfStat" = "scmap",
  "label_SCINA_ConfStat" = "SCINA",
  "label_scPred_ConfStat" = "scPred",
  "label_CHETAH_ConfStat" = "CHETAH",
  "label_scClassify_ConfStat" = "scClassify",
  "label_Seurat_ConfStat" = "Seurat",
  "ConfStat_VICTOR_label_singleR_NoReject" = "singleR_VICTOR",
  "ConfStat_VICTOR_label_scmap_NoReject" = "scmap_VICTOR",
  "ConfStat_VICTOR_label_SCINA_NoReject" = "SCINA_VICTOR",
  "ConfStat_VICTOR_label_scPred_NoReject" = "scPred_VICTOR",
  "ConfStat_VICTOR_label_CHETAH_NoReject" = "CHETAH_VICTOR",
  "ConfStat_VICTOR_label_scClassify_NoReject" = "scClassify_VICTOR",
  "ConfStat_VICTOR_label_Seurat_NoReject" = "Seurat_VICTOR"
)

# 计算 Accuracy 并合并到 Summary.df2
Summary.df2 <- lapply(names(method_mapping), function(method) {
  TargetCell_data %>%
    group_by(FileID, Actual_Cell_Type) %>%
    summarise(
      Value = sum(get(method) %in% c("TP", "TN"), na.rm = TRUE) /
        sum(get(method) %in% c("TP", "TN", "FP", "FN"), na.rm = TRUE),
      Method = method_mapping[[method]]
    )
}) %>%
  bind_rows() %>%
  ungroup()


# # 检查结果
# head(Summary.df2)


# # # 设置 Victors_Method 和 Set_Method
# # Victors_Method <- "singleR_VICTOR"  # 替换为您的实际 VICTORS 方法名称
# # Set_Method <- "singleR"  # 替换为您的实际 Set 方法名称
#
# # # 计算 Victors_Method 和 Set_Method 之间的差异，并找出最大的差异对应的 FileID 和 Actual_Cell_Type
# # max_diff <- Summary.df2 %>%
# #   filter(Method %in% c(Victors_Method, Set_Method)) %>%
# #   pivot_wider(names_from = Method, values_from = Value) %>%
# #   mutate(Difference = !!sym(Victors_Method) - !!sym(Set_Method)) %>%
# #   arrange(desc(Difference)) %>%
# #   slice(1) %>%
# #   select(FileID, Actual_Cell_Type, Difference)
#
# # # 显示结果
# # max_diff

# # 设置要获取的名次
# Set_TopN <- 3  # 设置为 2, 3 或其他值，获取对应名次

# 计算每个 Actual_Cell_Type 在 TargetCell_data 中的出现次数
type_count <- TargetCell_data %>%
  group_by(Actual_Cell_Type) %>%
  summarise(Count = n())

# 计算 Victors_Method 和 Set_Method 之间的差异，并找到指定名次的情况
max_diff <- Summary.df2 %>%
  filter(Method %in% c(Victors_Method, Set_Method)) %>%
  pivot_wider(names_from = Method, values_from = Value) %>%
  mutate(Difference = !!sym(Victors_Method) - !!sym(Set_Method)) %>%
  left_join(type_count, by = "Actual_Cell_Type") %>%
  arrange(desc(Difference), desc(Count)) %>%
  slice(Set_TopN) %>%  # 获取指定的第 N 名
  select(FileID, Actual_Cell_Type, Difference, Count)

# # 显示结果
# max_diff


max_value_fileID <- max_diff$FileID

################################################################################

# 第四步：使用FileID筛选TargetCell_data
TargetCell_data <- TargetCell_data %>%
  filter(FileID %in% max_value_fileID)

if(Set_Method == "SCINA"){
  TargetCell_data <- TargetCell_data %>% filter(label_SCINA_ConfStat != "Other")
}
if(Set_Method == "scmap"){
  TargetCell_data <- TargetCell_data %>% filter(label_scmap_ConfStat != "Other")
}

# Reshape the data for visualization
TargetCell_long <- TargetCell_data %>%
  select(Actual_Cell_Type, label_singleR_NoReject, label_singleR, label_scmap_NoReject, label_scmap, label_SCINA_NoReject, label_SCINA, label_scPred_NoReject, label_scPred,
         label_CHETAH_NoReject, label_CHETAH, label_scClassify_NoReject, label_scClassify, label_Seurat_NoReject, label_Seurat) %>%
  pivot_longer(cols = -Actual_Cell_Type, names_to = "Labeling_Method", values_to = "Predicted_Cell_Type")



# Modify labeling methods and order them
TargetCell_long$Labeling_Method <- gsub("label_", "", TargetCell_long$Labeling_Method)
TargetCell_long$Labeling_Method <- factor(TargetCell_long$Labeling_Method,
                                          levels = c("singleR_NoReject", "singleR", "scmap_NoReject", "scmap", "SCINA_NoReject", "SCINA", "scPred_NoReject", "scPred",
                                                     "CHETAH_NoReject", "CHETAH","scClassify_NoReject", "scClassify","Seurat_NoReject", "Seurat"))



Platform_Query <- TargetCell_data$Sample_Platform %>% unique()
Platform_Ref <-TargetCell_data$Ref_Platform %>% unique()


# Prepare and reshape the data for each method
prepare_data <- function(data, method, diag_methods, Set_CellType) {
  data %>%
    filter(Actual_Cell_Type == Set_CellType) %>%
    # filter(Mislabel_CellType == Set_CellType) %>%
    select(Predicted_Cell_Type = {{method}}, all_of(diag_methods)) %>%
    pivot_longer(cols = -Predicted_Cell_Type, names_to = "Diag_Method", values_to = "Class") %>%
    mutate(Method_Type = as.character(substitute(method)))
}

singleR_data <- prepare_data(TargetCell_data, label_singleR_NoReject, c("label_singleR_ConfStat", "ConfStat_VICTOR_label_singleR_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
scmap_data <- prepare_data(TargetCell_data, label_scmap_NoReject, c("label_scmap_ConfStat", "ConfStat_VICTOR_label_scmap_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
SCINA_data <- prepare_data(TargetCell_data, label_SCINA_NoReject, c("label_SCINA_ConfStat", "ConfStat_VICTOR_label_SCINA_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
scPred_data <- prepare_data(TargetCell_data, label_scPred_NoReject, c("label_scPred_ConfStat", "ConfStat_VICTOR_label_scPred_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
CHETAH_data <- prepare_data(TargetCell_data, label_CHETAH_NoReject, c("label_CHETAH_ConfStat", "ConfStat_VICTOR_label_CHETAH_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
scClassify_data <- prepare_data(TargetCell_data, label_scClassify_NoReject, c("label_scClassify_ConfStat", "ConfStat_VICTOR_label_scClassify_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)
Seurat_data <- prepare_data(TargetCell_data, label_Seurat_NoReject, c("label_Seurat_ConfStat", "ConfStat_VICTOR_label_Seurat_NoReject"), Set_CellType = max_diff$Actual_Cell_Type)



# Combine all method data
combined_data <- bind_rows(singleR_data, scmap_data, SCINA_data, scPred_data,
                           CHETAH_data,scClassify_data,Seurat_data) %>%
  group_by(Method_Type, Diag_Method, Predicted_Cell_Type, Class) %>%
  summarise(Count = n(), .groups = 'drop')

# Rename Diag_Method names and reorder them
combined_data <- combined_data %>%
  mutate(Diag_Method = case_when(
    Diag_Method == "ConfStat_VICTOR_label_SCINA_NoReject" ~ "SCINA_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_scPred_NoReject" ~ "scPred_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_scmap_NoReject" ~ "scmap_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_singleR_NoReject" ~ "singleR_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_CHETAH_NoReject" ~ "CHETAH_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_scClassify_NoReject" ~ "scClassify_VICTOR",
    Diag_Method == "ConfStat_VICTOR_label_Seurat_NoReject" ~ "Seurat_VICTOR",
    Diag_Method == "label_SCINA_ConfStat" ~ "SCINA",
    Diag_Method == "label_scPred_ConfStat" ~ "scPred",
    Diag_Method == "label_scmap_ConfStat" ~ "scmap",
    Diag_Method == "label_singleR_ConfStat" ~ "singleR",
    Diag_Method == "label_CHETAH_ConfStat" ~ "CHETAH",
    Diag_Method == "label_scClassify_ConfStat" ~ "scClassify",
    Diag_Method == "label_Seurat_ConfStat" ~ "Seurat",
    TRUE ~ Diag_Method
  )) %>%
  mutate(Diag_Method = factor(Diag_Method, levels = c("singleR", "singleR_VICTOR", "scmap", "scmap_VICTOR", "SCINA", "SCINA_VICTOR", "scPred", "scPred_VICTOR",
                                                      "CHETAH","CHETAH_VICTOR", "scClassify", "scClassify_VICTOR", "Seurat","Seurat_VICTOR")))



## Filter data for specific cell types
# combined_data2 <- combined_data %>%
#   filter(Predicted_Cell_Type %in% c("Gamma cell", "Epsilon cell", "Delta cell", "Beta cell")) %>%
#   filter(Class %in% c("TN", "FP", "TP", "FN"))

if(Set_RmUnassign){
  combined_data <- combined_data %>%
    filter(!Predicted_Cell_Type %in% c("Unassign"))
}


combined_data <- combined_data %>%
  filter(Diag_Method %in% c(Set_Method,Victors_Method))

#### 修改combined_data for 作圖美觀 ####

unique_cell_types <- unique(data_cross_platform$Mislabel_CellType)
unique_cell_types <- unique_cell_types[!unique_cell_types == "None"]

# 生成所有可能的 Diag_Method 和 Predicted_Cell_Type 组合
all_combinations <- expand.grid(Diag_Method = c(Set_Method, Victors_Method),
                                Predicted_Cell_Type = unique_cell_types,
                                Class = unique(combined_data$Class), # 确保包含所有类别
                                stringsAsFactors = FALSE)

# 为缺失的组合创建新行，这里假设Class为NA的情况应该改为具体的Class类别
missing_entries <- all_combinations %>%
  mutate(Method_Type = "method", Count = 0) %>% # 添加缺失的条目
  anti_join(combined_data, by = c("Diag_Method", "Predicted_Cell_Type", "Class"))

# 将缺失的条目合并到 combined_data 中
combined_data <- rbind(combined_data, missing_entries)

# 替换 Diag_Method 中的 "_VICTORS" 为 "_VICTOR"
combined_data$Diag_Method <- gsub("_VICTORS", "_VICTOR", combined_data$Diag_Method)

# Remove "Unassign"
combined_data <- combined_data[!combined_data$Predicted_Cell_Type == "Unassign",]
combined_data <- combined_data %>% filter(!is.na(Class))

#### Bubble Plot ####
# Load necessary libraries
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

if(Set_Ref_State == "with" ){
  Set_Title_End <- paste0(" when reference include ",Set_Obs_CellType)
}else if (Set_Ref_State == "Comp"){
  Set_Title_End <- paste0(" when reference include all cell type")
}else{
  Set_Title_End <- paste0(" when reference lack ",Set_Obs_CellType)
}


# Define the color scheme for classes
color_Class <- list(
  "TP" = "#b58b2a",
  "TN" = "#e8bc56",
  "FN" = "#368a5b",
  "FP" = "#73bd94",
  "Other" = "#919191"
)


# 设定y轴的范围，这里假设为0到最大值的一定比例或具体数值
y_min <- 0
# y_max <- 717 # 示例，根据你的数据调整

if (Set_Dataset %in% c("GSE132044_SamePlatform", "GSE132044_CrossPlatform")) {
  y_max <- 717
}else if(Set_Dataset == "scRNAseqPanc" ){
  y_max <- 326
}else if(Set_Dataset == "HLCA_core" ){
  y_max <- 717

}



## Plot a Bar Chart
library(ggplot2)

# 绘制条形图并有条件地添加文本标签（仅非零值）
ggplot(combined_data, aes(x = Predicted_Cell_Type, y = Count, fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_text(aes(label = ifelse(Count == 0, "", Count), # 如果Count为0，则标签为空字符串，否则为Count值
                y = ifelse(Count == 0, 0, Count + 0.03 * (y_max - y_min))), # 为0的标签不显示
            position = position_dodge(width = 0.7),
            vjust = 0,angle = 25,
            size = 3.5) +
  scale_y_continuous(limits = c(NA, y_max)) +
  scale_fill_manual(values = unlist(color_Class)) +
  facet_wrap(~ Diag_Method, ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(title = paste0("Annotation on ", max_diff$Actual_Cell_Type,"\n", Set_Title_End),
       subtitle = paste("Query:", Platform_Query," Reference:",Platform_Ref),
       x = "Predicted Cell Type", y = "Count", fill = "Class") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        aspect.ratio = 1) -> Plot_Bar

Plot_Bar



filtered_dataframe_FP <- combined_data %>% filter(Class == 'FP')



#### Export ####
# Name_Note <- paste0(Name_Note,"_",Set_Obs_CellType,Set_Title_End)
Name_Note <- paste0(Set_Method,"_",max_value_fileID,"_",Set_Obs_CellType,Set_Title_End)

# Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_ExportFolder <- paste0(folder_path,"Diagnosis_Bar/") # Name_ExportFolder <- paste0("Export_IntegAll")
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
# Name_ExportFolder <- paste0(Name_ExportFolder,"/",Name_time_wo_micro)
# if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder

# Name_Export <- paste0(Name_time_wo_micro,"_",Name_Note)
Name_Export <- paste0(Set_Dataset,"_",Set_Method,"_",max_value_fileID,"_Top",Set_TopN)


## Export PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export,"_Check_Bar.pdf"),
    width = 7, height = 9) #  width = 17, height = 17)

print(Plot_Bar)


dev.off()

## Export TSV
try(write_tsv(combined_data, paste0(Name_ExportFolder, "/", Name_Export, "_Check_Bar.tsv")))
