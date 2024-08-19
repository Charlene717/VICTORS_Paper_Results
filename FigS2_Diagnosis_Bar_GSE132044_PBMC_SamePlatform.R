##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("tidyr")) install.packages("tidyr"); library(tidyr)


if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

#### Load data ####
load("D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/Fig3_Fig4_GSE132044_PBMC/20240213020617_GSE132044_Sum_SVDLRglmnet_ROC.RData")
folder_path <- "D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/FigS2_GSE132044_PBMC_SamePlatform_Bar/"

# all_data$Actual_Cell_Type %>% unique()
all_data_ori <- all_data


#### Set Parameter ####
# Set_Obs_CellType <- "Alpha cell" # "Beta cell" # "Alpha cell" # "Acinar cell"

Set_Ref_State <- "lack" # "with" #"lack" #"Comp"
Set_RmUnassign <- FALSE

# 设置参数以指定要比较的方法
Set_Method <- "singleR" # 示例，可以改为其他方法，如"singleR", "scmap", "SCINA", "scPred"等
Victors_Method <- paste(Set_Method, "VICTORS", sep = "_")


Summary.df <- data_same_platform
if(Set_Method == "scPred"){
  Summary.df <- Summary.df %>% filter(!FileID == "20231211182309SZHNDI")
}

#### Data prepocessing ####
# 處理數據集，避免重複的數據並選擇相關列
selected_data <- all_data %>%
  select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType, Actual_Cell_Type, contains("DiagPara"), contains("label")) %>%
  distinct()

# Replace 'unknown' with 'Unassign' in the label_SCINA_NoReject column
selected_data <- selected_data %>%
  mutate(label_SCINA_NoReject = ifelse(label_SCINA_NoReject == "unknown", "Unassign", label_SCINA_NoReject))

# 筛选Sample_Platform和Ref_Platform相等的行
selected_data <- selected_data %>%
  filter(Sample_Platform == Ref_Platform |
           (grepl("10x", Sample_Platform) & grepl("10x", Ref_Platform)))

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
# 注意，这里的selected_data应该是你要从中筛选数据的DataFrame
# 这一步需要明确selected_data的来源，这里假设它已经定义好了
if(Set_Ref_State == "with" ){
  TargetCell_data <- selected_data %>%
    filter(Mislabel_CellType != Set_Obs_CellType) %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
} else if (Set_Ref_State == "Comp"){
  TargetCell_data <- selected_data %>%
    filter(Mislabel_CellType == "None") %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
} else {
  TargetCell_data <- selected_data %>%
    filter(Mislabel_CellType == Set_Obs_CellType) %>%
    filter(Actual_Cell_Type == Set_Obs_CellType)
}


# 第三步：找到对应Mislabel_CellType且Method为Set_Method的Value最小的FileID
min_value_fileID <- Summary.df %>%
  filter(Metric == "Accuracy", Mislabel_CellType %in% Set_Obs_CellType, Method == Set_Method) %>%
  arrange(Value) %>%
  slice(1) %>%
  pull(FileID)

# 第四步：使用FileID筛选TargetCell_data
TargetCell_data <- TargetCell_data %>%
  filter(FileID %in% min_value_fileID)

if(Set_Method == "SCINA"){
  TargetCell_data <- TargetCell_data %>% filter(label_SCINA_DiagPara != "Other")
}
if(Set_Method == "scmap"){
  TargetCell_data <- TargetCell_data %>% filter(label_scmap_DiagPara != "Other")
}

# Reshape the data for visualization
TargetCell_long <- TargetCell_data %>%
  select(Actual_Cell_Type, label_singleR_NoReject, label_singleR, label_scmap_NoReject, label_scmap, label_SCINA_NoReject, label_SCINA, label_scPred_NoReject, label_scPred) %>%
  pivot_longer(cols = -Actual_Cell_Type, names_to = "Labeling_Method", values_to = "Predicted_Cell_Type")

# Modify labeling methods and order them
TargetCell_long$Labeling_Method <- gsub("label_", "", TargetCell_long$Labeling_Method)
TargetCell_long$Labeling_Method <- factor(TargetCell_long$Labeling_Method, levels = c("singleR_NoReject", "singleR", "scmap_NoReject", "scmap", "SCINA_NoReject", "SCINA", "scPred_NoReject", "scPred"))



Platform_Query <- TargetCell_data$Sample_Platform %>% unique()
Platform_Ref <-TargetCell_data$Ref_Platform %>% unique()

# # Filter data for Set_Obs_CellType
# TargetCell_data <- selected_data %>%
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

singleR_data <- prepare_data(TargetCell_data, label_singleR_NoReject, c("label_singleR_DiagPara", "DiagPara_label_singleR_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
scmap_data <- prepare_data(TargetCell_data, label_scmap_NoReject, c("label_scmap_DiagPara", "DiagPara_label_scmap_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
SCINA_data <- prepare_data(TargetCell_data, label_SCINA_NoReject, c("label_SCINA_DiagPara", "DiagPara_label_SCINA_NoReject_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)
scPred_data <- prepare_data(TargetCell_data, label_scPred_NoReject_Annot, c("label_scPred_DiagPara_Annot", "DiagPara_label_scPred_NoReject_Annot_SVGLRglmnet_ROC"), Set_CellType = Set_Obs_CellType, Set_CellType_Reverse = Set_Ref_State)

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
y_max <- 717 # 示例，根据你的数据调整

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
  labs(title = paste0("Annotation on ", Set_Obs_CellType,"\n", Set_Title_End),
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
Name_Note <- paste0(Set_Method,"_",min_value_fileID,"_",Set_Obs_CellType,Set_Title_End)

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_ExportFolder <- folder_path # Name_ExportFolder <- paste0("Export_IntegAll")
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
# Name_ExportFolder <- paste0(Name_ExportFolder,"/",Name_time_wo_micro)
# if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
Name_Export <- paste0(Name_time_wo_micro,"_",Name_Note)


## Export PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export,"_Check_Bar.pdf"),
    width = 7, height = 9) #  width = 17, height = 17)

print(Plot_Bar)


dev.off()

## Export TSV
try(write_tsv(combined_data, paste0(Name_ExportFolder, "/", Name_Export, "_Check_Bar.tsv")))
