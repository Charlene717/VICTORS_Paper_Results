# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)

##### Set parameter ######
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))



################################################################################

# 获取R版本信息
r_version_info <- R.Version()

# 将R版本信息转换为数据框
r_version_df <- as.data.frame(as.list(r_version_info))

# 获取所有已安装的包信息
installed_packages <- installed.packages()

# 将已安装的包信息转换为数据框
installed_packages_df <- as.data.frame(installed_packages)

# 合并R版本信息和包信息
# 这里只是示例，你可以根据需要决定是否需要合并
# combined_info <- list(R_Version = r_version_df, Installed_Packages = installed_packages_df)

# 保存R版本信息为CSV文件
write.csv(r_version_df, paste0(Name_FileID,"_R_Version_Info.csv"), row.names = FALSE)

# 保存已安装的包信息为CSV文件
write.csv(installed_packages_df, paste0(Name_FileID,"_Installed_Packages_Info.csv"), row.names = FALSE)

# 打印保存信息
cat("R version information saved to R_Version_Info.csv\n")
cat("Installed packages information saved to Installed_Packages_Info.csv\n")


################################################################################

# 获取当前R会话信息
session_info <- sessionInfo()

# 提取R版本信息
r_version_info <- data.frame(
  R_version = ifelse(length(session_info$R.version$version.string) == 0, NA, session_info$R.version$version.string),
  Platform = ifelse(length(session_info$R.version$platform) == 0, NA, session_info$R.version$platform),
  Arch = ifelse(length(session_info$R.version$arch) == 0, NA, session_info$R.version$arch),
  OS = ifelse(length(session_info$R.version$os) == 0, NA, session_info$R.version$os),
  System = ifelse(length(session_info$R.version$system) == 0, NA, session_info$R.version$system),
  Date = ifelse(length(session_info$R.version$date) == 0, NA, session_info$R.version$date),
  stringsAsFactors = FALSE
)

# 提取加载的包信息
loaded_packages <- session_info$otherPkgs
loaded_packages_info <- do.call(rbind, lapply(loaded_packages, function(pkg) {
  data.frame(
    Package = pkg$Package,
    Version = pkg$Version,
    Priority = ifelse(is.null(pkg$Priority), NA, pkg$Priority),
    Built = pkg$Built,
    stringsAsFactors = FALSE
  )
}))

# 提取附加的包信息
attached_packages <- session_info$loadedOnly
attached_packages_info <- do.call(rbind, lapply(attached_packages, function(pkg) {
  data.frame(
    Package = pkg$Package,
    Version = pkg$Version,
    Priority = ifelse(is.null(pkg$Priority), NA, pkg$Priority),
    Built = pkg$Built,
    stringsAsFactors = FALSE
  )
}))

# 保存R版本信息为CSV文件
write.csv(r_version_info, paste0(Name_FileID,"_sessionInfo_R_Version_Info.csv"), row.names = FALSE)

# 保存加载的包信息为CSV文件
write.csv(loaded_packages_info, paste0(Name_FileID,"_sessionInfo_Loaded_Packages_Info.csv"), row.names = FALSE)

# 保存附加的包信息为CSV文件
write.csv(attached_packages_info, paste0(Name_FileID,"_sessionInfo_Attached_Packages_Info.csv"), row.names = FALSE)

# 打印保存信息
cat("R version information saved to R_Version_Info.csv\n")
cat("Loaded packages information saved to Loaded_Packages_Info.csv\n")
cat("Attached packages information saved to Attached_Packages_Info.csv\n")
