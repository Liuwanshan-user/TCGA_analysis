################################################################################
#                    TCGA 子宫内膜癌数据分析 (最终修复版)
#         聚类热图 + 相关性热图 + 差异检验 + 箱线图 + OS/DFS生存分析
#
#  版本: FINAL (2024)
#  
#  主要功能:
#  1. 聚类热图 - 基因表达聚类分析（样本分组）
#  2. 相关性热图 - 基因间相关性分析（带显著性标记）
#  3. 差异分析 - Tumor vs Normal 差异表达检验
#  4. 箱线图 - 基因表达可视化（组合图 + 单基因图）
#  5. OS生存分析 - 总生存期分析
#  6. DFS生存分析 - 无病生存期分析
#
#  核心修复内容:
#  ✅ 修复1: OS和DFS生存曲线Risk table显示错误的问题
#     - 每个基因使用独立的临时数据框（temp_data_os / temp_data_dfs）
#     - 统一基因列名为gene_expression，避免动态引用冲突
#     - 直接使用公式而非as.formula()，避免环境问题
#  
#  ✅ 修复2: DFS分析完全独立计算，不再依赖OS数据
#     - DFS从头开始独立准备数据（表达数据 + 疾病状态 + 临床数据）
#     - DFS使用独立的患者集（可能小于OS，因为需要疾病状态数据）
#     - DFS独立计算中位数和分组，不受OS影响
#
#  ✅ 修复3: 每次循环后清理临时变量，避免内存泄漏和数据污染
#
#  使用说明:
#  1. 准备 Gene list.csv 文件（包含Gene列）
#  2. 确保网络连接（需要下载TCGA数据）
#  3. 运行脚本，结果将保存在 TCGA_UCEC_results/ 目录
#
################################################################################

# ==============================
# 安装和加载必要的包
# ==============================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_bioc <- c("TCGAbiolinks", "SummarizedExperiment", "edgeR")
for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

packages_cran <- c("survival", "survminer", "ggplot2", "dplyr", "tidyr",
                   "pheatmap", "ggpubr", "tibble", "readr", "Hmisc")
for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(tibble)
library(readr)
library(Hmisc)

# ==============================
# 设置工作环境
# ==============================
set.seed(1)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 读取基因列表
gene_df <- read_csv("Gene list.csv")
gene_list <- gene_df$Gene

# TCGA项目设置
project <- "TCGA-UCEC"  # 子宫内膜癌
cancer_type <- "UCEC"

# 创建结果输出目录
output_dir <- paste0("TCGA_", cancer_type, "_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
#                           第一部分：数据预处理
# ══════════════════════════════════════════════════════════════════════════════

cat("====== 步骤1: 下载TCGA数据 ======\n")

# 下载表达数据
query_exp <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp)
data_exp <- GDCprepare(query_exp)

exp_matrix <- assay(data_exp)
gene_info <- rowData(data_exp)

# 下载临床数据
clinical <- GDCquery_clinic(project = project, type = "clinical")

cat("样本数量:", ncol(exp_matrix), "\n")

# ==============================
# 步骤2: 提取目标基因
# ==============================

cat("\n====== 步骤2: 提取目标基因 ======\n")

gene_mapping <- data.frame(
  ensembl_id = rownames(exp_matrix),
  gene_symbol = gene_info$gene_name
)

target_genes <- gene_mapping %>% filter(gene_symbol %in% gene_list)

cat("找到基因:", nrow(target_genes), "/", length(gene_list), "\n")

missing_genes <- setdiff(gene_list, target_genes$gene_symbol)
if (length(missing_genes) > 0) {
  cat("未找到:", paste(missing_genes, collapse = ", "), "\n")
}

exp_target <- exp_matrix[target_genes$ensembl_id, , drop = FALSE]
rownames(exp_target) <- target_genes$gene_symbol

# ==============================
# 步骤3: 样本分组与归一化
# ==============================

cat("\n====== 步骤3: 样本分组与归一化 ======\n")

sample_type <- substr(colnames(exp_target), 14, 15)
sample_group <- ifelse(as.numeric(sample_type) < 10, "Tumor", "Normal")

cat("肿瘤样本:", sum(sample_group == "Tumor"), "\n")
cat("正常样本:", sum(sample_group == "Normal"), "\n")

# 归一化 (log2 CPM)
dge <- DGEList(counts = exp_target)
dge <- calcNormFactors(dge)
exp_norm <- cpm(dge, log = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
#                           第二部分：聚类热图
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤4: 绑制聚类热图 ======\n")

# 样本注释
annotation_col <- data.frame(
  Group = sample_group,
  row.names = colnames(exp_norm)
)

annotation_colors <- list(
  Group = c(Tumor = "#E41A1C", Normal = "#377EB8")
)

# 聚类热图 (只保存png)
png(file.path(output_dir, "clustering_heatmap.png"), width = 1400, height = 800, res = 100)
pheatmap(exp_norm,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         scale = "row",
         clustering_method = "ward.D2",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_colnames = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = paste0("Gene Expression Clustering Heatmap (", cancer_type, ")"),
         fontsize_row = 10)
dev.off()

cat("聚类热图已保存\n")

# ══════════════════════════════════════════════════════════════════════════════
#                        第三部分：基因相关性聚类热图（带显著性）
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤5: 绑制相关性聚类热图 ======\n")

# 计算基因间相关性和p值
cor_result <- rcorr(t(exp_norm), type = "pearson")
gene_cor <- cor_result$r  # 相关系数矩阵
gene_pval <- cor_result$P  # p值矩阵

# 创建显著性标记矩阵
sig_matrix <- matrix("", nrow = nrow(gene_cor), ncol = ncol(gene_cor))
rownames(sig_matrix) <- rownames(gene_cor)
colnames(sig_matrix) <- colnames(gene_cor)

sig_matrix[gene_pval < 0.001] <- "***"
sig_matrix[gene_pval >= 0.001 & gene_pval < 0.01] <- "**"
sig_matrix[gene_pval >= 0.01 & gene_pval < 0.05] <- "*"

# 创建显示文本：相关系数 + 显著性标记
display_matrix <- matrix(paste0(sprintf("%.2f", gene_cor), sig_matrix),
                         nrow = nrow(gene_cor), ncol = ncol(gene_cor))
rownames(display_matrix) <- rownames(gene_cor)
colnames(display_matrix) <- colnames(gene_cor)

# 对角线设为空
diag(display_matrix) <- ""

# 相关性聚类热图 (只保存png)
png(file.path(output_dir, "correlation_heatmap.png"), width = 1000, height = 900, res = 100)
pheatmap(gene_cor,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         clustering_method = "ward.D2",
         display_numbers = display_matrix,
         fontsize_number = 7,
         fontsize_row = 10,
         fontsize_col = 10,
         main = paste0("Gene Correlation Heatmap (", cancer_type, ")\n* p<0.05, ** p<0.01, *** p<0.001"))
dev.off()

# 保存相关性结果表格
cor_df <- as.data.frame(as.table(gene_cor))
colnames(cor_df) <- c("Gene1", "Gene2", "Correlation")
pval_df <- as.data.frame(as.table(gene_pval))
colnames(pval_df) <- c("Gene1", "Gene2", "P_value")
cor_results <- merge(cor_df, pval_df, by = c("Gene1", "Gene2"))
cor_results <- cor_results %>%
  filter(as.character(Gene1) < as.character(Gene2)) %>%  # 去除重复和对角线
  mutate(
    Correlation = round(Correlation, 4),
    P_value = signif(P_value, 4),
    FDR = signif(p.adjust(P_value, method = "BH"), 4),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(P_value)

write.csv(cor_results, file.path(output_dir, "gene_correlation_results.csv"), row.names = FALSE)

cat("相关性热图已保存\n")
cat("显著性标记: * p<0.05, ** p<0.01, *** p<0.001\n")

# ══════════════════════════════════════════════════════════════════════════════
#                     第四部分：差异表达统计检验 + 箱线图
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤6: 差异表达统计检验 ======\n")

# 准备长格式数据
exp_long <- as.data.frame(t(exp_norm)) %>%
  mutate(Sample = rownames(.), Group = sample_group) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "Gene", values_to = "Expression")

# 对每个基因进行统计检验
deg_results <- data.frame()

for (gene in unique(exp_long$Gene)) {
  
  gene_data <- exp_long %>% filter(Gene == gene)
  tumor_exp <- gene_data$Expression[gene_data$Group == "Tumor"]
  normal_exp <- gene_data$Expression[gene_data$Group == "Normal"]
  
  # Wilcoxon检验
  wilcox_test <- wilcox.test(tumor_exp, normal_exp)
  
  # t检验
  t_test <- t.test(tumor_exp, normal_exp)
  
  # 计算统计量
  mean_tumor <- mean(tumor_exp, na.rm = TRUE)
  mean_normal <- mean(normal_exp, na.rm = TRUE)
  logFC <- mean_tumor - mean_normal
  
  deg_results <- rbind(deg_results, data.frame(
    Gene = gene,
    Mean_Tumor = round(mean_tumor, 3),
    Mean_Normal = round(mean_normal, 3),
    logFC = round(logFC, 3),
    FC = round(2^logFC, 3),
    Wilcoxon_P = signif(wilcox_test$p.value, 4),
    T_test_P = signif(t_test$p.value, 4)
  ))
}

# FDR校正
deg_results$Wilcoxon_FDR <- signif(p.adjust(deg_results$Wilcoxon_P, method = "BH"), 4)
deg_results$T_test_FDR <- signif(p.adjust(deg_results$T_test_P, method = "BH"), 4)

# 显著性标记
deg_results$Significance <- ifelse(
  deg_results$Wilcoxon_FDR < 0.05,
  ifelse(deg_results$logFC > 0, "Up", "Down"),
  "NS"
)

deg_results <- deg_results %>% arrange(Wilcoxon_P)

write.csv(deg_results, file.path(output_dir, "DEG_statistics.csv"), row.names = FALSE)

cat("差异统计完成！\n")
print(deg_results)

# ==============================
# 步骤7: 箱线图
# ==============================

cat("\n====== 步骤7: 绑制箱线图 ======\n")

exp_long$Group <- factor(exp_long$Group, levels = c("Normal", "Tumor"))

# 组合箱线图
p_all <- ggplot(exp_long, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7) +
  stat_compare_means(aes(group = Group), method = "wilcox.test", 
                     label = "p.signif", label.y.npc = 0.95) +
  scale_fill_manual(values = c("Normal" = "#377EB8", "Tumor" = "#E41A1C")) +
  labs(title = paste0("Gene Expression: Normal vs Tumor (", cancer_type, ")"),
       x = NULL, y = "Expression (log2 CPM)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "boxplot_all_genes.pdf"), p_all, 
       width = max(10, length(unique(exp_long$Gene)) * 0.8), height = 6)
ggsave(file.path(output_dir, "boxplot_all_genes.png"), p_all, 
       width = max(10, length(unique(exp_long$Gene)) * 0.8), height = 6, dpi = 300)

# 单基因箱线图
pdf(file.path(output_dir, "boxplot_individual_genes.pdf"), width = 5, height = 5)

for (gene in unique(exp_long$Gene)) {
  
  gene_data <- exp_long %>% filter(Gene == gene)
  
  p_single <- ggplot(gene_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    stat_compare_means(method = "wilcox.test", label.x = 1.5, label.y.npc = 0.95) +
    scale_fill_manual(values = c("Normal" = "#377EB8", "Tumor" = "#E41A1C")) +
    labs(title = gene, x = NULL, y = "Expression (log2 CPM)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none"
    )
  
  print(p_single)
}

dev.off()

cat("箱线图已保存\n")

# ══════════════════════════════════════════════════════════════════════════════
#                           第五部分：OS生存分析 (已修复)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤8: 准备OS生存分析数据 ======\n")

# 只用肿瘤样本
tumor_idx <- sample_group == "Tumor"
exp_tumor <- exp_norm[, tumor_idx]

# 患者ID
patient_id <- substr(colnames(exp_tumor), 1, 12)

# 同一患者多样本取平均
exp_df <- as.data.frame(t(exp_tumor))
exp_df$patient_id <- patient_id

exp_avg_os <- exp_df %>%
  group_by(patient_id) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()

rownames(exp_avg_os) <- exp_avg_os$patient_id
exp_avg_os$patient_id <- NULL

# 合并临床数据
clinical$patient_id <- clinical$submitter_id

survival_data_os <- merge(
  data.frame(patient_id = rownames(exp_avg_os), exp_avg_os, check.names = FALSE),
  clinical[, c("patient_id", "vital_status", "days_to_death", 
               "days_to_last_follow_up")],
  by = "patient_id"
)

# 计算OS
survival_data_os$OS_time <- ifelse(
  survival_data_os$vital_status == "Dead",
  as.numeric(survival_data_os$days_to_death),
  as.numeric(survival_data_os$days_to_last_follow_up)
)
survival_data_os$OS_status <- ifelse(survival_data_os$vital_status == "Dead", 1, 0)

survival_data_os <- survival_data_os %>% filter(!is.na(OS_time) & OS_time > 0)

cat("OS生存分析样本数:", nrow(survival_data_os), "\n")

# ==============================
# 步骤9: OS生存分析 (已修复Risk table问题)
# ==============================

cat("\n====== 步骤9: OS生存分析 ======\n")

os_results <- data.frame()
available_genes_os <- intersect(colnames(exp_avg_os), colnames(survival_data_os))

pdf(file.path(output_dir, "OS_survival_curves.pdf"), width = 8, height = 7)

for (gene in available_genes_os) {
  
  # 🔥 核心修复：创建完全独立的临时数据框
  temp_data_os <- survival_data_os[, c("patient_id", gene, "OS_time", "OS_status", "vital_status")]
  
  # 重命名基因列为gene_expression，避免列名冲突
  colnames(temp_data_os)[colnames(temp_data_os) == gene] <- "gene_expression"
  
  # 基于OS数据集计算中位数和分组
  median_exp_os <- median(temp_data_os$gene_expression, na.rm = TRUE)
  temp_data_os$exp_group <- factor(
    ifelse(temp_data_os$gene_expression > median_exp_os, "High", "Low"),
    levels = c("Low", "High")
  )
  
  # 🔥 关键：直接使用公式，明确指定数据来源
  fit_os <- survfit(Surv(OS_time, OS_status) ~ exp_group, data = temp_data_os)
  diff_os <- survdiff(Surv(OS_time, OS_status) ~ exp_group, data = temp_data_os)
  p_value_os <- 1 - pchisq(diff_os$chisq, 1)
  
  # Cox模型
  cox_model_os <- coxph(Surv(OS_time, OS_status) ~ gene_expression, data = temp_data_os)
  cox_sum_os <- summary(cox_model_os)
  
  os_results <- rbind(os_results, data.frame(
    Gene = gene,
    HR = round(cox_sum_os$conf.int[1, 1], 3),
    HR_95CI_lower = round(cox_sum_os$conf.int[1, 3], 3),
    HR_95CI_upper = round(cox_sum_os$conf.int[1, 4], 3),
    LogRank_P = signif(p_value_os, 4),
    Cox_P = signif(cox_sum_os$coefficients[1, 5], 4)
  ))
  
  # 🔥 关键：明确传入fit和data
  km_plot_os <- ggsurvplot(
    fit_os,
    data = temp_data_os,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = c("#2E9FDF", "#E7B800"),
    title = paste0("OS: ", gene, " (", cancer_type, ")"),
    xlab = "Time (days)",
    ylab = "Overall Survival",
    legend.title = gene,
    legend.labs = c("Low", "High"),
    ggtheme = theme_bw()
  )
  
  print(km_plot_os)
  
  # 清理临时变量
  rm(temp_data_os, fit_os, diff_os, cox_model_os)
}

dev.off()

os_results$FDR <- signif(p.adjust(os_results$LogRank_P, method = "BH"), 4)
os_results <- os_results %>% arrange(LogRank_P)

write.csv(os_results, file.path(output_dir, "OS_survival_results.csv"), row.names = FALSE)

cat("OS分析完成！\n")
print(os_results)

# ══════════════════════════════════════════════════════════════════════════════
#                      第六部分：DFS生存分析 (完全重新计算)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤10: 准备DFS生存分析数据 (独立计算) ======\n")

# 🔥 关键修复：DFS完全独立准备数据，不依赖OS

# 从 colData 获取样本信息
sample_info <- as.data.frame(colData(data_exp))

# 查找疾病状态相关字段
cat("查找疾病状态字段...\n")

all_cols <- colnames(sample_info)
disease_cols <- grep("disease|tumor|response", all_cols, ignore.case = TRUE, value = TRUE)
cat("找到的相关字段:", paste(disease_cols, collapse = ", "), "\n")

# 优先查找的字段
priority_fields <- c(
  "paper_tumor_status",
  "follow_ups_disease_response",
  "tumor_status",
  "disease_response"
)

# 找到第一个存在且有"WITH TUMOR"或"TUMOR"值的字段
status_field <- NULL
for (field in c(priority_fields, disease_cols)) {
  if (field %in% colnames(sample_info)) {
    vals <- sample_info[[field]]
    has_tumor <- sum(grepl("TUMOR|tumor|WITH|with", vals, ignore.case = TRUE), na.rm = TRUE)
    if (has_tumor > 0) {
      status_field <- field
      cat("\n选择字段:", field, "\n")
      cat("字段值分布:\n")
      print(table(vals, useNA = "ifany"))
      break
    }
  }
}

if (!is.null(status_field)) {
  
  # 🔥 步骤1: 独立准备DFS的表达数据
  tumor_idx_dfs <- sample_group == "Tumor"
  exp_tumor_dfs <- exp_norm[, tumor_idx_dfs]
  
  # 🔥 步骤2: 提取样本的疾病状态信息
  tumor_info <- sample_info[tumor_idx_dfs, ]
  tumor_info$sample_id <- rownames(tumor_info)
  tumor_info$patient_id <- substr(tumor_info$sample_id, 1, 12)
  tumor_info$disease_status <- tumor_info[[status_field]]
  
  # 识别有肿瘤/进展的样本
  tumor_info$has_tumor <- grepl("^WT|With Tumor$", tumor_info$disease_status)
  
  cat("\n肿瘤样本中疾病状态分布:\n")
  print(table(tumor_info$disease_status, useNA = "ifany"))
  
  # 🔥 步骤3: 对每个患者，合并表达数据和疾病状态
  # 如果同一患者有多个样本，取平均表达值，并优先选择有进展的状态
  
  # 合并表达数据与疾病状态
  exp_dfs_df <- as.data.frame(t(exp_tumor_dfs))
  exp_dfs_df$sample_id <- rownames(exp_dfs_df)
  exp_dfs_df <- merge(exp_dfs_df, 
                      tumor_info[, c("sample_id", "patient_id", "disease_status", "has_tumor")],
                      by = "sample_id")
  
  # 按患者ID汇总
  patient_dfs_data <- exp_dfs_df %>%
    group_by(patient_id) %>%
    summarise(
      across(all_of(rownames(exp_target)), mean, na.rm = TRUE),  # 表达量取平均
      has_tumor = max(has_tumor, na.rm = TRUE),  # 只要有一个样本有肿瘤就算有
      disease_status = first(disease_status[has_tumor == max(has_tumor)])  # 优先取有肿瘤的状态
    ) %>%
    as.data.frame()
  
  rownames(patient_dfs_data) <- patient_dfs_data$patient_id
  
  cat("\n按患者汇总后的数据:\n")
  cat("患者数:", nrow(patient_dfs_data), "\n")
  cat("has_tumor分布:\n")
  print(table(patient_dfs_data$has_tumor, useNA = "ifany"))
  
  # 🔥 步骤4: 独立合并临床数据
  survival_data_dfs <- merge(
    patient_dfs_data,
    clinical[, c("patient_id", "vital_status", "days_to_death", "days_to_last_follow_up")],
    by = "patient_id"
  )
  
  # 🔥 步骤5: 计算DFS时间和状态
  survival_data_dfs$DFS_time <- ifelse(
    survival_data_dfs$vital_status == "Dead",
    as.numeric(survival_data_dfs$days_to_death),
    as.numeric(survival_data_dfs$days_to_last_follow_up)
  )
  
  # DFS状态：有进展/复发 OR 死亡 = 1
  survival_data_dfs$DFS_status <- ifelse(
    survival_data_dfs$has_tumor | survival_data_dfs$vital_status == "Dead",
    1, 0
  )
  
  # 过滤无效数据
  survival_data_dfs <- survival_data_dfs %>% filter(!is.na(DFS_time) & DFS_time > 0)
  
  cat("\nDFS生存分析样本数:", nrow(survival_data_dfs), "\n")
  
  # 详细统计对比
  cat("\n====== DFS vs OS 详细对比 ======\n")
  cat("OS样本数:", nrow(survival_data_os), "\n")
  cat("DFS样本数:", nrow(survival_data_dfs), "\n")
  
  # 找出共同患者进行对比
  common_patients <- intersect(survival_data_os$patient_id, survival_data_dfs$patient_id)
  cat("共同患者数:", length(common_patients), "\n")
  
  if (length(common_patients) > 0) {
    os_common <- survival_data_os %>% filter(patient_id %in% common_patients)
    dfs_common <- survival_data_dfs %>% filter(patient_id %in% common_patients)
    
    cat("\n共同患者中:\n")
    cat("  OS事件数（仅死亡）:", sum(os_common$OS_status), "\n")
    cat("  DFS事件数（进展+死亡）:", sum(dfs_common$DFS_status), "\n")
    cat("  新增DFS事件（有肿瘤但未死）:", 
        sum(dfs_common$DFS_status) - sum(os_common$OS_status), "\n")
  }
  
  # 检查是否有差异
  if (sum(survival_data_dfs$DFS_status) > length(common_patients) * 0.1) {  # 至少有10%的事件
    
    cat("\n开始DFS分析...\n")
    
    dfs_results <- data.frame()
    available_genes_dfs <- intersect(rownames(exp_target), colnames(survival_data_dfs))
    
    pdf(file.path(output_dir, "DFS_survival_curves.pdf"), width = 8, height = 7)
    
    for (gene in available_genes_dfs) {
      if (gene %in% colnames(survival_data_dfs)) {
        
        # 🔥 核心修复：创建完全独立的临时数据框
        temp_data_dfs <- survival_data_dfs[, c("patient_id", gene, "DFS_time", "DFS_status", 
                                               "has_tumor", "disease_status", "vital_status")]
        
        # 重命名基因列为gene_expression，避免列名冲突
        colnames(temp_data_dfs)[colnames(temp_data_dfs) == gene] <- "gene_expression"
        
        # 基于DFS数据集独立计算中位数和分组
        median_exp_dfs <- median(temp_data_dfs$gene_expression, na.rm = TRUE)
        temp_data_dfs$exp_group <- factor(
          ifelse(temp_data_dfs$gene_expression > median_exp_dfs, "High", "Low"),
          levels = c("Low", "High")
        )
        
        # 🔥 关键：直接使用公式，明确指定数据来源
        fit_dfs <- survfit(Surv(DFS_time, DFS_status) ~ exp_group, data = temp_data_dfs)
        diff_dfs <- survdiff(Surv(DFS_time, DFS_status) ~ exp_group, data = temp_data_dfs)
        p_value_dfs <- 1 - pchisq(diff_dfs$chisq, 1)
        
        # Cox模型
        cox_model_dfs <- coxph(Surv(DFS_time, DFS_status) ~ gene_expression, data = temp_data_dfs)
        cox_sum_dfs <- summary(cox_model_dfs)
        
        dfs_results <- rbind(dfs_results, data.frame(
          Gene = gene,
          HR = round(cox_sum_dfs$conf.int[1, 1], 3),
          HR_95CI_lower = round(cox_sum_dfs$conf.int[1, 3], 3),
          HR_95CI_upper = round(cox_sum_dfs$conf.int[1, 4], 3),
          LogRank_P = signif(p_value_dfs, 4),
          Cox_P = signif(cox_sum_dfs$coefficients[1, 5], 4)
        ))
        
        # 🔥 关键：明确传入fit和data，使用完全独立的变量
        km_plot_dfs <- ggsurvplot(
          fit_dfs,
          data = temp_data_dfs,
          pval = TRUE,
          pval.method = TRUE,
          conf.int = TRUE,
          risk.table = TRUE,
          risk.table.height = 0.25,
          palette = c("#2E9FDF", "#E7B800"),
          title = paste0("DFS: ", gene, " (", cancer_type, ")"),
          xlab = "Time (days)",
          ylab = "Disease-Free Survival",
          legend.title = gene,
          legend.labs = c("Low", "High"),
          ggtheme = theme_bw()
        )
        
        print(km_plot_dfs)
        
        # 清理临时变量
        rm(temp_data_dfs, fit_dfs, diff_dfs, cox_model_dfs)
      }
    }
    
    dev.off()
    
    dfs_results$FDR <- signif(p.adjust(dfs_results$LogRank_P, method = "BH"), 4)
    dfs_results <- dfs_results %>% arrange(LogRank_P)
    
    write.csv(dfs_results, file.path(output_dir, "DFS_survival_results.csv"), row.names = FALSE)
    
    cat("\nDFS分析完成！\n")
    print(dfs_results)
    
  } else {
    cat("\n注意: DFS事件数过少，跳过DFS分析\n")
  }
  
} else {
  cat("\n未找到有效的疾病状态字段，跳过DFS分析\n")
}

# ══════════════════════════════════════════════════════════════════════════════
#                              分析完成
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("══════════════════════════════════════════\n")
cat("              分析完成！\n")
cat("══════════════════════════════════════════\n")
cat("\n结果目录:", output_dir, "\n")

cat("\n输出文件:\n")
cat("  [热图]\n")
cat("    - clustering_heatmap.png        聚类热图\n")
cat("    - correlation_heatmap.png       相关性聚类热图（带显著性）\n")
cat("    - gene_correlation_results.csv  相关性统计结果\n")
cat("  [差异分析]\n")
cat("    - DEG_statistics.csv            差异统计结果\n")
cat("    - boxplot_all_genes.pdf/png     组合箱线图\n")
cat("    - boxplot_individual_genes.pdf  单基因箱线图\n")
cat("  [生存分析]\n")
cat("    - OS_survival_results.csv       OS生存分析结果\n")
cat("    - OS_survival_curves.pdf        OS生存曲线 (已修复Risk table)\n")
if (exists("dfs_results") && nrow(dfs_results) > 0) {
  cat("    - DFS_survival_results.csv      DFS生存分析结果\n")
  cat("    - DFS_survival_curves.pdf       DFS生存曲线 (已修复Risk table)\n")
}

cat("\n====== 结果摘要 ======\n")

cat("\n差异表达 (Wilcoxon FDR < 0.05):\n")
sig_deg <- deg_results %>% filter(Wilcoxon_FDR < 0.05)
if (nrow(sig_deg) > 0) {
  print(sig_deg[, c("Gene", "logFC", "FC", "Wilcoxon_P", "Wilcoxon_FDR", "Significance")])
} else {
  cat("  无显著差异基因\n")
}

cat("\nOS生存 (LogRank P < 0.05):\n")
sig_os <- os_results %>% filter(LogRank_P < 0.05)
if (nrow(sig_os) > 0) {
  print(sig_os)
} else {
  cat("  无显著基因\n")
}

if (exists("dfs_results") && nrow(dfs_results) > 0) {
  cat("\nDFS生存 (LogRank P < 0.05):\n")
  sig_dfs <- dfs_results %>% filter(LogRank_P < 0.05)
  if (nrow(sig_dfs) > 0) {
    print(sig_dfs)
  } else {
    cat("  无显著基因\n")
  }
}

cat("\n分析完毕！\n")
cat("\n🔥 已完全修复:\n")
cat("   1. OS和DFS生存曲线Risk table显示错误\n")
cat("   2. DFS分析现在完全独立计算，不再依赖OS数据\n")
cat("   3. DFS使用独立的患者集和表达数据分组\n")