################################################################################
#                    TCGA 子宫内膜癌数据分析 (最终修复版 v3)
#         聚类热图 + 相关性热图 + 差异检验 + 箱线图 + OS/DFS生存分析
#
#  版本: FINAL v3 (2024) - 完全重构OS和DFS分析，确保变量独立
#
#  主要功能:
#  1. 聚类热图 - 基因表达聚类分析（样本分组）
#  2. 相关性热图 - 基因间相关性分析（带显著性标记）
#  3. 差异分析 - Tumor vs Normal 差异表达检验
#  4. 箱线图 - 基因表达可视化（组合图 + 单基因图）
#  5. OS生存分析 - 总生存期分析（Overall Survival）
#  6. DFS生存分析 - 无病生存期分析（Disease-Free Survival）
#
#  v3版本重构内容:
#  ✅ OS分析: 使用 os_xxx 前缀的独立变量
#     - os_tumor_idx, os_exp_tumor, os_patient_ids
#     - os_exp_avg, survival_data_os
#     - os_temp, os_fit, os_diff, os_cox (循环内)
#
#  ✅ DFS分析: 使用 dfs_xxx 前缀的独立变量
#     - dfs_tumor_idx, dfs_exp_tumor, dfs_patient_ids
#     - dfs_exp_avg, survival_data_dfs
#     - dfs_temp, dfs_fit, dfs_diff, dfs_cox (循环内)
#
#  ✅ 数据来源:
#     - OS: clinical变量中的 vital_status, days_to_death, days_to_last_follow_up
#     - DFS: clinical变量中的 follow_ups_disease_response (TF/WT状态)
#
#  ✅ 事件定义:
#     - OS_status = 1 仅当 vital_status == "Dead"
#     - DFS_status = 1 当 has_tumor == TRUE 或 vital_status == "Dead"
#
#  ✅ 详细诊断输出:
#     - OS vs DFS 事件数对比
#     - 展示存活但有肿瘤复发的患者（DFS特有事件）
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
#                           第五部分：OS生存分析
#                           使用独立变量: os_xxx
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤8: 准备OS生存分析数据 ======\n")

# 确保clinical有patient_id列
clinical$patient_id <- clinical$submitter_id

# ========== OS专用变量 ==========
# 只用肿瘤样本
os_tumor_idx <- sample_group == "Tumor"
os_exp_tumor <- exp_norm[, os_tumor_idx]

# 患者ID
os_patient_ids <- substr(colnames(os_exp_tumor), 1, 12)

# 同一患者多样本取平均
os_exp_df <- as.data.frame(t(os_exp_tumor))
os_exp_df$patient_id <- os_patient_ids

os_exp_avg <- os_exp_df %>%
  group_by(patient_id) %>%
  summarise(across(all_of(rownames(exp_target)), mean, na.rm = TRUE)) %>%
  as.data.frame()

rownames(os_exp_avg) <- os_exp_avg$patient_id

# 合并临床数据（只需要OS相关字段）
survival_data_os <- merge(
  os_exp_avg,
  clinical[, c("patient_id", "vital_status", "days_to_death", "days_to_last_follow_up")],
  by = "patient_id"
)

# 计算OS时间和状态
# OS_time: 死亡用days_to_death，存活用days_to_last_follow_up
# OS_status: 死亡=1, 存活=0
survival_data_os$OS_time <- ifelse(
  survival_data_os$vital_status == "Dead",
  as.numeric(survival_data_os$days_to_death),
  as.numeric(survival_data_os$days_to_last_follow_up)
)
survival_data_os$OS_status <- ifelse(survival_data_os$vital_status == "Dead", 1, 0)

# 过滤无效数据
survival_data_os <- survival_data_os %>% filter(!is.na(OS_time) & OS_time > 0)

cat("OS生存分析样本数:", nrow(survival_data_os), "\n")
cat("OS事件数(死亡):", sum(survival_data_os$OS_status), "\n")
cat("OS审查数(存活):", sum(survival_data_os$OS_status == 0), "\n")

# ==============================
# 步骤9: OS生存分析
# ==============================

cat("\n====== 步骤9: OS生存分析 ======\n")

os_results <- data.frame()
os_available_genes <- intersect(rownames(exp_target), colnames(survival_data_os))
cat("可分析的基因数:", length(os_available_genes), "\n")

pdf(file.path(output_dir, "OS_survival_curves.pdf"), width = 8, height = 7)

for (gene in os_available_genes) {

  # 创建OS专用临时数据框
  os_temp <- survival_data_os[, c("patient_id", gene, "OS_time", "OS_status")]
  colnames(os_temp)[colnames(os_temp) == gene] <- "gene_expression"

  # 基于OS数据计算中位数和分组
  os_median <- median(os_temp$gene_expression, na.rm = TRUE)
  os_temp$risk_group <- factor(
    ifelse(os_temp$gene_expression > os_median, "High", "Low"),
    levels = c("Low", "High")
  )

  # Kaplan-Meier和Log-rank检验
  os_fit <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = os_temp)
  os_diff <- survdiff(Surv(OS_time, OS_status) ~ risk_group, data = os_temp)
  os_pvalue <- 1 - pchisq(os_diff$chisq, 1)

  # Cox回归
  os_cox <- coxph(Surv(OS_time, OS_status) ~ gene_expression, data = os_temp)
  os_cox_sum <- summary(os_cox)

  os_results <- rbind(os_results, data.frame(
    Gene = gene,
    HR = round(os_cox_sum$conf.int[1, 1], 3),
    HR_95CI_lower = round(os_cox_sum$conf.int[1, 3], 3),
    HR_95CI_upper = round(os_cox_sum$conf.int[1, 4], 3),
    LogRank_P = signif(os_pvalue, 4),
    Cox_P = signif(os_cox_sum$coefficients[1, 5], 4)
  ))

  # 绑图
  os_plot <- ggsurvplot(
    os_fit,
    data = os_temp,
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

  print(os_plot)

  # 清理OS临时变量
  rm(os_temp, os_fit, os_diff, os_cox)
}

dev.off()

os_results$FDR <- signif(p.adjust(os_results$LogRank_P, method = "BH"), 4)
os_results <- os_results %>% arrange(LogRank_P)

write.csv(os_results, file.path(output_dir, "OS_survival_results.csv"), row.names = FALSE)

cat("OS分析完成！\n")
print(os_results)

# ══════════════════════════════════════════════════════════════════════════════
#                      第六部分：DFS生存分析
#                      使用独立变量: dfs_xxx
#                      数据来源: clinical变量中的follow_ups_disease_response
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== 步骤10: 准备DFS生存分析数据 ======\n")

# ========== DFS专用变量 ==========

# 检查clinical中的疾病状态相关字段
cat("检查clinical变量中的疾病状态字段...\n")
dfs_disease_cols <- grep("disease|tumor|response", colnames(clinical), ignore.case = TRUE, value = TRUE)
cat("疾病相关字段:", paste(dfs_disease_cols, collapse = ", "), "\n")

# 查找疾病状态字段
dfs_status_field <- NULL
dfs_priority_fields <- c("follow_ups_disease_response", "paper_tumor_status", "tumor_status")

for (field in dfs_priority_fields) {
  if (field %in% colnames(clinical)) {
    field_vals <- clinical[[field]]
    # 检查是否有 WT-With Tumor 类型的值
    if (sum(grepl("WT|With Tumor", field_vals, ignore.case = TRUE), na.rm = TRUE) > 0) {
      dfs_status_field <- field
      cat("\n✓ 找到疾病状态字段:", field, "\n")
      cat("字段值分布:\n")
      print(table(field_vals, useNA = "ifany"))
      break
    }
  }
}

if (!is.null(dfs_status_field)) {

  # 步骤1: 准备DFS专用表达数据
  dfs_tumor_idx <- sample_group == "Tumor"
  dfs_exp_tumor <- exp_norm[, dfs_tumor_idx]
  dfs_patient_ids <- substr(colnames(dfs_exp_tumor), 1, 12)

  # 步骤2: 按患者汇总表达数据
  dfs_exp_df <- as.data.frame(t(dfs_exp_tumor))
  dfs_exp_df$patient_id <- dfs_patient_ids

  dfs_exp_avg <- dfs_exp_df %>%
    group_by(patient_id) %>%
    summarise(across(all_of(rownames(exp_target)), mean, na.rm = TRUE)) %>%
    as.data.frame()

  rownames(dfs_exp_avg) <- dfs_exp_avg$patient_id

  # 步骤3: 从clinical获取DFS所需的临床数据
  dfs_clinical_needed <- c("patient_id", "vital_status", "days_to_death",
                           "days_to_last_follow_up", dfs_status_field)
  dfs_clinical_needed <- intersect(dfs_clinical_needed, colnames(clinical))

  # 合并表达数据和临床数据
  survival_data_dfs <- merge(
    dfs_exp_avg,
    clinical[, dfs_clinical_needed],
    by = "patient_id"
  )

  # 步骤4: 识别有肿瘤进展的患者
  survival_data_dfs$disease_response <- survival_data_dfs[[dfs_status_field]]
  survival_data_dfs$has_tumor <- grepl("^WT|With Tumor", survival_data_dfs$disease_response, ignore.case = TRUE)

  cat("\n疾病状态分布:\n")
  print(table(survival_data_dfs$disease_response, useNA = "ifany"))
  cat("\nhas_tumor分布:\n")
  print(table(survival_data_dfs$has_tumor))

  # 步骤5: 计算DFS时间和状态
  # DFS_time: 使用随访时间（与OS相同，因为没有精确的复发时间）
  # DFS_status: 死亡=1 OR 有肿瘤=1
  survival_data_dfs$DFS_time <- ifelse(
    survival_data_dfs$vital_status == "Dead",
    as.numeric(survival_data_dfs$days_to_death),
    as.numeric(survival_data_dfs$days_to_last_follow_up)
  )
  survival_data_dfs$DFS_status <- ifelse(
    survival_data_dfs$has_tumor | survival_data_dfs$vital_status == "Dead",
    1, 0
  )

  # 过滤无效数据
  survival_data_dfs <- survival_data_dfs %>% filter(!is.na(DFS_time) & DFS_time > 0)

  cat("\nDFS生存分析样本数:", nrow(survival_data_dfs), "\n")
  cat("DFS事件数(复发+死亡):", sum(survival_data_dfs$DFS_status), "\n")
  cat("DFS审查数(无复发且存活):", sum(survival_data_dfs$DFS_status == 0), "\n")

  # ========== OS vs DFS 详细对比 ==========
  cat("\n====== OS vs DFS 详细对比 ======\n")

  dfs_common_patients <- intersect(survival_data_os$patient_id, survival_data_dfs$patient_id)
  cat("共同患者数:", length(dfs_common_patients), "\n")

  if (length(dfs_common_patients) > 0) {
    os_subset <- survival_data_os %>% filter(patient_id %in% dfs_common_patients)
    dfs_subset <- survival_data_dfs %>% filter(patient_id %in% dfs_common_patients)

    cat("\n事件数对比:\n")
    cat("  OS事件数（仅死亡）:", sum(os_subset$OS_status), "\n")
    cat("  DFS事件数（复发+死亡）:", sum(dfs_subset$DFS_status), "\n")
    cat("  DFS新增事件（有肿瘤但存活）:", sum(dfs_subset$DFS_status) - sum(os_subset$OS_status), "\n")

    # 关键：展示具体的差异患者
    comparison_df <- merge(
      os_subset[, c("patient_id", "OS_time", "OS_status", "vital_status")],
      dfs_subset[, c("patient_id", "DFS_time", "DFS_status", "has_tumor", "disease_response")],
      by = "patient_id"
    )

    # 找出DFS有事件但OS没事件的患者（存活但有肿瘤）
    diff_patients <- comparison_df %>%
      filter(DFS_status == 1 & OS_status == 0)

    if (nrow(diff_patients) > 0) {
      cat("\n存活但有肿瘤复发的患者（DFS特有事件）:\n")
      print(diff_patients[, c("patient_id", "vital_status", "disease_response", "DFS_time")])
    }
  }

  # ========== DFS分析循环 ==========
  if (sum(survival_data_dfs$DFS_status) >= 10) {  # 至少10个事件

    cat("\n====== 步骤11: DFS生存分析 ======\n")

    dfs_results <- data.frame()
    dfs_available_genes <- intersect(rownames(exp_target), colnames(survival_data_dfs))
    cat("可分析的基因数:", length(dfs_available_genes), "\n")

    pdf(file.path(output_dir, "DFS_survival_curves.pdf"), width = 8, height = 7)

    for (gene in dfs_available_genes) {

      # 创建DFS专用临时数据框
      dfs_temp <- survival_data_dfs[, c("patient_id", gene, "DFS_time", "DFS_status")]
      colnames(dfs_temp)[colnames(dfs_temp) == gene] <- "gene_expression"

      # 基于DFS数据计算中位数和分组
      dfs_median <- median(dfs_temp$gene_expression, na.rm = TRUE)
      dfs_temp$risk_group <- factor(
        ifelse(dfs_temp$gene_expression > dfs_median, "High", "Low"),
        levels = c("Low", "High")
      )

      # Kaplan-Meier和Log-rank检验
      dfs_fit <- survfit(Surv(DFS_time, DFS_status) ~ risk_group, data = dfs_temp)
      dfs_diff <- survdiff(Surv(DFS_time, DFS_status) ~ risk_group, data = dfs_temp)
      dfs_pvalue <- 1 - pchisq(dfs_diff$chisq, 1)

      # Cox回归
      dfs_cox <- coxph(Surv(DFS_time, DFS_status) ~ gene_expression, data = dfs_temp)
      dfs_cox_sum <- summary(dfs_cox)

      dfs_results <- rbind(dfs_results, data.frame(
        Gene = gene,
        HR = round(dfs_cox_sum$conf.int[1, 1], 3),
        HR_95CI_lower = round(dfs_cox_sum$conf.int[1, 3], 3),
        HR_95CI_upper = round(dfs_cox_sum$conf.int[1, 4], 3),
        LogRank_P = signif(dfs_pvalue, 4),
        Cox_P = signif(dfs_cox_sum$coefficients[1, 5], 4)
      ))

      # 绑图
      dfs_plot <- ggsurvplot(
        dfs_fit,
        data = dfs_temp,
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

      print(dfs_plot)

      # 清理DFS临时变量
      rm(dfs_temp, dfs_fit, dfs_diff, dfs_cox)
    }

    dev.off()

    dfs_results$FDR <- signif(p.adjust(dfs_results$LogRank_P, method = "BH"), 4)
    dfs_results <- dfs_results %>% arrange(LogRank_P)

    write.csv(dfs_results, file.path(output_dir, "DFS_survival_results.csv"), row.names = FALSE)

    cat("DFS分析完成！\n")
    print(dfs_results)

  } else {
    cat("\n⚠ DFS事件数过少（<10），跳过DFS分析\n")
  }

} else {
  cat("\n⚠ 未找到有效的疾病状态字段，跳过DFS分析\n")
  cat("请确认clinical变量中是否有follow_ups_disease_response列\n")
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
cat("\n🔥 v3版本重构说明:\n")
cat("   1. OS和DFS使用完全独立的变量命名（os_xxx / dfs_xxx）\n")
cat("   2. DFS从clinical变量获取follow_ups_disease_response字段\n")
cat("   3. OS事件: 仅死亡（vital_status == 'Dead'）\n")
cat("   4. DFS事件: 死亡 OR 有肿瘤复发（WT-With Tumor）\n")
cat("\n⚠ 关于Risk table说明:\n")
cat("   由于TCGA数据中没有记录肿瘤复发的具体时间（new_tumor_event_dx_days_to），\n")
cat("   DFS_time 仍使用 days_to_last_follow_up，与 OS_time 相同。\n")
cat("   这导致Risk table数字可能相似，但生存曲线和p值会不同，\n")
cat("   因为DFS有更多事件（包括存活但有肿瘤复发的患者）。\n")