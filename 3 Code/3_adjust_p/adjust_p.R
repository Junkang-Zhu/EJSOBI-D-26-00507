# 读取数据（第一列为SampleID，第二列为p值）
p_data <- read.table("p.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 查看数据结构
colnames(p_data)
head(p_data)

# 提取p值列（假设第二列是p值）
p_values <- p_data[, 2]

# 对p值进行FDR校正（Benjamini-Hochberg方法）
p_adjusted <- p.adjust(p_values, method = "BH")

# 将校正后的p值添加到数据框中
p_data$p_adjusted <- p_adjusted

# 导出结果
write.table(p_data, file = "p_adjusted.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 汇报校正结果
sig_before <- sum(p_values < 0.05, na.rm = TRUE)
sig_after <- sum(p_adjusted < 0.05, na.rm = TRUE)
cat(sprintf("校正前显著数: %d, 校正后显著数: %d\n", sig_before, sig_after))
