# 读取数据
otu <- read.table("otu.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 提取并删除Family列
Family <- otu$Family
otu_numeric <- otu[, -1]   # 删除第二列

# 按科合并（相同科的数值相加）
otu_merged <- rowsum(otu_numeric, group = Family)

# 计算相对丰度（每列除以该列总和）
otu_rel <- otu_merged / rep(colSums(otu_merged), each = nrow(otu_merged))

# 验证：每列总和应为1
colSums(otu_rel)

# 输出结果
write.table(otu_rel, file = "otu2.txt", sep = "\t", quote = FALSE, col.names = NA)



# 基于之前跑出来的 otu_rel（行=科，列=样品，值为相对丰度）
# 计算每个科在所有样品中的总平均相对丰度
otu_rel$Overall_mean <- rowMeans(otu_rel)

# 按总平均相对丰度降序排序
otu_sorted <- otu_rel[order(otu_rel$Overall_mean, decreasing = TRUE), ]

# 提取前10个科的名称
top10_families <- rownames(otu_sorted)[1:10]

# 创建新数据框：前10的科保留，其余合并为 "Others"
otu_top10 <- otu_sorted[, !colnames(otu_sorted) %in% "Overall_mean"]  # 删除Overall_mean列

# 分离前10和非前10
otu_top <- otu_top10[rownames(otu_top10) %in% top10_families, ]
otu_others <- otu_top10[!rownames(otu_top10) %in% top10_families, ]

# 将非前10的科合并为一行 "Others"
others_row <- colSums(otu_others)
otu_final <- rbind(otu_top, Others = others_row)

# 按总平均值重新排序（可选，让Others在最后）
otu_final$Overall_mean <- rowMeans(otu_final)
otu_final <- otu_final[order(otu_final$Overall_mean, decreasing = TRUE), ]
otu_final <- otu_final[, !colnames(otu_final) %in% "Overall_mean"]

# 输出结果
write.table(otu_final, file = "otu3.txt", sep = "\t", quote = FALSE, col.names = NA)