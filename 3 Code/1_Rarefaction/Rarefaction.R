rm(list=ls())#clear Global Environment
library(openxlsx)
library(vegan)
asv=read.xlsx("1.xlsx",1,rowNames = TRUE)

#求和查看每个样本的和
colSums(asv)
#抽平并把结果赋值给“otu_Flattening”
otu_Flattening =as.data.frame(t(rrarefy(t(asv),min(colSums(asv)))))
#查看抽平后的每个样本的和
colSums(otu_Flattening)
#将抽平后的otu结果表（otu_Flattening）保存到该工作目录下，并命名为“otu_Results.xlsx”
write.table(otu_Flattening, file="2.csv",sep =",", quote=FALSE)

