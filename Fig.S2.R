# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"

# load GSE146476 data----------
exprs=read.table(paste0(PathToData,"GSE146476_Intestinal-region-specific-Wnt-signalling_gene_counts.txt"))
rownames(exprs)=exprs$V1
colnames(exprs)=exprs[1,]
exprs=exprs[-1,]
exprs=exprs[,c(10:12,22:24)]
genes=rownames(exprs)
exprs=apply(exprs,2,as.numeric)
rownames(exprs)=genes
exprs.WT=data.frame(exprs)

# DE CO vs SI--------
class=c("CO","CO","CO","SIP","SIP","SIP")
col.annot=data.frame(class)
rownames(col.annot)=colnames(exprs.WT)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs.WT,
                              colData = col.annot,
                              design= ~ class)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("class","CO","SIP"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks=sort(ranks,decreasing = T)

res$padj[is.na(res$padj)]=1
de=res[res$padj<0.0001,]
de=de[abs(de$log2FoldChange)>2,]


# heatmap--------
exprs.z=log2(exprs.WT[rownames(de),]+1)
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=( exprs.z[i,]-mean(as.numeric( exprs.z[i,])))/sd( exprs.z[i,])
}
breaks=seq(-1.5,1.5,0.01)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pdf("Fig.S2A.DE.Heatmap.pdf",height = 5,width = 3)
pheatmap::pheatmap(exprs.z,breaks =breaks,color = col,show_rownames = F)
dev.off()

# load organoids data-----------
load(paste0(PathToData,"/organoids.RData")) 
exprs=data$data$quantile
m=apply(exprs,1,max)
exprs=exprs[m>5,]
genes=intersect(rownames(exprs.z),rownames(exprs))
exprs=exprs[genes,]


# PCA------------
library(scales)
pc=prcomp(t(exprs))
pch=rep(19,ncol(exprs))
pch[data$col.annot$origin=="SIP"]=17
pdf("Fig.S2B.PCA.pdf",4,4)
plot(pc$x[,1:2],pch=pch,col=alpha(data$col.annot$color,0.4))
legend("topleft",c("CO","SIP"),pch=unique(pch))
dev.off()

