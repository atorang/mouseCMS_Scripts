# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load and filter data----------
load(paste0(PathToData,"/organoids.RData")) 
data$data$raw=data$data$raw[,order(data$col.annot$genotype)]
data$data$quantile=data$data$quantile[,order(data$col.annot$genotype)]
data$col.annot=data$col.annot[order(data$col.annot$genotype),]

exprs=data$data$quantile
m=apply(exprs,1,function(x){sum(x>=6)})
exprs=exprs[m>=2,]

s1=apply(exprs[,data$col.annot$batch=="1"],1,mean)
s2=apply(exprs[,data$col.annot$batch=="2"],1,mean)
exprs=exprs[abs(s2-s1)<1.5,]

iq=apply(exprs,1,IQR)
iq=sort(iq,decreasing = T)
exprs=exprs[names(iq)[1:1000],]

exprs.z=exprs
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
}

# heatmaps parameters-----------
library(ComplexHeatmap)
col_fun =circlize:: colorRamp2(seq(-1.5,1.5,3/6), colorRampPalette(c("#1483E4","white","#FF8C00"))(7) )
col.origin=c("steelblue3","salmon");names(col.origin)=c("CO","SIP")
data$col.annot$A="wt"
data$col.annot$A[setdiff(grep("A",data$col.annot$genotype),grep("A2",data$col.annot$genotype))]="mut"
data$col.annot$B="wt"
data$col.annot$B[grep("B",data$col.annot$genotype)]="mut"
data$col.annot$K="wt"
data$col.annot$K[grep("K",data$col.annot$genotype)]="mut"
data$col.annot$P="wt"
data$col.annot$P[grep("P",data$col.annot$genotype)]="mut"
data$col.annot$N="wt"
data$col.annot$N[grep("N",data$col.annot$genotype)]="mut"
data$col.annot$A2="wt"
data$col.annot$A2[grep("A2",data$col.annot$genotype)]="mut"
data$col.annot$S="wt"
data$col.annot$S[grep("S",data$col.annot$genotype)]="mut"
col.gene=c("gray","red");names(col.gene)=c("wt","mut")

# heatmap-----------
Heatmap(as.matrix(exprs.z),
        column_split =data$col.annot$genotype,
        cluster_column_slices = T,
        cluster_columns = T,
        column_dend_reorder = F,
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize=10),
        column_title = "%s",
        column_title_rot=90,
        row_km_repeats=100, 
        show_row_names = F,
        show_column_names = F,
        column_labels = data$col.annot$genotype,
        column_names_gp = gpar(fontsize=4),
        show_parent_dend_line = FALSE,
        row_title_rot = 0,
        col=col_fun,
        column_gap = unit(0.7, "mm"),
        border = F,
        top_annotation = HeatmapAnnotation(df=data$col.annot[,c(5,14:20)],
                                           col = list(origin=col.origin,
                                                      A=col.gene,B=col.gene,
                                                      K=col.gene,P=col.gene,
                                                      N=col.gene,A2=col.gene,
                                                      S=col.gene),
                                           simple_anno_size = unit(0.3, "cm"),
                                           show_legend = c(T,
                                                           #T,
                                                           T,F,F,F,F,F,F),
                                           gap = unit(0.1, "cm")),
       name="Z-Score", show_heatmap_legend = T,
       heatmap_legend_param = list(title="Z-Score",at=c(-2,0,2)),
       show_row_dend = F, show_column_dend = T)
