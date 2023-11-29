# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# define function-----------
iClassification.class.data.to.subset.cols <- function (data, col.name, values){
  .check.data.init(data, "iClassification.class.data.to.subset.cols")
  data.subset = data$col.annot[,col.name] %in% values
  if(sum(data.subset)==1){ #cast data values as matrices in case number of cols = 1 (otherwise the matrices become vectors)
    if(nrow(data$row.annot)==1){
      stop(paste("[iClassification.class.data.to.subset.cols]: can't subset to a single column if data has only 1 row!", sep=""), call. = TRUE)	
    }
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = matrix(data$data[[i]][,data.subset], ncol=1)				
    }		
    # subset the annotation		
    data$col.annot = data$col.annot[data.subset,]				
  }else{	
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = data$data[[i]][,data.subset]				
    }		
    # subset the annotation		
    data$col.annot = data$col.annot[data.subset,]				
  }
  
  return(data)
}

# CMS colors---------
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74"); names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# load data-----------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "genotype", setdiff(data$col.annot$genotype,c("AK","AKPN"))) 

organoids= setdiff(data$col.annot$genotype,c("A"))
# load cell type signatures by Haber complete list-----------
sig=data.frame(readxl::read_xlsx(paste0(PathToData,"/Haber.Combine.All.xlsx")))
for (i in 1:ncol(sig)) {
  x=unique(sig[,i])
  sig[,i]=NA
  sig[1:length(x),i]=x
  
}
sig1=unlist(sig)

sig2=data.frame(gene=sig1,celltype=c(rep("Enteroendocrine",nrow(sig)),
                                     rep("Enterocyte.Immature.Distal",nrow(sig)),
                                     rep("Enterocyte.Immature.Proximal",nrow(sig)),
                                     rep("Enterocyte.Mature.Distal",nrow(sig)),
                                     rep("Enterocyte.Mature.Proximal",nrow(sig)),
                                     rep("Enterocyte",nrow(sig)),
                                     rep("Enterocyte.Progenitor.Early",nrow(sig)),
                                     rep("Enterocyte.Progenitor.Late",nrow(sig)),
                                     rep("Goblet",nrow(sig)),
                                     rep("Paneth",nrow(sig)),
                                     rep("Stem",nrow(sig)),
                                     rep("TA.G2",nrow(sig)),
                                     rep("Tuft",nrow(sig))
))
sig2=sig2[!is.na(sig2$gene),]
sig2$celltype[duplicated(sig2$gene)]
sig2$gene[duplicated(sig2$gene)]=NA
sig2=sig2[!is.na(sig2$gene),]
rownames(sig2)=sig2$gene
table(sig2$celltype)


sig.complete=sig2[-which(sig2$celltype%in%c("Enterocyte.Immature.Distal","Enterocyte.Immature.Proximal","Enterocyte.Progenitor.Early","TA.G2")),]
rm(sig1,sig)

# define saving dataframes-------------------------
pval.haber.complete=data.frame(matrix(rep(0,length(unique(sig.complete$celltype))*length(organoids)),
                                      nrow=length(unique(sig.complete$celltype)),ncol = length(organoids)))
rownames(pval.haber.complete)=unique(sig.complete$celltype)
colnames(pval.haber.complete)=organoids
padj.haber.complete=NES.haber.complete=pval.haber.complete


# DE and GSEA-------------------------
for (Smple in organoids) {
  print(Smple)
  data.m = iClassification.class.data.to.subset.cols(data, "genotype", c(Smple, "A")) 
  exprs.m=as.matrix(data.m$data$raw)
  exprs=exprs.m
  rownames(exprs)=data.m$row.annot[rownames(exprs.m),"gene"]
  # DE for Smple vs A-------
  exprs1=exprs
  annot=data.m$col.annot$genotype
  annot[annot==Smple]="Smple"
  
  col.annot=data.frame(annot=annot)
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = exprs1,
                                colData = col.annot,
                                design= ~ annot)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds,name = "annot_Smple_vs_A")
  ranks=res$log2FoldChange 
  ranks[is.na(ranks)]=0
  names(ranks)=rownames(res)
  ranks=sort(ranks,decreasing = T)
  

  # GSEA for Haber celltypes in Smple vs WT------
  library(fgsea)
  library(ggplot2)
  pw=list()
  for (i in 1:length(unique(sig.complete$celltype))) {
    pw[[i]]=sig.complete$gene[sig.complete$celltype==unique(sig.complete$celltype)[i]]
    names(pw)[i]=unique(sig.complete$celltype)[i]
  }
  set.seed(42)
  fgseaRes <- data.frame(fgsea(pw, ranks, minSize=5, maxSize = 1000))
  rownames(fgseaRes)=fgseaRes$pathway
  pval.haber.complete[,Smple]=fgseaRes[rownames(pval.haber.complete),"pval"]
  padj.haber.complete[,Smple]=fgseaRes[rownames(pval.haber.complete),"padj"]
  NES.haber.complete[,Smple]=fgseaRes[rownames(pval.haber.complete),"NES"]
  

}
# plots------------
library(ComplexHeatmap)

CMS=data$col.annot[,c("genotypeCMS","genotype")]
CMS=CMS[!duplicated(CMS$genotype),]
rownames(CMS)=CMS$genotype
colnames(CMS)=c("CMS","Genotype")
CMS=CMS[colnames(pval.haber.complete),]

pval=pval.haber.complete
pval[pval<0.0001]=0.0001
pval=as.matrix(pval)
ES=as.matrix(NES.haber.complete)
ES=ES[c("Enterocyte","Enterocyte.Mature.Proximal","Enterocyte.Mature.Distal",
        "Enterocyte.Progenitor.Late","Enteroendocrine",
        "Paneth","Stem","Tuft","Goblet"),]
pval=pval[c("Enterocyte","Enterocyte.Mature.Proximal","Enterocyte.Mature.Distal",
            "Enterocyte.Progenitor.Late","Enteroendocrine",
            "Paneth","Stem","Tuft","Goblet"),]
col_fun = circlize::colorRamp2(c(-2,-0.9,0.9,2), c("darkblue","white","white","#B51000"))
h1=Heatmap(ES, name = "Percentage", col = col_fun, rect_gp = gpar(type = "none"), 
           width = ncol(ES)*unit(10, "mm"), 
           height = nrow(ES)*unit(10, "mm"),
           cell_fun = function(jh, ih, x, y, width, height, fill) {
             grid.rect(x = x, y = y, width = width, height = height, 
                       gp = gpar(col = NA, fill = col_fun(ES[ih, jh])))
             if(ih > 0&jh>0) {
               grid.circle(x = x, y = y, r = abs(log10(pval)[ih, jh])/9*unit.c(width) , 
                           gp = gpar(fill = col_fun(ES[ih, jh]), col = "black"))
             }},
           show_heatmap_legend = F,
           cluster_rows = F, cluster_columns =T,cluster_column_slices = F,
           top_annotation = HeatmapAnnotation(CMS=as.matrix(CMS$CMS),
                                              col = list(CMS=cms.cols),
                                              # gp = gpar(col = "black"),
                                              simple_anno_size = unit(0.5, "cm"),
                                              show_legend = c(T),
                                              gap = unit(0.1, "cm")),
           column_split = CMS$CMS,
           border = TRUE)
lgd=Legend(title = "Normalized Enrichment Score",col_fun = col_fun,title_position = "leftcenter-rot")
draw(h1,heatmap_legend_list=lgd)


