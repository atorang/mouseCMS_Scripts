# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load data----------
load(file=paste0(PathToData,"/scRNAseq.colon.GSE125970.RData"))

# load Metabolic pathways---------
library(fgsea)
library(ggplot2)
res=list()
de.res=list()
pw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))
pw=pw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
        "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
        "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
        "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
        "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
        "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
        "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
        "KEGG_LINOLEIC_ACID_METABOLISM")]
# Finding differentially expressed features (This will take long time to run !!!!!!!!!)----------------------------
library("Seurat")
celltypes=setdiff(unique(annot.colon$CellType),"Progenitor")
for (cell in celltypes) {
  print(cell)
  group1=names(Idents(data ))[annot.colon[names(Idents(data )),"CellType"]==cell]
  progen=names(Idents(data ))[annot.colon[names(Idents(data )),"CellType"]=="Progenitor"]
  group2=setdiff(names(Idents(data )),c(group1,progen))
  de <- FindMarkers(data, ident.1 =group1,  ident.2 =group2, min.pct = 0,thresh.test =0,logfc.threshold = 0)
  head(de, n = 5)
  de.res[[cell]]=de
  logfc=de$avg_log2FC
  names(logfc)=rownames(de)
  
  ranks=logfc
  ranks=sort(ranks,decreasing = T)
  
  set.seed(42)
  fgseaRes <- fgsea(pw, ranks, minSize=5, maxSize = 2000)
  
  enriched=which(fgseaRes$pval<0.05)
  enriched=enriched[order(fgseaRes$padj[fgseaRes$pval<0.05])]
  enriched=enriched[!is.na(enriched)]
    res[[cell]]=fgseaRes

}

# summary result------------------
path=res$Enterocyte$pathway
res.matrix=data.frame(row.names = path,Goblet=rep(0,length(path)))
p.matrix=data.frame(row.names = path,Goblet=rep(0,length(path)))
for (i in 1:length(res)) {
  x=data.frame(res[[i]])
  rownames(x)=x$pathway
  res.matrix[,i]=x[path,"NES"]
  p.matrix[,i]=x[path,"pval"]
  colnames(res.matrix)[i]=colnames(p.matrix)[i]=names(res)[i]
}
res.matrix=res.matrix[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                        "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
                        "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
                        "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
                        "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
                        "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
                        "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
                        "KEGG_LINOLEIC_ACID_METABOLISM"),]
p.matrix=p.matrix[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                    "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
                    "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
                    "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
                    "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
                    "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
                    "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
                    "KEGG_LINOLEIC_ACID_METABOLISM"),]
rownames(res.matrix)=rownames(p.matrix)=c("Amino sugar, nucleotide sugar","Pentose, Glucuronate","Fructose, mannose","Galactose","Glutamine","Glutathione","Nitrogen",
                                          "Glycerophospholipid","Fatty acid","Starch, sucrose","Tyrosine",
                                          "Arachidonic acid","Linoeic acid")
# plot ----------
col_fun = circlize::colorRamp2(c(-2,-0.9,0.9,2), c("darkblue","white","white","#B51000"))
h1=Heatmap(as.matrix(res.matrix), name = "Percentage", col = col_fun, rect_gp = gpar(type = "none"), 
           width = ncol(res.matrix)*unit(10, "mm"), 
           height = nrow(res.matrix)*unit(10, "mm"),
           cell_fun = function(jh, ih, x, y, width, height, fill) {
             grid.rect(x = x, y = y, width = width, height = height, 
                       gp = gpar(col = NA, fill = col_fun(res.matrix[ih, jh])))
             if(ih > 0&jh>0) {
               grid.circle(x = x, y = y, r = unit(abs(log10(p.matrix)[ih, jh])/9*unit.c(width),"mm") , 
                           gp = gpar(fill = NA, col = "black"))
             }},
           show_heatmap_legend = FALSE,
           cluster_rows = FALSE, cluster_columns = FALSE)
lgd=Legend(title = "Normalized Enrichment Score",col_fun = col_fun,title_position = "leftcenter-rot")
draw(h1,heatmap_legend_list=lgd)
