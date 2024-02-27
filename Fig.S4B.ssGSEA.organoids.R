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
# load data----------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP")
# load signatures----------
library(GSVA)
library(fgsea)
pw=gmtPathways(paste0(my.data.dir,"General/msigdb.v7.4.symbols.gmt"))
sigs.entero=names(pw)[grep("ENTEROCYTE",names(pw))]
pw=pw[c("KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
        "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
        "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
        "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
        "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
        "KEGG_LINOLEIC_ACID_METABOLISM",sigs.entero)]
sigs.entero=c("Duodenal early immature enterocyte",
              "Duodenal late immature enterocyte",
              "Duodenal mature enterocyte",
              "SI_24W_C3 enterocyte progenitor type1",
              "SI_24W_C4 enterocyte progenitor type2",
              "Colon_24W_C10 enterocyte")
names(pw)=c("Fructose, mannose","Galactose","Glutamine","Glutathione","Nitrogen",
            "Glycerophospholipid","Fatty acid","Starch, sucrose","Tyrosine",
            "Arachidonic acid","Linoeic acid",sigs.entero)

# orthologs---------
exprs.Alex=as.matrix(data$data$raw)
mmm=apply(exprs.Alex,1,mean)
exprs.Alex=exprs.Alex[log2(mmm+1)>0,]
exprs=exprs.Alex

load(paste0(PathToData,"/orth.v105.RData"))

rownames(exprs) = orth$Gene.name.1[match(rownames(exprs), orth$Gene.stable.ID)]
exprs = exprs[!is.na(rownames(exprs)),]
exprs = exprs[apply(exprs,1,function(x) sum(is.na(x))==0),]
exprs=exprs[!duplicated(rownames(exprs)),]

colnames(exprs)=gsub("F1.|.SIP","",colnames(exprs))
# ssGSEA------------
gsva.es <- gsva(exprs, pw, verbose=FALSE,method="ssgsea") #"gsva", "ssgsea", "zscore", "plage"

gsva.es.z=gsva.es
for (i in 1:nrow(gsva.es.z)) {
  gsva.es.z[i,]=(gsva.es.z[i,]-mean(gsva.es.z[i,]))/sd(gsva.es.z[i,])
}
x=apply(gsva.es.z, 2,sum)
breaks=seq(-2,2,0.01)
col=colorRampPalette(c("goldenrod1","white","darkorchid4"))(length(breaks))

library(grid) 
library(pheatmap)
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))

pdf("Fig.S4B.ssGSEA.organoids.pdf",width = 14,height = 4)
pheatmap::pheatmap(gsva.es.z[,order(x)],breaks=breaks,col=col,cluster_cols = F)
dev.off()
