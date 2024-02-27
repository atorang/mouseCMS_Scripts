# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# define function--------
iClassification.best.gene = function(data, data.type, feature="probeset", na.rm=F) {#In repeated genes, those with max mean will be kept
  .check.data.init(data, "iClassification.best.gene")
  dup.genes = names(table(data$row.annot[,"gene"])[table(data$row.annot[,"gene"])>1]) #Find genes with more than 1 probes
  dup.probes.exprs = apply(data$data[[data.type]][data$row.annot[,"gene"] %in% dup.genes,],1,mean)
  remove.probes = apply(as.matrix(dup.genes),1, function(x) {  
    probes = data$row.annot[which(data$row.annot[,"gene"] == x),feature]
    probes[-which.max(dup.probes.exprs[probes])] })
  remove.probes = unlist(remove.probes)
  remove.probes = c(remove.probes,data$row.annot[is.na(data$row.annot[,"gene"]),feature])
  data = iClassification.class.data.to.subset.rows(data, feature, data$row.annot[!data$row.annot[,feature] %in% remove.probes,feature])
  return(data)
}
# load data----------
load(paste0(PathToData, "/GSE100478.RData"))
data=iClassification.best.gene(data,"rma")
exprs=data$data$rma
rownames(exprs)=data$row.annot[,"gene"]
colnames(exprs)=data$col.annot$CellLine
# load signatures----------
library(GSVA)
library(fgsea)
pw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))
sigs.entero=names(pw)[grep("ENTEROCYTE",names(pw))]
pw=pw[sigs.entero]
sigs.entero=c("Duodenal early immature enterocyte",
              "Duodenal late immature enterocyte",
              "Duodenal mature enterocyte",
              "SI_24W_C3 enterocyte progenitor type1",
              "SI_24W_C4 enterocyte progenitor type2",
              "Colon_24W_C10 enterocyte")
names(pw)= sigs.entero

# ssGSEA------------
gsva.es <- gsva(exprs, pw, verbose=FALSE,method="ssgsea") #"gsva", "ssgsea", "zscore", "plage"
x=apply(gsva.es, 2,sum)
breaks=seq(-0.5,0.5,0.01)

gsva.es.z=gsva.es
for (i in 1:nrow(gsva.es.z)) {
  gsva.es.z[i,]=(gsva.es.z[i,]-mean(gsva.es.z[i,]))/sd(gsva.es.z[i,])
}

x=apply(gsva.es.z, 2,sum)
breaks=seq(-2,2,0.01)
col=colorRampPalette(c("goldenrod1","white","darkorchid4"))(length(breaks))

library(grid) 
library(pheatmap)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

pdf("Fig.6B.ssGSEA.CellLines.pdf",width = 8,height = 2.5)
pheatmap(gsva.es.z[,order(x)],
         breaks=breaks,col=col,cluster_cols = F,border_color = "NA")
dev.off()
