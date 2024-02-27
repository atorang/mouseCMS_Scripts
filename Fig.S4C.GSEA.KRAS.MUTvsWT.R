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
# define GSEplot function--------
# This is an editted version of plotEnrichment
GSEplot<-function (pathway, stats, gseaParam = 1, ticksSize = 0.2) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    # geom_point(color = "black",size = 0.1) + 
    # geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    # geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0,colour = "black") +
    geom_vline(xintercept = 0,colour = "black") +
    geom_line(color = "darkblue",size = 5*ticksSize) + #Enrichment Color
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y =element_text(colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.ticks = element_blank())+
    # scale_y_continuous(breaks=5,limits=c(min(bottoms), max(tops)))+
    geom_segment(data = data.frame(x = pathway),
                 mapping = aes(x = x, y = -diff+ min(bottoms)-diff/8, xend = x, yend = min(bottoms)-diff/8),
                 size = ticksSize) + 
    theme(panel.border = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "", y = "Enrichment score")
  g
}

# load data---------
load(paste0(PathToData, "/Synapse.RData"))
data=iClassification.class.data.to.subset.cols(data,"CMS","CMS3")
data=iClassification.class.data.to.subset.cols(data,"KRAS",c(0,1))
exprs=data$data$normal
exprs=exprs[!duplicated(data$row.annot$normal$gene),]
rownames(exprs)=data$row.annot$normal[rownames(exprs),"gene"]
# load pathways-----------
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
names(pw)=sigs.entero

# DE, kras mut vs wt--------
library(limma)
class=as.character(data$col.annot$KRAS)
design = model.matrix(~ 0+ class)
colnames(design) = gsub("clas","",colnames(design))
contrast.matrix = makeContrasts(s1-s0,
                                levels=design)
fit.class = contrasts.fit(lmFit(exprs, design), contrast.matrix)
fit.class = eBayes(fit.class)
tt=topTable(fit.class, number=nrow(fit.class), coef=1)

# GSEA------
library(fgsea)
library(ggplot2)

ranks=tt$logFC
names(ranks)=rownames(tt)
ranks=sort(ranks,decreasing = T)

fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway

pdf("Fig.S4C.GSEA.KRAS.MUTvsWT.pdf",height =3,4,family = "sans", pointsize=8) #total size 6.85*11

for(i in rownames(fgseaRes)[order(fgseaRes$padj)]){
  print(GSEplot(pw[[i]], ranks) + 
          labs(title=paste0("CMS3, mutant KRAS vs WT","\n",
                            "NES=",round(fgseaRes[i,"NES"],2),
                            ", PVal=",round(fgseaRes[i,"pval"],4),
                            "\n",i))
  )  
}
dev.off()
