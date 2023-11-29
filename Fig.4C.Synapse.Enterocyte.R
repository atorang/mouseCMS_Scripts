# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# define function--------
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
load(paste0(PathToData,"/Synapse.RData"))
exprs=data$data$normal
exprs=exprs[!duplicated(data$row.annot$gene),]
rownames(exprs)=data$row.annot[rownames(exprs),"gene"]

exprs=data.frame(exprs)
col.annot=data.frame(row.names = colnames(exprs),
                     CMS= c(data$col.annot$CMS))


# DE, CMS3 vs CMS1,2,4 using limma--------
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",1:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",1:4)]
annot[annot=="CMS3"]="Smple"
annot[annot!="Smple"]="tumor"

library(limma)
design <- model.matrix(~ 0+(annot))
colnames(design) = substr(colnames(design),nchar(colnames(design))-4,nchar(colnames(design)))
contrast.matrix = makeContrasts(Smple-tumor,levels=design)
fit = contrasts.fit(lmFit(exprs1, design), contrast.matrix)
fit <- eBayes(fit)
tt1=topTable(fit,number=nrow(exprs1), coef=1)
logfc=tt1$logFC
names(logfc)=rownames(tt1)

ranks=logfc
ranks=sort(ranks,decreasing = T)
cms3.rank=ranks
# GSEA-------------------------
library(fgsea)
library(ggplot2)
pw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))
pw=pw[names(pw)[grep("ENTEROCYTE", names(pw))]]
names(pw)=c("Duodenal early immature enterocytes","Duodenal late immature enterocytes",
            "Duodenal mature enterocytes","SI enterocyte progenitor 1",
            "SI enterocyte progenitor 2","CO enterocyte")

cms3.fgseaRes <- data.frame(fgsea(pw, cms3.rank, minSize=5, maxSize = 1000))
rownames(cms3.fgseaRes)=cms3.fgseaRes$pathway

# plot-------------
for (i in names(pw)) {
  plot(GSEplot(pw[[i]], cms3.rank) + 
         labs(title=paste0("CMS3 vs CMS1,2,4","\n",
                           "NES=",round(cms3.fgseaRes[i,"NES"],2),
                           ", PVal=",round(cms3.fgseaRes[i,"pval"],4),
                           "\n",i)))
}

