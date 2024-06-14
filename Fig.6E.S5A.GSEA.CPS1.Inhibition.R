# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load pathways-----------
library(fgsea)
pw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))
pw=pw[c("KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
        "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
        "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
        "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
        "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
        "KEGG_LINOLEIC_ACID_METABOLISM",
        "KEGG_CELL_CYCLE","KEGG_DNA_REPLICATION")]

names(pw)=c("Fructose, mannose","Galactose","Glutamine","Glutathione","Nitrogen",
            "Glycerophospholipid","Fatty acid","Starch, sucrose","Tyrosine",
            "Arachidonic acid","Linoeic acid",
            "Cell cycle (KEGG)","DNA replication (KEGG)")

# Ref: Wang Y, Song W, Wang J, Wang T, Xiong X, Qi Z, Fu W, Yang X, Chen YG. 
# Single-cell transcriptome analysis reveals differential nutrient absorption functions in human intestine. 
# J Exp Med. 2020 Feb 3;217(2):e20191130. doi: 10.1084/jem.20191130. PMID: 31753849; PMCID: PMC7041720.
# Table S2
sig=data.frame(readxl::read_excel(paste0(PathToData,"jem_20191130_tables2.xlsx"),sheet = "rectum"))
pw2=list()
for (i in 1:length(unique(sig$cluster))) {
  pw2[[i]]=sig$gene[sig$cluster==unique(sig$cluster)[i]]
  names(pw2)[i]=unique(sig$cluster)[i]
}
pw=c(pw,pw2)
# GSEplot function--------
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
load(paste0(PathToData, "/CellLine.GSE248558.RData"))
exprs=data$data$raw
exprs=exprs[!is.na(data$row.annot$gene)&!duplicated(data$row.annot$gene),]
rownames(exprs)=data$row.annot[rownames(exprs),"gene"]

mmm=apply(exprs,1,max)
exprs=exprs[mmm>0,]
# DE-------
library(DESeq2)
col.annot=data$col.annot
for (i in 1:ncol(col.annot)) {
  col.annot[,i]=as.factor(col.annot[,i])
}
dds <- DESeqDataSetFromMatrix(countData = exprs,
                              colData = col.annot,
                              design= ~ treatment+CellLine)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds,contrast = c("treatment","100umH3B120","control"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks=sort(ranks,decreasing = T)
# GSEA------
library(fgsea)
library(ggplot2)
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway

pdf("Fig.6E.S5A.GSEA.CPS1.Inhibition.pdf",height =3,4,family = "sans", pointsize=8) #total size 6.85*11
for(i in rownames(fgseaRes)[order(fgseaRes$padj)]){
  print(GSEplot(pw[[i]], ranks) + 
          labs(title=paste0("100umH3B120 vs control","\n",
                            "NES=",round(fgseaRes[i,"NES"],2),
                            ", PVal=",round(fgseaRes[i,"pval"],4),
                            "\n",i))
  )  
}
dev.off()
