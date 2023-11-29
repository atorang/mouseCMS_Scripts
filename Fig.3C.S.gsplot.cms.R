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


# This is an edited version of plotEnrichment
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
    geom_hline(yintercept = 0,colour = "black") +
    geom_vline(xintercept = 0,colour = "black") +
    geom_line(color = "darkblue",size = 5*ticksSize) + #Enrichment Color
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y =element_text(colour = "black", margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.ticks = element_blank())+
    geom_segment(data = data.frame(x = pathway),
                 mapping = aes(x = x, y = -diff+ min(bottoms)-diff/8, xend = x, yend = min(bottoms)-diff/8),
                 size = ticksSize) + 
    theme(panel.border = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "", y = "Enrichment score")
  g
}

# load pw----------
library(fgsea)
library(ggplot2)
synapse.pw=gmtPathways(paste0(PathToData,"/Synapse.genesets.gmt"))
allpw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))

pw=c((synapse.pw["WNT_FLIER"]),
     (allpw[c("BIOCARTA_WNT_PATHWAY","HALLMARK_WNT_BETA_CATENIN_SIGNALING",
              "KEGG_WNT_SIGNALING_PATHWAY","REACTOME_SIGNALING_BY_WNT",
              "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY")]),
     (synapse.pw["MYC_TARGETS_ZELLER"]),
     (allpw[c("PID_MYC_PATHWAY",
              "BIOCARTA_SRCRPTP_PATHWAY",
              "REACTOME_TRANSLATION")]),
     
     (allpw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
              "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
              "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
              "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
              "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
              "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
              "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
              "KEGG_LINOLEIC_ACID_METABOLISM")]),
     
     (synapse.pw[c("MESENCH_LOBODA","EMT_CORE_GENES")]),
     (allpw[c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "KEGG_TGF_BETA_SIGNALING_PATHWAY",
              "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
              "GOBP_RESPONSE_TO_WOUNDING")]),
     (synapse.pw["STROMAL_ESTIMATE"]),             
     (allpw[c("REACTOME_INTEGRIN_SIGNALING",
              "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS")])
     
)

names(pw)=gsub("MESENCH_LOBODA","Mesenchymal",names(pw))
names(pw)=gsub("STROMAL_ESTIMATE","Stromal infiltration",names(pw))
names(pw)=gsub("WNT_FLIER","WNT targets(Flier)",names(pw))
names(pw)=gsub("BIOCARTA_WNT_PATHWAY","WNT targets(Biocarta)",names(pw))
names(pw)=gsub("HALLMARK_WNT_BETA_CATENIN_SIGNALING","WNT targets(Hallmark)",names(pw))
names(pw)=gsub("KEGG_WNT_SIGNALING_PATHWAY","WNT targets(Kegg)",names(pw))
names(pw)=gsub("REACTOME_SIGNALING_BY_WNT","WNT targets(Reactome)",names(pw))
names(pw)=gsub("GOBP_CANONICAL_WNT_SIGNALING_PATHWAY","WNT targets(GOBP)",names(pw))
names(pw)=gsub("MYC_TARGETS_ZELLER","MYC targets (Zeller)",names(pw))
names(pw)=gsub("PID_MYC_PATHWAY","MYC targets(PID)",names(pw))
names(pw)=gsub("EMT_CORE_GENES","EMT activation",names(pw))
names(pw)=gsub("KEGG_TGF_BETA_SIGNALING_PATHWAY","TGFB activation",names(pw))
names(pw)=gsub("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","Matrix remodeling",names(pw))
names(pw)=gsub("GOBP_RESPONSE_TO_WOUNDING","Wound response",names(pw))
names(pw)=gsub("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE","Immune response",names(pw))
names(pw)=gsub("REACTOME_PD__SIGNALING","PD1 activation",names(pw))
names(pw)=gsub("GOBP_COMPLEMENT_ACTIVATION","Complement activation",names(pw))

names(pw)=gsub("BIOCARTA_SRCRPTP_PATHWAY","SRC",names(pw))
names(pw)=gsub("KEGG_JAK_STAT_SIGNALING_PATHWAY","JAK STAT",names(pw))
names(pw)=gsub("REACTOME_TRANSLATION","Transation ribosome",names(pw))
names(pw)=gsub("REACTOME_INTEGRIN_SIGNALING","Integrin B3",names(pw))
names(pw)=gsub("REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS","VEGF, VEGFR",names(pw))

names(pw)=gsub("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM","Amino sugar, nucleotide sugar",names(pw))
names(pw)=gsub("KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS","Pentose, glucuronate",names(pw))
names(pw)=gsub("KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM","Fructose, mannose",names(pw))
names(pw)=gsub("KEGG_GALACTOSE_METABOLISM","Galactose",names(pw))
names(pw)=gsub("GOBP_GLUTAMINE_METABOLIC_PROCESS","Glutamine",names(pw))
names(pw)=gsub("KEGG_GLUTATHIONE_METABOLISM","Glutathione",names(pw))
names(pw)=gsub("KEGG_NITROGEN_METABOLISM","Nitrogen",names(pw))
names(pw)=gsub("KEGG_GLYCEROPHOSPHOLIPID_METABOLISM","Glycerophospholipid",names(pw))
names(pw)=gsub("KEGG_FATTY_ACID_METABOLISM","Fatty acid",names(pw))
names(pw)=gsub("KEGG_STARCH_AND_SUCROSE_METABOLISM","Starch, sucrose",names(pw))
names(pw)=gsub("KEGG_ARACHIDONIC_ACID_METABOLISM","Arachidonic acid",names(pw))
names(pw)=gsub("KEGG_LINOLEIC_ACID_METABOLISM","Linoleic acid",names(pw))
names(pw)=gsub("KEGG_TYROSINE_METABOLISM","Tyrosine",names(pw))


# ---------------SIP------------
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 
data$data$raw=data$data$raw[,order(data$col.annot$genotype)]
data$data$quantile=data$data$quantile[,order(data$col.annot$genotype)]
data$col.annot=data$col.annot[order(data$col.annot$genotype),]
data$col.annot$CMS=data$col.annot$genotypeCMS

data$data$raw=data$data$raw[,order(data$col.annot$CMS)]
data$data$quantile=data$data$quantile[,order(data$col.annot$CMS)]
data$col.annot=data$col.annot[order(data$col.annot$CMS),]

# orthologs---------
exprs.m=as.matrix(data$data$raw)
mmm=apply(exprs.m,1,mean)
exprs.m=exprs.m[log2(mmm+1)>0,]
exprs=exprs.m

load(paste0(PathToData,"/orth.v105.RData"))

rownames(exprs) = orth$Gene.name.1[match(rownames(exprs), orth$Gene.stable.ID)]
exprs = exprs[!is.na(rownames(exprs)),]
exprs = exprs[apply(exprs,1,function(x) sum(is.na(x))==0),]
exprs=exprs[!duplicated(rownames(exprs)),]



colnames(exprs)=gsub("F1.|.SIP|.CO","",colnames(exprs))

# DE, CMS2 vs CMS3,4 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS2"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms2=ranks=sort(ranks,decreasing = T)

# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms2=fgseaRes

# DE, CMS3 vs CMS2,4 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS3"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms3=ranks=sort(ranks,decreasing = T)



# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms3=fgseaRes

# DE, CMS4 vs CMS2,3 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS4"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms4=ranks=sort(ranks,decreasing = T)



# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms4=fgseaRes

# GSEA plots--------------

for(i in names(pw)){
  print(GSEplot(pw[[i]], ranks.cms2) + 
          labs(title=paste0("CMS2 vs CMS3,4","\n",
                            "NES=",round(fgseaRes.cms2[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms2[i,"pval"],4),
                            "\n",i))
  )
  print(GSEplot(pw[[i]], ranks.cms3) + 
          labs(title=paste0("CMS3 vs CMS2,4","\n",
                            "NES=",round(fgseaRes.cms3[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms3[i,"pval"],4),
                            "\n",i))
  )
  print(GSEplot(pw[[i]], ranks.cms4) + 
          labs(title=paste0("CMS4 vs CMS2,3","\n",
                            "NES=",round(fgseaRes.cms4[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms4[i,"pval"],4),
                            "\n",i))
  )
}

# ---------------CO------------
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "CO") 
data$data$raw=data$data$raw[,order(data$col.annot$genotype)]
data$data$quantile=data$data$quantile[,order(data$col.annot$genotype)]
data$col.annot=data$col.annot[order(data$col.annot$genotype),]
data$col.annot$CMS=data$col.annot$genotypeCMS

data$data$raw=data$data$raw[,order(data$col.annot$CMS)]
data$data$quantile=data$data$quantile[,order(data$col.annot$CMS)]
data$col.annot=data$col.annot[order(data$col.annot$CMS),]

# orthologs---------
exprs.m=as.matrix(data$data$raw)
mmm=apply(exprs.m,1,mean)
exprs.m=exprs.m[log2(mmm+1)>0,]
exprs=exprs.m

load(paste0(PathToData,"/orth.v105.RData"))

rownames(exprs) = orth$Gene.name.1[match(rownames(exprs), orth$Gene.stable.ID)]
exprs = exprs[!is.na(rownames(exprs)),]
exprs = exprs[apply(exprs,1,function(x) sum(is.na(x))==0),]
exprs=exprs[!duplicated(rownames(exprs)),]



colnames(exprs)=gsub("F1.|.SIP|.CO","",colnames(exprs))

# DE, CMS2 vs CMS3,4 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS2"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms2=ranks=sort(ranks,decreasing = T)

# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms2=fgseaRes

# DE, CMS3 vs CMS2,4 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS3"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms3=ranks=sort(ranks,decreasing = T)



# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms3=fgseaRes

# DE, CMS4 vs CMS2,3 using Deseq2--------
col.annot=data$col.annot
rownames(col.annot)=colnames(exprs)
exprs1=exprs[,col.annot$CMS%in%paste0("CMS",2:4)]
annot=col.annot$CMS[col.annot$CMS%in%paste0("CMS",2:4)]
annot[annot=="CMS4"]="Smple"
annot[annot!="Smple"]="tumor"
col.annot=col.annot[col.annot$CMS%in%paste0("CMS",2:4),]
col.annot$annot=annot

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs1,
                              colData = col.annot,
                              design= ~ batch+annot)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,contrast=c("annot","Smple","tumor"))
ranks=res$log2FoldChange 
ranks[is.na(ranks)]=0
names(ranks)=rownames(res)
ranks.cms4=ranks=sort(ranks,decreasing = T)



# GSEA-------------------------
fgseaRes <- data.frame(fgsea(pw, ranks, minSize=3, maxSize = 1000))
rownames(fgseaRes)=fgseaRes$pathway
fgseaRes.cms4=fgseaRes

# GSEA plots--------------

for(i in names(pw)){
  print(GSEplot(pw[[i]], ranks.cms2) + 
          labs(title=paste0("CMS2 vs CMS3,4","\n",
                            "NES=",round(fgseaRes.cms2[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms2[i,"pval"],4),
                            "\n",i))
  )
  print(GSEplot(pw[[i]], ranks.cms3) + 
          labs(title=paste0("CMS3 vs CMS2,4","\n",
                            "NES=",round(fgseaRes.cms3[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms3[i,"pval"],4),
                            "\n",i))
  )
  print(GSEplot(pw[[i]], ranks.cms4) + 
          labs(title=paste0("CMS4 vs CMS2,3","\n",
                            "NES=",round(fgseaRes.cms4[i,"NES"],2),
                            ", PVal=",round(fgseaRes.cms4[i,"pval"],4),
                            "\n",i))
  )
}
