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
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 

data$data$raw=data$data$raw[,order(data$col.annot$genotype)]
data$data$quantile=data$data$quantile[,order(data$col.annot$genotype)]
data$col.annot=data$col.annot[order(data$col.annot$genotype),]

exprs=data$data$quantile
# orthologs---------
exprs.m=as.matrix(data$data$quantile)
mmm=apply(exprs.m,1,function(x){sum(x>0)})
exprs.m=exprs.m[mmm>=3,]
exprs=exprs.m
dim(exprs)

load(paste0(PathToData,"/orth.v105.RData")) 

rownames(exprs) = orth$Gene.name.1[match(rownames(exprs), orth$Gene.stable.ID)]
exprs = exprs[!is.na(rownames(exprs)),]
exprs = exprs[apply(exprs,1,function(x) sum(is.na(x))==0),]
exprs=exprs[!duplicated(rownames(exprs)),]


colnames(exprs)=gsub("F1.|.SIP","",colnames(exprs))
# Labeling mouse samples-------------
library(GSVA)
synapse.pw=gmtPathways(paste0(PathToData,"/Synapse.genesets.gmt"))
allpw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))

#include no sig with positive score in two subtypes
pw=c(synapse.pw["WNT_FLIER"],
     allpw[c("BIOCARTA_WNT_PATHWAY","KEGG_WNT_SIGNALING_PATHWAY")],
     synapse.pw["MYC_TARGETS_ZELLER"],
     allpw[c("PID_MYC_PATHWAY",
             "BIOCARTA_SRCRPTP_PATHWAY",
             "REACTOME_TRANSLATION")],
     
     allpw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
             "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
             "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
             "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
             "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
             "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
             "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
             "KEGG_LINOLEIC_ACID_METABOLISM")],
     
     synapse.pw[c("MESENCH_LOBODA","EMT_CORE_GENES")],
     allpw[c("KEGG_TGF_BETA_SIGNALING_PATHWAY",
             "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
             "GOBP_RESPONSE_TO_WOUNDING")],
     synapse.pw["STROMAL_ESTIMATE"],             
     allpw[c("REACTOME_INTEGRIN_SIGNALING",
             "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS"   
     )]
)

CMS.lab=c(rep("CMS2",7),rep("CMS3",13),rep("CMS4",8))
# ssgsea----------
res=c()
genotype=c()
MAX=c()
SecMAX=c()
probScores=c()
listA=c("A","AKP","AP","AP2","AKP2","AKPS")
listB=c("APN","AKPN","BPN","BPNA2","KPN","KPNA2")
listC=c("K","KP","BP","AK")
for (smpl in listA) {
  #smpl="A"
  exprs2=exprs[,data$col.annot$genotype%in%c(smpl,listB,listC)]
  gsva.es <- gsva(exprs2, pw, verbose=FALSE,method="zscore") #"gsva", "ssgsea", "zscore"
  final1=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x>0.5)/length(x)})[,-1]
  final2=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x<(-0.5))/length(x)})[,-1]
  final=final1-final2
  rownames(final)=paste0("CMS",2:4)
  
  genom=gsub("\\..*","",colnames(final)[-1])
  final=final[,-1]
  final=aggregate(t(final), by=list(genom), FUN = mean)
  rownames(final)=final$Group.1
  m=apply(final[,-1],1,max)
  m2=apply(final[,-1],1,function(x){max(x[-which(x==max(x))])})
  
  x=colnames(final)[final[smpl,]==m[smpl]]
  if(length(x)==1){res=c(res,x)}else{res=c(res,NA)}
  genotype=c(genotype,smpl)
  MAX=c(MAX,m[smpl])
  SecMAX=c(SecMAX,m2[smpl])
  probScores=rbind(probScores,final[smpl,])
  rownames(probScores)[nrow(probScores)]=smpl
}
for (smpl in listB) {
  exprs2=exprs[,data$col.annot$genotype%in%c(smpl,listA,listC)]
  gsva.es <- gsva(exprs2, pw, verbose=FALSE,method="zscore") #"gsva", "ssgsea", "zscore"
  final1=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x>0.5)/length(x)})[,-1]
  final2=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x<(-0.5))/length(x)})[,-1]
  final=final1-final2
  rownames(final)=paste0("CMS",2:4)
  
  genom=gsub("\\..*","",colnames(final)[-1])
  final=final[,-1]
  final=aggregate(t(final), by=list(genom), FUN = mean)
  rownames(final)=final$Group.1
  m=apply(final[,-1],1,max)
  m2=apply(final[,-1],1,function(x){max(x[-which(x==max(x))])})
  
  x=colnames(final)[final[smpl,]==m[smpl]]
  if(length(x)==1){res=c(res,x)}else{res=c(res,NA)}
  genotype=c(genotype,smpl)
  MAX=c(MAX,m[smpl])
  SecMAX=c(SecMAX,m2[smpl])
  probScores=rbind(probScores,final[smpl,])
  rownames(probScores)[nrow(probScores)]=smpl
}
for (smpl in listC) {
  exprs2=exprs[,data$col.annot$genotype%in%c(smpl,listA,listB)]
  gsva.es <- gsva(exprs2, pw, verbose=FALSE,method="zscore") #"gsva", "ssgsea", "zscore"
  final1=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x>0.5)/length(x)})[,-1]
  final2=aggregate(gsva.es, by=list(CMS.lab), FUN = function(x){sum(x<(-0.5))/length(x)})[,-1]
  final=final1-final2
  rownames(final)=paste0("CMS",2:4)
  
  genom=gsub("\\..*","",colnames(final)[-1])
  final=final[,-1]
  final=aggregate(t(final), by=list(genom), FUN = mean)
  rownames(final)=final$Group.1
  m=apply(final[,-1],1,max)
  m2=apply(final[,-1],1,function(x){max(x[-which(x==max(x))])})
  
  x=colnames(final)[final[smpl,]==m[smpl]]
  if(length(x)==1){res=c(res,x)}else{res=c(res,NA)}
  genotype=c(genotype,smpl)
  MAX=c(MAX,m[smpl])
  SecMAX=c(SecMAX,m2[smpl])
  probScores=rbind(probScores,final[smpl,])
  rownames(probScores)[nrow(probScores)]=smpl
}
res=data.frame(names = genotype,CMS=res,Max=MAX,SecMAX=SecMAX,dif=(MAX-SecMAX))
res$v=res$dif/res$Max
writexl::write_xlsx(data.frame(samples=rownames(res),res),"First.CMS.Stratification.xlsx")

