# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# Functions---------
iClassification.best.entrezId = function(data, data.type, feature="probeset", na.rm=F) {#In repeated genes, those with max mean will be kept
  dup.genes = names(table(data$row.annot[,"entrezId"])[table(data$row.annot[,"entrezId"])>1]) #Find genes with more than 1 probes
  dup.probes.exprs = apply(data$data[[data.type]][data$row.annot[,"entrezId"] %in% dup.genes,],1,mean)
  remove.probes = apply(as.matrix(dup.genes),1, function(x) {  
    probes = data$row.annot[which(data$row.annot[,"entrezId"] == x),feature]
    probes[-which.max(dup.probes.exprs[probes])] })
  remove.probes = unlist(remove.probes)
  remove.probes = c(remove.probes,data$row.annot[is.na(data$row.annot[,"entrezId"]),feature])
  data = iClassification.class.data.to.subset.rows(data, feature, data$row.annot[!data$row.annot[,feature] %in% remove.probes,feature])
  return(data)
}

GEO.Affy<-function(GSEname,
                   directory=paste0(GSEname),
                   output.directory=NA,
                   output.name=paste0(GSEname,".Affy")){
  if (dir.exists(directory)) {stop(paste0("Error: The directory already exists!"))} 
  dir.create(directory,recursive = T)
  library(GEOquery)
  cel.dir = directory
  
  # load series and platform data from GEO
  gseMatrix<-getGEO(GSEname, GSEMatrix=T, destdir =cel.dir)
  gset <- getGEO(GSEname, GSEMatrix=F, destdir =cel.dir)
  gsmlist = GSMList(gset)
  header = lapply(gsmlist, function(x) c(x@header$title, x@header$characteristics_ch1,paste0("source: ",x@header$source_name_ch1),
                                         paste0("PatId: ",x@header$description)))
  
  x=c()
  for(i in 1:length(header)){x=c(x,(header[[i]]))}
  x=gsub("\\:.*","",x)
  n=sum(sort(table(x),decreasing = T)>1)
  annot = matrix(nrow=length(header), ncol=3+n) #ncol=4 + number of papameters
  colnames(annot) = c("GSM","Dataset", "SampleId", names(sort(table(x),decreasing = T))[1:n]) #number of papameters
  rownames(annot) = names(header)
  Check=function(x,info){
    if(length(gsub(x, "", info[grep(x, info)]))==0){
      return(NA)
    }else{
      return(gsub(x, "", info[grep(x, info)]))
    }
  }
  for(i in 1:length(header)){
    info = header[[i]]
    annot[i,1] = names(header)[i]           #GEO name
    annot[i,2] = GSEname                    #GEO name of dataset
    annot[i,3] = gsub("\\:.*","",info[1])   #Sample name
    for (j in 2:length(info)) {
      annot[i,gsub("\\: .*","",info[j])] = gsub(".*\\: ","",info[j])  #Other Col.annot (MAY NEED ADABTATION)
    }
  }
  annot=data.frame(annot,gseMatrix[[1]]@phenoData@data)
  annot = as.data.frame(annot, stringsAsFactors=F)
  
  # Download data 
  options(timeout = 2500)
  filePaths = getGEOSuppFiles(GSEname,baseDir = cel.dir)
  untar(paste0(cel.dir,"/",GSEname,"/",GSEname,"_RAW.tar"),exdir=cel.dir)
  cel.files = list.files(cel.dir)[grep("CEL.gz",list.files(cel.dir))]
  gsm = sub("^(GSM\\d+).*", "\\1", cel.files)
  library(oligo)
  batch = read.celfiles(paste0(cel.dir,"/",cel.files[match(annot$GSM, gsm )])) 
  eset = rma(batch) 
  rma = exprs(eset)
  
  #Delete files
  file.remove(paste0(cel.dir,"/",cel.files))
  unlink(paste0(cel.dir,"/",GSEname), recursive=TRUE)
  
  # row.annot
  genes=as.character(rownames(rma))
  library(annotate)
  library("hgu133plus2.db")
  
  row.annot=select(hgu133plus2.db, genes,c("SYMBOL","ENTREZID","ENSEMBL","CHR"))
  colnames(row.annot)<-c("probeset","gene","entrezId","ensemblId","chromosome")
  row.annot<-row.annot[!duplicated(row.annot$probeset),]
  rownames(row.annot)<-row.annot$probeset
  
  row.annot<-row.annot[!duplicated(row.annot$probeset),]
  rownames(row.annot)<-row.annot$probeset
  row.annot=row.annot[rownames(rma),]
  rownames(row.annot)<-rownames(rma)
  row.annot$probeset<-rownames(rma)
  row.annot=data.frame(row.annot,gseMatrix[[1]]@featureData@data[rownames(rma),])
  
  # save data
  data = list()
  data$data$rma = rma
  data$col.annot = annot
  data$row.annot = row.annot
  colnames(data$data$rma) = data$col.annot$GSM
  data$init = T
  if(is.na(output.directory)){output.directory=directory}
  save(data, file=paste(output.directory,"/",output.name,".RData",sep=""))
}

#---------------------prepare all datasets-------------
# CMS labels and clinical info from Guinney et al.-------------
#download data from Synapse platform (doi:10.7303/syn2623706)
CMS.dir=paste0(PathToDatar,"ColorectalCancerSubtypingConsortium/data/mergedPhenotype/cms_labels_public_all.txt")
cms=read.table(CMS.dir,sep="\t",header = T)

clinical.dir=paste0(PathToDatar,"ColorectalCancerSubtypingConsortium/data/mergedPhenotype/clinical_molecular_public_all.txt")
clinical=read.table(clinical.dir,sep="\t",header = T)
# KFSYSCC-------
KFSYSCC.dir=paste0(PathToData,"ColorectalCancerSubtypingConsortium/data/geneExpression/KFSYSCC")
setwd(KFSYSCC.dir)
rma=read.table("kfsyscc_expression.tsv",header = T)
fRMA=read.table("kfsyscc_frma_expression.tsv",header = T)
rownames(rma)<-rma[,1]
rma<-rma[,-1]
rma=rma[,colnames(fRMA)]
rma=rma[,colnames(fRMA)[16:ncol(fRMA)]]

# col.annot
annot = matrix(nrow=ncol(fRMA), ncol=4)
colnames(annot) = c("SampleId","dataset","rmaSamples","frmaSamples")
rownames(annot) = colnames(fRMA)

annot[,1]=colnames(fRMA)
annot[,2]="KFSYSCC"
annot[16:ncol(fRMA),3]=1
annot[1:15,3]=0
annot[,4]=1
annot = as.data.frame(annot, stringsAsFactors=F)
annot[,3]=as.numeric(annot[,3])
annot[,4]=as.numeric(annot[,4])

# row.annot
library(biomaRt)
mart=useMart("ensembl",dataset="hsapiens_gene_ensembl")
sum(rownames(fRMA)!=rownames(rma))
setdiff(rownames(fRMA),rownames(rma))
fRMA=fRMA[rownames(rma),]

sum(rownames(fRMA)!=rownames(rma))
genes=as.character(rownames(fRMA))
View(listAttributes(mart))
platform="affy_hg_u133_plus_2" 
row.annot1=getBM(attributes = c(platform,"hgnc_symbol","entrezgene_id","ensembl_gene_id","gene_biotype","description","source",
                                "chromosome_name","start_position","end_position","strand"),
                 filters    = platform,
                 values     = as.vector(genes), 
                 mart       = mart)
colnames(row.annot1)<-c("probeset","gene","entrezId","ensemblId","geneBiotype","description","source",
                        "chromosome","start","end","strand")
row.annot1<-row.annot1[!duplicated(row.annot1$probeset),]
rownames(row.annot1)<-row.annot1$probeset
row.annot1=row.annot1[rownames(fRMA),]
rownames(row.annot1)<-rownames(fRMA) 
row.annot1$probeset=rownames(row.annot1)

sum(row.annot1$probeset!=rownames(rma))

# save data
data = list()
data$data$rma = rma
data$data$fRMA = fRMA
data$col.annot = annot
data$row.annot = row.annot1
data$init = T
row.annot1=data$row.annot
genes=rownames(data$data$rma)

library(annotate)
library("hgu133plus2.db")

row.annot=select(hgu133plus2.db, genes,c("SYMBOL","ENTREZID","ENSEMBL","CHR"))
colnames(row.annot)<-c("probeset","gene","entrezId","ensemblId","chromosome")
row.annot<-row.annot[!duplicated(row.annot$probeset),]
rownames(row.annot)<-row.annot$probeset
NAs=which(is.na(row.annot$gene))
NAs=row.annot$probeset[NAs]

data$row.annot=row.annot
save(data, file=paste(PathToData, "Synapse/KFSYSCC.Affy.RData",sep=""))

# PETACC3-------
PETACC3.dir=paste0(PathToData,"CRCSubtypingConsortium/datasets/PETACC3")
setwd(PETACC3.dir)
load("Normals_Processed/tumors_and_normals.probeset.level.matrix.RData")
rma=data.frame(t(X))

# col.annot
annot = matrix(nrow=ncol(rma), ncol=2)
colnames(annot) = c("SampleId","dataset")
rownames(annot) = colnames(rma)

annot[,1]=colnames(rma)
annot[,2]="PETACC3"
annot = as.data.frame(annot, stringsAsFactors=F)

# row.annot
row.annot=read.table("ArrayAnnotation/petacc.platform.annotation.tsv",header = T,sep="\t",fill = T,stringsAsFactors = F)
row.annot[51824,2]="---"
row.annot=row.annot[!duplicated(row.annot$ProbesetID),]
rownames(row.annot)=row.annot[,1]
genes=intersect(rownames(rma),rownames(row.annot))
row.annot=row.annot[genes,]
row.annot=data.frame(row.annot)
rma=rma[genes,]
colnames(row.annot)[c(1,2)]<-c("probeset","entrezId")

# save data
data = list()
data$data$rma = rma
data$col.annot = annot
data$row.annot = row.annot
data$init = T
save(data, file=paste(PathToData, "Synapse/PETACC3.Affy.RData",sep="")) #Specify the dataset name

# GEO data------------
GEO.Affy("GSE33113", directory = paste0(PathToData, "Synapse/GSE33113"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE39582", directory = paste0(PathToData, "Synapse/GSE39582"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE35896", directory = paste0(PathToData, "Synapse/GSE35896"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE13294", directory = paste0(PathToData, "Synapse/GSE13294"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE14333", directory = paste0(PathToData, "Synapse/GSE14333"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE17536", directory = paste0(PathToData, "Synapse/GSE17536"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE20916", directory = paste0(PathToData, "Synapse/GSE20916"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE2109", directory = paste0(PathToData, "Synapse/GSE2109"), output.directory = paste0(PathToData, "Synapse"))   
GEO.Affy("GSE23878", directory = paste0(PathToData, "Synapse/GSE23878"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE37892", directory = paste0(PathToData, "Synapse/GSE37892"), output.directory = paste0(PathToData, "Synapse")) 
GEO.Affy("GSE13067", directory = paste0(PathToData, "Synapse/GSE13067"), output.directory = paste0(PathToData, "Synapse")) 
#---------------------combine and batch remove-------------
# load data---------
GSEname.list=c("GSE39582", "GSE13294", "GSE14333", "GSE17536", "GSE20916", "GSE2109",
               "GSE23878", "GSE33113", "GSE35896", "GSE37892", "KFSYSCC",  "GSE13067")
all.data=list()
exprs=NULL
col.annot=NULL
for(i in GSEname.list){
  load(paste0(PathToData,"Synapse/",i,".Affy.RData"))
  if (i=="KFSYSCC"){
    x=data$col.annot
    x$GSM=NA
    x$Dataset=x$dataset
    x$SampleId=gsub(".CEL","",x$SampleId)
    x=x[,c("GSM","Dataset","SampleId")]
    rownames(x)=x$SampleId
    colnames(data$data$rma)=gsub(".CEL","",colnames(data$data$rma))
    x=x[colnames(data$data$rma),]
    data$col.annot=x
    rm(x)
  }
  all.data[[i]]=data
  if (is.null(exprs)) {
    exprs <- data$data$rma
    col.annot <- data$col.annot[,1:3]
  } else {
    exprs <- cbind(exprs, data$data$rma[rownames(exprs),])
    col.annot<-rbind(col.annot,data$col.annot[,1:3])
  }
}

# row.annot----------
row.annot=all.data[[1]]$row.annot
# col.annot----------
CMS.dir=paste0(PathToData,"Synapse/","ColorectalCancerSubtypingConsortium/data/mergedPhenotype/cms_labels_public_all.txt")
cms=read.table(CMS.dir,sep="\t",header = T)
rownames(cms)=toupper(cms$sample)
rownames(cms)[which(cms$dataset=="gse33113")]=
  col.annot$GSM[match(cms$sample[which(cms$dataset=="gse33113")],col.annot$SampleId)]
samples=intersect(rownames(col.annot),rownames(cms))

exprs=exprs[,samples]
col.annot=col.annot[samples,]
col.annot=cbind(col.annot,cms[rownames(col.annot),-(1:2)])


clinical.dir=paste0(PathToData,"Synapse/","ColorectalCancerSubtypingConsortium/data/mergedPhenotype/clinical_molecular_public_all.txt")
clinical=read.table(clinical.dir,sep="\t",header = T)
rownames(clinical)=toupper(clinical$sample)
rownames(clinical)[which(clinical$sample%in%col.annot$SampleId)]=
  col.annot$GSM[match(clinical$sample[which(clinical$sample%in%col.annot$SampleId)],col.annot$SampleId)]

col.annot=cbind(col.annot,clinical[rownames(col.annot),-(1:2)])

# remove batch effect------
library(preprocessCore)
rma.quantile.normalized=normalize.quantiles(as.matrix(exprs))
colnames(rma.quantile.normalized)=colnames(exprs)
rownames(rma.quantile.normalized)=rownames(exprs)
library(sva)
rma.combat=ComBat(as.matrix(rma.quantile.normalized),batch = as.numeric(as.factor(col.annot$Dataset)))

# save data-------
data = list()
data$data$rma = exprs
data$data$rma.Quantile.comBat = data.frame(rma.combat)
data$col.annot = col.annot
data$row.annot = row.annot
data$init = T
save(data, file=paste(PathToData, "Synapse/Synapse.RData",sep=""))