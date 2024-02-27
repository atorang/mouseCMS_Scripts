# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load count data from fastq-----
Reverse_stranded= read.table(paste0(PathToData,"/rawCounts1/output_Reverse_stranded.txt"),stringsAsFactors = F)
names= read.table(paste0(PathToData,"/rawCounts1/Colnames.txt"),stringsAsFactors = F)

#Gene annotation: download ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
gtf= read.table(paste0(PathToData,"/gencode.v44.primary_assembly.annotation.gtf"), sep="\t",stringsAsFactors = F)
gtf=gtf[as.character(gtf$V3)=="gene",]


#rownames
rownames(Reverse_stranded)=Reverse_stranded$V1


#Remove column 1
Reverse_stranded=Reverse_stranded[,-1]


#Remove 4 rows
Reverse_stranded=Reverse_stranded[-c(1:4),]


#colnames
colnames(Reverse_stranded)=gsub("\\_.*","",gsub("ReadsPerGene.out.tab","",gsub("trim.|-","",t(names))))

rownames(gtf)=rownames(Reverse_stranded)


# define raw---------
raw=Reverse_stranded
batch=rep(1,ncol(raw))
names(batch)=colnames(raw)

raw=raw[!duplicated(gsub("\\..*","",rownames(raw))),]
rownames(raw)=gsub("\\..*","",rownames(raw))

gtf=gtf[!duplicated(gsub("\\..*","",rownames(gtf))),]
rownames(gtf)=gsub("\\..*","",rownames(gtf))

# Normalization--------
library(preprocessCore)
quantile=normalize.quantiles(as.matrix(log2(raw+1)))
rownames(quantile)=row.names(raw)
colnames(quantile)=colnames(raw)

# row.annot----------
row.annot=gtf

#gene names
gtf1=gtf[setdiff(grep("gene_name ",gtf$V9 , fixed=TRUE),grep("gene_name ENSG00",gtf$V9 , fixed=TRUE)),]
names=gtf1$V9
names=gsub(".*\\gene_name ","",names)
names=gsub("\\;.*","",names)
genes=rownames(gtf1)
row.annot$gene=NA
row.annot[genes,"gene"]=names


#hgnc Id
gtf1=gtf[grepl("hgnc_id HGNC:",gtf$V9 , fixed=TRUE),]
names=gtf1$V9
names=gsub(".*\\hgnc_id HGNC:","",names)
names=gsub("\\;.*","",names)
genes=rownames(gtf1)
row.annot$hgncId=NA
row.annot[genes,"hgncId"]=names


#ensembl Id
gtf1=gtf[grepl("gene_id ",gtf$V9 , fixed=TRUE),]
names=gtf1$V9
names=gsub(".*\\gene_id ","",names)
names=gsub("\\..*","",names)
genes=rownames(gtf1)
row.annot$ensemblId=NA
row.annot[genes,"ensemblId"]=names


#transcript Id
gtf1=gtf[grepl("gene_id ",gtf$V9 , fixed=TRUE),]
names=gtf1$V9
names=gsub(".*\\gene_id ","",names)
names=gsub("\\;.*","",names)
genes=rownames(gtf1)
row.annot$transcriptId=NA
row.annot[genes,"transcriptId"]=names

#gene type
gtf1=gtf[grepl("gene_type ",gtf$V9 , fixed=TRUE),]
names=gtf1$V9
names=gsub(".*\\gene_type ","",names)
names=gsub("\\;.*","",names)
genes=rownames(gtf1)
row.annot$geneType=NA
row.annot[genes,"geneType"]=names

#havana gene id
gtf1=gtf[grepl("havana_gene ",gtf$V9 , fixed=TRUE),]
names=gtf1$V9
names=gsub(".*\\havana_gene ","",names)
names=gsub("\\;.*","",names)
genes=rownames(gtf1)
row.annot$havanaId =NA
row.annot[genes,"havanaId"]=names

#entrez Id
library(annotate)
library('org.Hs.eg.db')
columns(org.Hs.eg.db)
row.annot1=select(org.Hs.eg.db,rownames(raw),keytype='ENSEMBL',c("ENTREZID",'SYMBOL',"ENSEMBL"))
row.annot1<-row.annot1[!duplicated(row.annot1$ENSEMBL),]
rownames(row.annot1)<-row.annot1$ENSEMBL
sum(is.na(row.annot1$ENTREZID))
sum(rownames(row.annot1)!=rownames(row.annot))
row.annot$entrezId=row.annot1$ENTREZID
rm(row.annot1,gtf1)



colnames(row.annot)[1:9]=c("chromosome","source","a","start","end","a2","strand","a3","description")
row.annot=row.annot[,order(colnames(row.annot))]
row.annot=row.annot[,-c(1:3)]
row.annot=row.annot[,c("gene","entrezId","ensemblId","hgncId","havanaId","transcriptId","geneType","description","source",
                       "chromosome","start","end","strand")]

# col.annot----
library(readxl)
col.annot=data.frame(row.names = colnames(raw),samplId=colnames(raw),
                     Lane=batch,species="human",sequencer="",stringsAsFactors = F,
                     lowDept=0,sampleType="Sanger Cell Line",dateId="2023")

col.annot$CellLine=gsub("n.*","",colnames(raw))
col.annot$treatment=gsub(".*n1|.*n2|.*n3|e|plus","",colnames(raw))
col.annot$duplicates=gsub("n1|n2|n3|e","",colnames(raw))

# save data-----
data=list()
data$data$raw=data.frame(raw)
data$data$quantile=data.frame(data.frame(quantile))
data$col.annot=data.frame(col.annot)
data$row.annot=data.frame(row.annot)
data$init=TRUE
save(data, file=paste0(PathToData,"/CellLine.GSE248558.RData"))
