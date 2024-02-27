# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load count data from fastq-----
unstranded= read.table(paste0(PathToData,"/output_unstranded.txt"),stringsAsFactors = F)
names= read.table(paste0(PathToData,"/Colnames.txt"),stringsAsFactors = F)


#Gene annotation:download: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
gtf= read.table(paste0(PathToData,"/gencode.vM25.primary_assembly.annotation.gtf"), sep="\t",stringsAsFactors = F)
gtf=gtf[as.character(gtf$V3)=="gene",]


#rownames
rownames(unstranded)=unstranded$V1


#Remove column 1
unstranded=unstranded[,-1]


#Remove 4 rows
unstranded=unstranded[-c(1:4),]


#colnames
colnames(unstranded)=gsub("\\_.*","",gsub(".fastq.gzReadsPerGene.out.tab","",gsub("trim.|-","",t(names))))


rownames(gtf)=rownames(unstranded)

# define raw---------
raw=data.frame(unstranded)
colnames(raw)=gsub("X","",colnames(raw))

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
gtf1=gtf[grepl("gene_name ",gtf$V9 , fixed=TRUE),]
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


colnames(row.annot)[1:9]=c("chromosome","source","a","start","end","a2","strand","a3","description")
row.annot=row.annot[,order(colnames(row.annot))]
row.annot=row.annot[,-c(1:3)]
row.annot=row.annot[,c("gene","ensemblId","hgncId","havanaId","transcriptId","geneType","description","source",
                       "chromosome","start","end","strand")]


# col.annot----
library(readxl)
annot=data.frame(read_xlsx(paste0(PathToData,"/annot.organoids.xlsx")))
rownames(annot)=annot$FileName
col.annot=annot[colname(raw),]

colnames(raw)=colnames(quantile)=rownames(col.annot)=col.annot$sampleName
# save data-----
data=list()
data$data$raw=data.frame(raw)
data$data$quantile=data.frame(data.frame(quantile))
data$col.annot=data.frame(col.annot)
data$row.annot=data.frame(row.annot)
data$init=TRUE
save(data, file=paste0(PathToData,"/organoids.RData"))
