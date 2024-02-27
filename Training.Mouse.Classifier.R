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

# CMS colors---------
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74"); names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# define function trainer---------
trainer<-function(exprs,data,train.repeats,genelist.length,n.fold, W=1,weight,bootstraps,alpha){
  Total=list()
  C=T
  for (w in W) {
    if(w==1&C){b=0;C=F}else{b=1}
    genelists = matrix(nrow=genelist.length , ncol=train.repeats)
    
    design = cbind(data$col.annot$CMS,data$col.annot$CMS,data$col.annot$CMS); colnames(design) = c( "CMS2", "CMS3", "CMS4")
    design[design[,1]!="CMS2",1]=0
    design[design[,2]!="CMS3",2]=0
    design[design[,3]!="CMS4",3]=0
    design[design!=0]=1
    design = apply(design,2,as.numeric)
    design.cms2 = cbind(Intercept=1,design[,"CMS2"])
    design.cms3 = cbind(Intercept=1,design[,"CMS3"])
    design.cms4 = cbind(Intercept=1,design[,"CMS4"])
    sample.results =  data.frame(cbind(class=rep(0,nrow(data$col.annot)), "predictions"=0,"errors"=0,"error.rate"=0), 
                                 row.names=rownames(data$col.annot))
    sample.results$class = data$col.annot$CMS
    sample.results = cbind(sample.results, design)
    sample.results$CMS2=sample.results$CMS3=sample.results$CMS4=0
    results = vector(length=train.repeats)
    
    for(i in seq(train.repeats)){
      exprs.train = exprs[,bootstraps[i,]]
      test.samples = setdiff(1:ncol(exprs), bootstraps[i,])
      exprs.test = exprs[,test.samples]
      # select genes by elastic net
      model=tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                               family="multinomial", thresh = 1e-07,
                               type.multinomial="grouped", nfolds=5,alpha=alpha,
                               # type.measure = "class" #,trace.it=1,
                               weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) tryCatch(cv.glmnet(t(exprs.train),data$col.annot$CMS[bootstraps[i,]],
                                           family="multinomial", thresh = 1e-07,
                                           type.multinomial="grouped", nfolds=5,alpha=alpha,
                                           # type.measure = "class" #,trace.it=1,
                                           weights = weight[bootstraps[i,]]
      ) ,
      error=function(e) model))))))))
      
      conf=data.frame(confusion.glmnet (model,newx=t(exprs.test),newy=data$col.annot$CMS[test.samples], family="multinomial"))
      rownames(conf)=paste0(conf[,1],conf[,2])
      correct=conf[c("CMS2CMS2","CMS3CMS3","CMS4CMS4"),"Freq"]
      correct[is.na(correct)]=0
      results[i] =sum(correct)/sum(conf$Freq)
      coef=predict(model, type = "coef")
      probes=coef$CMS2@i[-1]
      prediction=predict(model, newx=t(exprs.test), interval ="prediction",type="class")
      genelists[1:length(probes),i] = probes
      
      
      # update sample.results
      for(j in 1:ncol(exprs.test)){
        sample = colnames(exprs.test)[j]
        if(prediction[j]!=data$col.annot$CMS[test.samples[j]]){ # false prediction
          sample.results[sample,"errors"] = as.numeric(sample.results[sample,"errors"]) + 1
        }
        sample.results[sample,prediction[j]] = sample.results[sample,prediction[j]] + 1
        sample.results[sample,"predictions"] = sample.results[sample,"predictions"] + 1
        sample.results[sample,"error.rate"] = sample.results[sample,"errors"]/sample.results[sample,"predictions"]
      }
      print(paste("i=",i,"   ","w=",w ,"    ","b=",b,sep=""))
    }
    sample.results = sample.results[order(sample.results$error.rate),]
    Total[[paste0(b,"Weight",w)]]$sample.results=sample.results
    Total[[paste0(b,"Weight",w)]]$results=results
    Total[[paste0(b,"Weight",w)]]$genelists=genelists
  }
  return(Total)
}
# define alpha.optimizer------
#if optimization needed run: alpha=T
alpha.optimizer<-function(alpha=F,file.Name='alpha',exprs,bootstraps,weight,data,alpha.values=seq(0,1,0.1)){
  if(alpha){
    i=1
    exprs.train = exprs[,bootstraps[i,]]
    test.samples = setdiff(1:ncol(exprs), bootstraps[i,])
    exprs.test = exprs[,test.samples]
    # select genes by elastic net
    pdf(paste0(file.Name,'.pdf'))
    par(mfrow=c(2,2))
    for (j in alpha.values) {
      model=cv.glmnet(t(exprs.train),c(data$col.annot$CMS[bootstraps[i,]]),
                      family="multinomial", thresh = 1e-07,
                      type.multinomial="grouped", nfolds=5,alpha=j,
                      trace.it=0,weights = weight
      )
      plot(model,main=paste0(j))
      print(j)
    }
    dev.off()
  }
}
# load data----------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 

CMS=data.frame(readxl::read_excel(paste0(PathToData,"/First.CMS.Stratification.xlsx")))
CMS$v=CMS$dif/CMS$Max   #ratio of difference
CMS=CMS[CMS$Max>=0.5,]
CMS=CMS[CMS$SecMAX<0.2,]
data$col.annot$CMS=NA

data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS2"]]="CMS2"
data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS3"]]="CMS3"
data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS4"]]="CMS4"

table(data$col.annot$CMS)
data = iClassification.class.data.to.subset.cols(data, "CMS", c("CMS2","CMS3","CMS4")) 

# filteration--------
exprs=data$data$quantile
S=apply(exprs, 1, sd)
exprs=exprs[S>0.5,]
S=apply(exprs,1,function(x){sum(x>=5)})
exprs=exprs[S>=3,]
exprs.Train=exprs
rm(S,exprs)
data.Train=data
cms=data.Train$col.annot$CMS
names(cms)=rownames(data.Train$col.annot)
# normalization-----------------
exprs.tz=exprs.Train
exprs.tz=preprocessCore::normalize.quantiles(as.matrix(exprs.tz))
colnames(exprs.tz)=colnames(exprs.Train)
rownames(exprs.tz)=rownames(exprs.Train)
# setting hyperparameters and alpha optimization------------
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
library(glmnet)
exprs=exprs.tz
data=list()
data$col.annot=data.frame(row.names = names(cms),CMS=cms)
train.repeats=1000
genelist.length=200
n.fold=5
W=c(1)
weight=rep(1,length(data$col.annot$CMS))

bootstraps = t(apply(matrix(1:train.repeats, ncol=1),1, 
                     function(x) sample(1:nrow(data$col.annot),
                                        (nrow(data$col.annot)/n.fold)*(n.fold-1))))  # resampling without replacement

#optimize Alpha (the coef for Lasso and Ridge) one time run
alpha=F   #if optimization needed run: alpha=T
alpha.optimizer(alpha,file.Name='alpha1',exprs,bootstraps,weight,data)
alpha=j=0.95 #manually change to the best value by generated figure above

# train model and select genes-------
start.time = proc.time()[3]
Total=trainer(exprs,data,train.repeats,genelist.length,n.fold, W=1,weight,bootstraps,alpha)
cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))

save.image("geneSelection1.RData")

# result and gene selection----------
res=matrix(nrow=train.repeats,ncol=length(W))
Av=rep(0,length(W))
for (i in 1:length(W)) {
  res[,i]=Total[[i]]$results
  Av[i]=1-sum(Total[[i]]$sample.results$errors)/sum(Total[[i]]$sample.results$predictions)
}
names(Av)<-names(Total); Av   #99.3%
 
gene.freqs = table(Total[[1]]$genelists[!is.na(Total[[1]]$genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
length(gene.freqs)
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] #315 genes selected

# normalize after gene selection-----------------
genes=names(gene.freqs) 
genes=genes[gene.freqs>=100] #78 genes
exprs.tz=exprs.Train[genes,]
exprs.tz=preprocessCore::normalize.quantiles(as.matrix(exprs.tz))
colnames(exprs.tz)=colnames(exprs.Train)
rownames(exprs.tz)=genes
#---------------------step2---------------
# setting hyperparameters and alpha optimization------------
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
library(glmnet)
exprs=exprs.tz
data=list()
data$col.annot=data.frame(row.names = names(cms),CMS=cms)
train.repeats=1000
genelist.length=200
n.fold=5
W=c(1)
weight=rep(1,length(data$col.annot$CMS))

bootstraps = t(apply(matrix(1:train.repeats, ncol=1),1, 
                     function(x) sample(1:nrow(data$col.annot),
                                        (nrow(data$col.annot)/n.fold)*(n.fold-1))))  # resampling without replacement

#optimize Alpha (the coef for Lasso and Ridge) one time run
alpha=F   #if optimization needed run: alpha=T
alpha.optimizer(alpha,file.Name='alpha2',exprs,bootstraps,weight,data)
alpha=j=0.8 #manually change to the best value by generated figure above

# train model and select genes-------
start.time = proc.time()[3]
Total=trainer(exprs,data,train.repeats,genelist.length,n.fold, W=1,weight,bootstraps,alpha)
cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))
save.image("geneSelection2.RData")

# result and gene selection----------
res=matrix(nrow=train.repeats,ncol=length(W))
Av=rep(0,length(W))
for (i in 1:length(W)) {
  res[,i]=Total[[i]]$results
  Av[i]=1-sum(Total[[i]]$sample.results$errors)/sum(Total[[i]]$sample.results$predictions)
}
names(Av)<-names(Total); Av      #100%
 
gene.freqs = table(Total[[1]]$genelists[!is.na(Total[[1]]$genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
length(gene.freqs)
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] #67 genes selected

# normalize after gene selection-----------------
genes=names(gene.freqs) 
genes=genes[gene.freqs>100] #54
pheatmap(exprs[genes,])

exprs.tz=exprs.Train[genes,]

exprs.tz=preprocessCore::normalize.quantiles(as.matrix(exprs.tz))
colnames(exprs.tz)=colnames(exprs.Train)
rownames(exprs.tz)=genes
#---------------------step3---------------
# setting hyperparameters and alpha optimization------------
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
library(glmnet)
exprs=exprs.tz
data=list()
data$col.annot=data.frame(row.names = names(cms),CMS=cms)
train.repeats=1000
genelist.length=200
n.fold=5
W=c(1)
weight=rep(1,length(data$col.annot$CMS))

bootstraps = t(apply(matrix(1:train.repeats, ncol=1),1, 
                     function(x) sample(1:nrow(data$col.annot),
                                        (nrow(data$col.annot)/n.fold)*(n.fold-1))))  # resampling without replacement

#optimize Alpha (the coef for Lasso and Ridge) one time run
alpha=F   #if optimization needed run: alpha=T
alpha.optimizer(alpha,file.Name='alpha3',exprs,bootstraps,weight,data)
alpha=j=0.8 #manually change to the best value by generated figure above

# train model and select genes-------
start.time = proc.time()[3]
Total=trainer(exprs,data,train.repeats,genelist.length,n.fold, W=1,weight,bootstraps,alpha)
cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))
save.image("geneSelection3.RData")

# result and gene selection----------
res=matrix(nrow=train.repeats,ncol=length(W))
Av=rep(0,length(W))
for (i in 1:length(W)) {
  res[,i]=Total[[i]]$results
  Av[i]=1-sum(Total[[i]]$sample.results$errors)/sum(Total[[i]]$sample.results$predictions)
}
names(Av)<-names(Total); Av         #100%

gene.freqs = table(Total[[1]]$genelists[!is.na(Total[[1]]$genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
length(gene.freqs)
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] #54 genes selected

# normalize after gene selection-----------------
genes=names(gene.freqs) 
genes=genes[gene.freqs>=100] #48
pheatmap(exprs[genes,])

exprs.tz=exprs.Train[genes,]
exprs.tz=preprocessCore::normalize.quantiles(as.matrix(exprs.tz))
colnames(exprs.tz)=colnames(exprs.Train)
rownames(exprs.tz)=genes
#---------------------------------------train final classifier---------------
# Final classifier------------
load("geneSelection2.RData")
w=1
alpha=j=0.8 
exprs=exprs.tz

pheatmap(exprs)
model=cv.glmnet(as.matrix(t(exprs)),c(data$col.annot$CMS),
                family="multinomial", thresh = 1e-07,
                type.multinomial="grouped", nfolds=5,alpha=alpha,
                weights = weight
) 

plot(model)
confusion.glmnet (model,newx=t(exprs),newy=c(data$col.annot$CMS), family="multinomial")
prediction=predict(model, newx=t(exprs), interval ="prediction",type="response")
sort(apply(prediction, 1, max))
coef=predict(model, type = "coef")
probes=coef$CMS2@i[-1]
genelists = probes
genelists=rownames(exprs)[genelists]
genes=rownames(exprs)
data.classifier=data
rm(list=as.character(setdiff(ls(),c("model","genes","exprs.Train","exprs.tz","data.classifier","alpha"))))
save.image("Mouse.Classifier.RData")


