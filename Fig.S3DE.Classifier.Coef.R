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
iClassification.class.data.to.subset.rows = function (data, row.name, values){
  .check.data.init(data, "iClassification.class.data.to.subset.rows")
  data.subset = match(data$row.annot[,row.name], values, nomatch=0)
  if(length(which(data.subset>0))==1){ #cast annotation and values as matrices in case number of rows = 1 (otherwise the matrices become vectors)
    if(ncol(data$col.annot)==1){
      stop(paste("[iClassification.class.data.to.subset.rows]: can't subset to a single row if data has only 1 column!", sep=""), call. = TRUE)	
    }
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = matrix(data$data[[i]][data.subset>0,], nrow=1) 					
    }		
    # subset the annotation		
    data$row.annot = as.matrix(data$row.annot[data.subset>0,], nrow=1)				
  }else{
    # subset the data		
    for(i in 1:length(data$data)){
      data$data[[i]] = data$data[[i]][data.subset>0,] 					
    }		
    # subset the annotation		
    data$row.annot = data$row.annot[data.subset>0,]				
  }
  
  return(data)
}
# CMS colors---------
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74"); names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# load classifier-----------
# defined by script Training.Mouse.Classifier.R
library(glmnet)
load(paste0(PathToData,"/Mouse.Classifier.RData"))
#---------------------------SIP---------------
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 
data=iClassification.class.data.to.subset.rows(data,"ensemblId",genes)

exprs.test=data$data$quantile
exprs.test=exprs.test[genes,]
for (i in 1:ncol(exprs.test)) {
  comb=preprocessCore::normalize.quantiles(as.matrix(data.frame(exprs.tz,exprs.test[,i])))
  exprs.test[,i]=comb[,ncol(comb)]
}

prediction=predict(model, newx=t(exprs.test), interval ="prediction",type="response")
cls=predict(model, newx=t(exprs.test), interval ="prediction",type="class")
Max=apply(prediction, 1, max)
Predicted=cls
Predicted[Max<0.5]<-NA
results=data.frame(prediction,Nearest=cls,Predicted=Predicted,genotype=data$col.annot$genotype)


# model features and z-score-----------
x=coef(model)
x=data.frame(as.matrix(x$CMS2),as.matrix(x$CMS3),as.matrix(x$CMS4))[-1,]
S=apply(x, 1,sd)
x=x[S>0,]
genes.nonzero=rownames(x)
exprs.z=data$data$quantile[genes.nonzero,]
rownames(exprs.z)=data$row.annot[rownames(exprs.z),"gene"]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
}

# plot----------------
library(ComplexHeatmap)
col_fun =circlize:: colorRamp2(seq(-1.5,1.5,3/6), colorRampPalette(c("#1483E4","white","#FF8C00"))(7) )

h1={Heatmap(exprs.z,
            cluster_rows = T,
            row_km_repeats=100, 
            cluster_columns = T,
            column_split =factor(data$col.annot[,c(9)],
                                 levels = unique(data$col.annot$genotype)),
)}
o1 = row_order(h1)

h1.alt={Heatmap(exprs.z[o1,order(results$X1.1)],
                cluster_rows = F,
                row_km_repeats=100, 
                row_title="  ",
                row_title_gp = gpar(fontsize=8),
                row_names_gp = gpar(fontsize=8),
                cluster_columns =T,
                cluster_column_slices = F,
                show_column_dend = F,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize=8),
                column_title_rot=90,
                show_parent_dend_line = F,
                show_column_names = F,
                column_labels = data$col.annot$genotype,
                column_split =factor(results$X1.1[order(results$X1.1)],
                                     levels = unique(results$X1.1[order(results$X1.1)])),
                col=col_fun,
                row_gap = unit(0.3, "mm"), 
                column_gap = unit(1, "mm"),
                border = F,
                top_annotation = HeatmapAnnotation(
                  CMS=results$X1.1[order(results$X1.1)],
                  col = list(CMS=cms.cols),
                  simple_anno_size = unit(0.6, "cm"),
                  show_legend = c(F),
                  annotation_name_side = "left",
                  gap = unit(0.1, "cm")),
                name="Z-Score", show_heatmap_legend = F, 
                show_row_dend = F)}
draw(h1.alt,padding = unit(c(10, 2, 2, 2), "mm"))
#---------------------------CO----------------
# test classifier on all samples---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "CO") 
data=iClassification.class.data.to.subset.rows(data,"ensemblId",genes)

exprs.test=data$data$quantile
exprs.test=exprs.test[genes,]
for (i in 1:ncol(exprs.test)) {
  comb=preprocessCore::normalize.quantiles(as.matrix(data.frame(exprs.tz,exprs.test[,i])))
  exprs.test[,i]=comb[,ncol(comb)]
}

prediction=predict(model, newx=t(exprs.test), interval ="prediction",type="response")
cls=predict(model, newx=t(exprs.test), interval ="prediction",type="class")
Max=apply(prediction, 1, max)
Predicted=cls
Predicted[Max<0.5]<-NA
results=data.frame(prediction,Nearest=cls,Predicted=Predicted,genotype=data$col.annot$genotype)


# model features and z-score-----------
exprs.z=data$data$quantile[genes.nonzero,]
rownames(exprs.z)=data$row.annot[rownames(exprs.z),"gene"]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
}

# plot----------------
h2.alt={Heatmap(exprs.z[o1,order(results$X1.1)],
                cluster_rows = F,
                row_km_repeats=100, 
                row_title="  ",
                row_title_gp = gpar(fontsize=8),
                row_names_gp = gpar(fontsize=8),
                cluster_columns =T,
                cluster_column_slices = F,
                show_column_dend = F,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize=8),
                # column_title = "%s",
                column_title_rot=90,
                show_parent_dend_line = F,
                show_column_names = F,
                column_labels = data$col.annot$genotype,
                column_split =factor(results$X1.1[order(results$X1.1)],
                                     levels = unique(results$X1.1[order(results$X1.1)])),
                col=col_fun,
                row_gap = unit(0.3, "mm"), 
                column_gap = unit(1, "mm"),
                border = F,
                top_annotation = HeatmapAnnotation(
                  CMS=results$X1.1[order(results$X1.1)],
                  col = list(CMS=cms.cols),
                  simple_anno_size = unit(0.6, "cm"),
                  show_legend = c(F),
                  annotation_name_side = "left",
                  gap = unit(0.1, "cm")),
                name="Z-Score", show_heatmap_legend = F, 
                show_row_dend = F)}
draw(h2.alt,padding = unit(c(10, 2, 2, 2), "mm"))



