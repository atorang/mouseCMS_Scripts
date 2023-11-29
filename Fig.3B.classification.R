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

# load classifier-----------
if (!require("mouseCMS", quietly = TRUE)){
  devtools::install_github("atorang/mouseCMS")
}
 
library(mouseCMS)
# ---------------SIP------------
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 
# ordering the samples---------
Ord=c("A","AP","AP2","AK","AKP","AKP2","AKPS","AKPN","APN","K","KP","BP","BPN","BPNA2","KPN","KPNA2")
Ord=intersect(Ord,data$col.annot$genotype)
col.annot=data$col.annot
j=1
for (i in Ord) {
  x=which(data$col.annot$genotype==i)
  col.annot[j:(j+length(x)-1),]=data$col.annot[x,]
  j=j+length(x)
}
rownames(col.annot)=col.annot$sampleName
data$col.annot=col.annot
data$data$raw=data$data$raw[,rownames(data$col.annot)]
data$data$quantile=data$data$quantile[,rownames(data$col.annot)]
rm(i,j,x,col.annot)

# classification----------
exprs.test=data$data$quantile
results=mouseCMSclassifier(exprs.test, perform.log2=FALSE)
results=data.frame(results,genotype=data$col.annot$genotype)

res.number=table(paste0(data$col.annot$genotype),results$predictedCMS,useNA="always")
res.number=res.number[Ord,]
colnames(res.number)[4]="Unclassified"
res.prop=t(apply(res.number, 1, proportions))*100
# plot---------
col_fun = circlize::colorRamp2(c(0, 100), c("white","darkblue"))

h1=Heatmap(res.prop, name = "Percentage", col = col_fun, rect_gp = gpar(type = "none"), 
           width = ncol(res.prop)*unit(10, "mm"), 
           height = nrow(res.prop)*unit(10, "mm"),
           cell_fun = function(jh, ih, x, y, width, height, fill) {
             grid.rect(x = x, y = y, width = width, height = height, 
                       gp = gpar(col = "grey", fill = NA))
             if(ih > 0&jh>0) {
               grid.circle(x = x, y = y, r = abs(res.number[ih, jh])/75 , 
                           gp = gpar(fill = col_fun(res.prop[ih, jh]), col = NA))
               }},
           show_heatmap_legend = FALSE,
           cluster_rows = FALSE, cluster_columns = FALSE)
lgd=Legend(title = "percentage",col_fun = col_fun,title_position = "leftcenter-rot")
draw(h1,heatmap_legend_list=lgd)


# --------------CO------------
# load data---------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "CO") 
# ordering the samples---------
Ord=c("A","AP","AP2","AK","AKP","AKP2","AKPS","AKPN","APN","K","KP","BP","BPN","BPNA2","KPN","KPNA2")
Ord=intersect(Ord,data$col.annot$genotype)
col.annot=data$col.annot
j=1
for (i in Ord) {
  x=which(data$col.annot$genotype==i)
  col.annot[j:(j+length(x)-1),]=data$col.annot[x,]
  j=j+length(x)
}
rownames(col.annot)=col.annot$sampleName
data$col.annot=col.annot
data$data$raw=data$data$raw[,rownames(data$col.annot)]
data$data$quantile=data$data$quantile[,rownames(data$col.annot)]
rm(i,j,x,col.annot)

# classification----------
exprs.test=data$data$quantile
results=mouseCMSclassifier(exprs.test, perform.log2=FALSE)
results=data.frame(results,genotype=data$col.annot$genotype)

res.number=table(paste0(data$col.annot$genotype),results$predictedCMS,useNA="always")
res.number=res.number[Ord,]
colnames(res.number)[4]="Unclassified"
res.prop=t(apply(res.number, 1, proportions))*100
# plot---------
col_fun = circlize::colorRamp2(c(0, 100), c("white","darkblue"))

h1=Heatmap(res.prop, name = "Percentage", col = col_fun, rect_gp = gpar(type = "none"), 
           width = ncol(res.prop)*unit(10, "mm"), 
           height = nrow(res.prop)*unit(10, "mm"),
           cell_fun = function(jh, ih, x, y, width, height, fill) {
             grid.rect(x = x, y = y, width = width, height = height, 
                       gp = gpar(col = "grey", fill = NA))
             if(ih > 0&jh>0) {
               grid.circle(x = x, y = y, r = abs(res.number[ih, jh])/75 , 
                           gp = gpar(fill = col_fun(res.prop[ih, jh]), col = NA))
             }},
           show_heatmap_legend = FALSE,
           cluster_rows = FALSE, cluster_columns = FALSE)
lgd=Legend(title = "percentage",col_fun = col_fun,title_position = "leftcenter-rot")
draw(h1,heatmap_legend_list=lgd)

