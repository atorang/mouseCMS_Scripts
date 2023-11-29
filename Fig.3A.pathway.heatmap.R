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
# load data-----------
load(paste0(PathToData,"/organoids.RData")) 
data = iClassification.class.data.to.subset.cols(data, "origin", "SIP") 
data$data$raw=data$data$raw[,order(data$col.annot$genotype)]
data$data$quantile=data$data$quantile[,order(data$col.annot$genotype)]
data$col.annot=data$col.annot[order(data$col.annot$genotype),]

CMS=data.frame(readxl::read_excel(paste0(PathToData,"/First.CMS.Stratification.xlsx")))
CMS$v=CMS$dif/CMS$Max   #ratio of difference
CMS=CMS[CMS$Max>=0.5,]
CMS=CMS[CMS$SecMAX<0.2,]
data$col.annot$CMS=NA

data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS2"]]="CMS2"
data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS3"]]="CMS3"
data$col.annot$CMS[data$col.annot$genotype%in%CMS$samples[CMS$CMS=="CMS4"]]="CMS4"

data$data$raw=data$data$raw[,order(data$col.annot$CMS)]
data$data$quantile=data$data$quantile[,order(data$col.annot$CMS)]
data$col.annot=data$col.annot[order(data$col.annot$CMS),]
# define CMSs-------------------------
synapse.pw=gmtPathways(paste0(PathToData,"/Synapse.genesets.gmt"))
allpw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))

cms2=unique(c(unlist(synapse.pw["WNT_FLIER"]),
              unlist(allpw[c("BIOCARTA_WNT_PATHWAY")]),
              unlist(synapse.pw["MYC_TARGETS_ZELLER"]),
              unlist(allpw[c("PID_MYC_PATHWAY",
                             "BIOCARTA_SRCRPTP_PATHWAY",
                             "REACTOME_TRANSLATION")])
              ))

cms3=unique(unlist(allpw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                           "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
                           "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
                           "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
                           "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
                           "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
                           "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
                           "KEGG_LINOLEIC_ACID_METABOLISM")]))

cms4=unique(c(unlist(synapse.pw[c("MESENCH_LOBODA","EMT_CORE_GENES")]),
              unlist(allpw[c("KEGG_TGF_BETA_SIGNALING_PATHWAY",
                             "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
                             "GOBP_RESPONSE_TO_WOUNDING")]),
              unlist(synapse.pw["STROMAL_ESTIMATE"]),             
              unlist(allpw[c("REACTOME_INTEGRIN_SIGNALING",
                           "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS")])
            ))


pw=pw.cms=list(CMS2=cms2,CMS3=cms3,CMS4=cms4)

genes.path=c(unlist(synapse.pw["WNT_FLIER"]),
             unlist(allpw[c("BIOCARTA_WNT_PATHWAY")]),
             unlist(synapse.pw["MYC_TARGETS_ZELLER"]),
             unlist(allpw[c("PID_MYC_PATHWAY",
                            "BIOCARTA_SRCRPTP_PATHWAY",
                            "REACTOME_TRANSLATION")]),
             
             unlist(allpw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
                            "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
                            "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
                            "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
                            "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
                            "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
                            "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
                            "KEGG_LINOLEIC_ACID_METABOLISM")]),
             
             unlist(synapse.pw[c("MESENCH_LOBODA","EMT_CORE_GENES")]),
             unlist(allpw[c("KEGG_TGF_BETA_SIGNALING_PATHWAY",
               "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
               "GOBP_RESPONSE_TO_WOUNDING")]),
             unlist(synapse.pw["STROMAL_ESTIMATE"]),             
             unlist(allpw[c("REACTOME_INTEGRIN_SIGNALING",
               "REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS")])
             
)


# genes in pathways----------------
genes=unique(c(pw.cms$CMS2,pw.cms$CMS3,pw.cms$CMS4))
genes1=unlist(pw.cms)[duplicated(unlist(pw.cms))]
genes=setdiff(genes,genes1)

annot.row=unlist(pw.cms)[!duplicated(unlist(pw.cms))]
annot.row=data.frame(row.names = annot.row,CMS=substr(names(annot.row),1,4))
annot.row=data.frame(annot.row[genes,])
rownames(annot.row)=genes
colnames(annot.row)="CMS"
genes.path=data.frame(row.names = genes.path[!duplicated(genes.path)],path=names(genes.path)[!duplicated(genes.path)])
annot.row$pathway=gsub('[0-9]*',"",genes.path[rownames(annot.row),"path"])
table(tolower(annot.row$pathway))
table((annot.row$pathway))


annot.row$pathway=gsub("MESENCH_LOBODA","Mesenchymal",annot.row$pathway)
annot.row$pathway=gsub("STROMAL_ESTIMATE","Stromal infiltration",annot.row$pathway)
annot.row$pathway=gsub("WNT_FLIER","WNT targets",annot.row$pathway)
annot.row$pathway=gsub("BIOCARTA_WNT_PATHWAY","WNT targets",annot.row$pathway)
annot.row$pathway=gsub("MYC_TARGETS_ZELLER","MYC targets",annot.row$pathway)
annot.row$pathway=gsub("PID_MYC_PATHWAY","MYC targets",annot.row$pathway)
annot.row$pathway=gsub("EMT_CORE_GENES","EMT activation",annot.row$pathway)
annot.row$pathway=gsub("KEGG_TGF_BETA_SIGNALING_PATHWAY","TGFB activation",annot.row$pathway)
annot.row$pathway=gsub("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION","Matrix remodeling",annot.row$pathway)
annot.row$pathway=gsub("GOBP_RESPONSE_TO_WOUNDING","Wound response",annot.row$pathway)
annot.row$pathway=gsub("GOBP_ACTIVATION_OF_IMMUNE_RESPONSE","Immune response",annot.row$pathway)
annot.row$pathway=gsub("REACTOME_PD__SIGNALING","PD1 activation",annot.row$pathway)
annot.row$pathway=gsub("GOBP_COMPLEMENT_ACTIVATION","Complement activation",annot.row$pathway)

annot.row$pathway=gsub("BIOCARTA_SRCRPTP_PATHWAY","SRC",annot.row$pathway)
annot.row$pathway=gsub("KEGG_JAK_STAT_SIGNALING_PATHWAY","JAK STAT",annot.row$pathway)
annot.row$pathway=gsub("REACTOME_TRANSLATION","Transation ribosome",annot.row$pathway)
annot.row$pathway=gsub("REACTOME_INTEGRIN_SIGNALING","Integrin B3",annot.row$pathway)
annot.row$pathway=gsub("REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS","VEGF, VEGFR",annot.row$pathway)

annot.row$pathway=gsub("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM","Amino sugar, nucleotide sugar",annot.row$pathway)
annot.row$pathway=gsub("KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS","Pentose, glucuronate",annot.row$pathway)
annot.row$pathway=gsub("KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM","Fructose, mannose",annot.row$pathway)
annot.row$pathway=gsub("KEGG_GALACTOSE_METABOLISM","Galactose",annot.row$pathway)
annot.row$pathway=gsub("GOBP_GLUTAMINE_METABOLIC_PROCESS","Glutamine",annot.row$pathway)
annot.row$pathway=gsub("KEGG_GLUTATHIONE_METABOLISM","Glutathione",annot.row$pathway)
annot.row$pathway=gsub("KEGG_NITROGEN_METABOLISM","Nitrogen",annot.row$pathway)
annot.row$pathway=gsub("KEGG_GLYCEROPHOSPHOLIPID_METABOLISM","Glycerophospholipid",annot.row$pathway)
annot.row$pathway=gsub("KEGG_FATTY_ACID_METABOLISM","Fatty acid",annot.row$pathway)
annot.row$pathway=gsub("KEGG_STARCH_AND_SUCROSE_METABOLISM","Starch, sucrose",annot.row$pathway)
annot.row$pathway=gsub("KEGG_ARACHIDONIC_ACID_METABOLISM","Arachidonic acid",annot.row$pathway)
annot.row$pathway=gsub("KEGG_LINOLEIC_ACID_METABOLISM","Linoleic acid",annot.row$pathway)
annot.row$pathway=gsub("KEGG_TYROSINE_METABOLISM","Tyrosine",annot.row$pathway)
table(annot.row$pathway)
x=data.frame(table((annot.row$pathway)))
annot.row=annot.row[which(annot.row$pathway%in%x$Var1),]
# orthologs---------
exprs.m=as.matrix(data$data$quantile)
mmm=apply(exprs.m,1,function(x){sum(x>0)})
exprs.m=exprs.m[mmm>=3,]
exprs=exprs.m

load(paste0(PathToData,"/orth.v105.RData"))

rownames(exprs) = orth$Gene.name.1[match(rownames(exprs), orth$Gene.stable.ID)]
exprs = exprs[!is.na(rownames(exprs)),]
exprs = exprs[apply(exprs,1,function(x) sum(is.na(x))==0),]
exprs=exprs[!duplicated(rownames(exprs)),]

colnames(exprs)=gsub("F1.|.SIP","",colnames(exprs))
# heatmaps parameters-----------
genes=intersect(rownames(annot.row),rownames(exprs))
exprs.z=exprs[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
}
annot.row=annot.row[genes,]
cumsum (table(annot.row$pathway))

library(ComplexHeatmap)
col_fun =circlize:: colorRamp2(seq(-1.5,1.5,3/6), colorRampPalette(c("#1483E4","white","#FF8C00"))(7) )

loc=table(annot.row$pathway)
loc=loc[unique(annot.row$pathway)]
loc2=loc
loc2[1]=round(loc2[1]/2)
for(i in 2:length(loc)){
  loc2[i]=round(sum(loc[1:(i-1)])+loc[i]/2)
}
loc=loc2

# plot----------------
Heatmap(exprs.z[rownames(annot.row),],
        row_split =factor(annot.row$pathway,
                          levels = unique(annot.row$pathway)),
        row_dend_reorder =T,
        column_dend_reorder =T,
        cluster_row_slices = F,
        cluster_rows = T,
        row_km_repeats=100, 
        row_title="  ",
        row_title_gp = gpar(fontsize=10),
        left_annotation = rowAnnotation(foo=anno_mark(at = loc, labels = names(loc),
                                                      link_width = unit(7, "mm"),
                                                      side="left",labels_gp = gpar(fontsize=7),
                                                      padding = 0.2),
                                        CMS=annot.row[,1],col=list(CMS=cms.cols),
                                        show_legend = c(F)
                                        ),
        column_split =factor(data$col.annot$genotype,
                             levels = unique(data$col.annot$genotype)),
        cluster_column_slices = F,
        cluster_columns = T,
        column_title_side = "bottom",
        column_title_gp = gpar(fontsize=8),
        column_title = "%s",
        column_title_rot=90,
        show_row_names = F,
        show_column_names = F,
        column_labels = data$col.annot$genotype,
        column_names_gp = gpar(fontsize=4),
        show_parent_dend_line = FALSE,
        row_title_rot = 0,
        col=col_fun,
        row_gap = unit(0.3, "mm"), 
        column_gap = unit(0.5, "mm"),
        border = F,
        top_annotation = HeatmapAnnotation(CMS=data$col.annot$CMS,
                                           col = list(CMS=cms.cols),
                                           simple_anno_size = unit(0.6, "cm"),
                                           show_legend = c(F),
                                           annotation_name_side = "left",
                                           gap = unit(0.1, "cm")),
        name="Z-Score", show_heatmap_legend = T, 
        show_row_dend = F, show_column_dend = F)



