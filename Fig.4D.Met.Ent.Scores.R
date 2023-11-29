# ENTER YOUR DIRECTORY PATH---------
PathToData="ENTER YOUR DIRECTORY PATH"
# load data---------
load(paste0(PathToData,"/Synapse.High.RData"))
exprs=data$data$normal
exprs=exprs[!duplicated(data$row.annot$gene),]
rownames(exprs)=data$row.annot[rownames(exprs),"gene"]

cms=data$col.annot$CMS
names(cms)=colnames(exprs)

S=apply(exprs, 1, sd)
exprs=exprs[S!=0,]
S=apply(exprs, 1, max)
exprs=exprs[S>4,]

# CMS colors---------
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74"); names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# GSVA--------------------
library(GSVA)
library(fgsea)
pw=gmtPathways(paste0(PathToData,"/msigdb.v7.4.symbols.gmt"))
sigs.entero=names(pw)[grep("ENTEROCYTE",names(pw))]
pw=pw[c("KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM",
        "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
        "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM", "KEGG_GALACTOSE_METABOLISM",
        "GOBP_GLUTAMINE_METABOLIC_PROCESS","KEGG_GLUTATHIONE_METABOLISM",
        "KEGG_NITROGEN_METABOLISM","KEGG_GLYCEROPHOSPHOLIPID_METABOLISM",
        "KEGG_FATTY_ACID_METABOLISM","KEGG_STARCH_AND_SUCROSE_METABOLISM",
        "KEGG_TYROSINE_METABOLISM","KEGG_ARACHIDONIC_ACID_METABOLISM",
        "KEGG_LINOLEIC_ACID_METABOLISM",sigs.entero)]
sigs.entero=c("Duodenal early immature enterocyte",
              "Duodenal late immature enterocyte",
              "Duodenal mature enterocyte",
              "SI_24W_C3 enterocyte progenitor type1",
              "SI_24W_C4 enterocyte progenitor type2",
              "Colon_24W_C10 enterocyte")
names(pw)=c("Amino sugar, nucleotide sugar","Pentose, glucuronate",
            "Fructose, mannose","Galactose","Glutamine","Glutathione","Nitrogen",
            "Glycerophospholipid","Fatty acid","Starch, sucrose","Tyrosine",
            "Arachidonic acid","Linoeic acid",sigs.entero)


gsva.es <- gsva(as.matrix(exprs), pw, verbose=FALSE,method="zscore") #"gsva", "ssgsea", "zscore", "plage"

# score and heatmap-------
library(scales)
scor=apply(gsva.es[-c(14:19),], 2, sum)
ol=which(names(scor)%in%c("GSM929569", "GSM972311","GSM972046")) #outliers

library(ggplot2)
library(cowplot) 
cor.test(scor[-ol],gsva.es[14,-ol])

df=data.frame(x=scor[-ol], y=gsva.es[14,-ol],CMS=cms[-ol])
scatterPlot <-ggplot(df,aes(x=x, y=y)) + 
  geom_point(aes(color=CMS)) + 
  geom_smooth(method=lm, color="black")+
  geom_rug(aes(color=CMS)) +
  scale_color_manual(values = alpha(cms.cols,0.5)) + 
  labs(title=NULL,x="Methabolic Score", y = "Enterocyte Score")+
  theme_test() +
  theme(legend.position=c(0.15,0.85),plot.margin = unit(c(0,0,0,0), "mm"))

xdensity <- ggplot(df, aes(x, fill=CMS)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = cms.cols) + 
  scale_x_continuous(limits = c(-50,50),
                      expand = c(0,0))+
  theme_void() + 
  theme(legend.position = "none",plot.margin = unit(c(0,0,0,12), "mm"))

ydensity <- ggplot(df, aes(y, fill=CMS)) + 
  geom_density(alpha=.5) + 
  coord_flip() +
  scale_fill_manual(values = cms.cols) +
  scale_x_continuous(limits = c(-17,13),
                     expand = c(0,0))+
  theme_void() + 
  theme(legend.position = "none",plot.margin = unit(c(0,0,10,0), "mm"))


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 0.5), heights=c(0.5, 4))

