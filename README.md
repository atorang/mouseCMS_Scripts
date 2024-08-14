# mouseCMS_Scripts
Scripts for regenerating figures in the paper entitled:"Enterocyte-like differentiation defines the metabolic phenotype of CMS3 colorectal cancers and provides a therapeutic vulnerability."
Download data from https://www.synapse.org/#!Synapse:syn53040345/ to run the scripts.

## Data description:
* organoids.RData: in vitro mouse models generated in this study and available at GSE223323
* invivo.organoids.RData (in vivo mouse models generated in this study and available at GSE270929)
* CellLine.GSE248558.RData (CRC cell lines treated with CPS1 inhibition and untreated controls available at GSE248558)
* Synapse.RData (human gene expression datasets proflied using U133 Plus 2.0 platform collected in the work of Guinney et al.)
* Synapse.High.RData (human gene expression datasets proflied using U133 Plus 2.0 platform collected in the work of Guinney et al. + TCGA CRC samples)
* scRNAseq.colon.GSE125970.RData (single cell data of normal colon collected from GSE125970)
* msigdb.v7.4.symbols.gmt (gene signatures obtained from Broad Institute)
* Synapse.genesets.gmt (gene signatures obtained from work of Guinney et al.)
* jem_20191130_tables2.xlsx (gene signatures of human cell types developed by Wang et al.)
* Haber.Combine.All.xlsx (gene signatures of 	mouse SI cell types developed by Haber et al.)
* GSE146476_Intestinal-region-specific-Wnt-signalling_gene_counts.txt (count data directly obtained from publicly available data at GSE146476)


## Scripts to process raw fastq files available in GEO to count data:
* fastq2count.organoids.GSE223323.txt (GSE223323 and GSE270929)
* fastq2count.cellline.GSE248558.txt (GSE248558)

## Scripts to normalize/remove batch effect in count data and their final file names. 
In the final files both count data and normalized data are available:
* data.handling.organoids.R --> organoids.RData (GSE223323)
* data.handling.invivo.organoids.R --> invivo.organoids.RData (GSE270929)
* data.handling.cellLine.GSE248558.R --> CellLine.GSE248558.RData (GSE248558)

## Scripts for conducted analyses and figures:
### Human gene expression datasets:
* Scripts: data.handling.Synapse.R
### GSEA in patient samples: 
* Scripts: Fig.4C.Synapse.Enterocyte.R, Fig.S4C.GSEA.KRAS.MUTvsWT.R
* Datasets: Synapse.RData, msigdb.v7.4.symbols.gmt 
### GSEA in CRC cell lines: 
* Scripts: Fig.6E.S5A.GSEA.CPS1.Inhibition.R
* Datasets: CellLine.GSE248558.RData, msigdb.v7.4.symbols.gmt, jem_20191130_tables2.xlsx 
### GSEA in scRNAseq data: 
* Scripts: Fig.4B.CellType.Metabolism.R
* Datasets: scRNAseq.colon.GSE125970.RData, msigdb.v7.4.symbols.gmt
### GSEA in mouse organoids: 
* Scripts: Fig.4A.Genotype.Cell.All.R, Fig.3C.gsplot.cms.R
* Datasets: organoids.RData, Haber.Combine.All.xlsx, Synapse.genesets.gmt, 				msigdb.v7.4.symbols.gmt
### DE analysis in wild type colon and small intestine organoids:
* Scripts: Fig.S2.R
* Datasets: GSE146476_Intestinal-region-specific-Wnt-signalling_gene_counts.txt  			organoids.RData
### CMS classifier for mouse models: 
* Scripts: First.CMS.Stratification.R, Training.Mouse.Classifier.R, Fig.S3DE.Classifier.Coef.R, 		https://github.com/atorang/mouseCMS
* Datasets:  organoids.RData, Synapse.genesets.gmt, msigdb.v7.4.symbols.gmt 
### CMS stratification of mouse organoids: 
* Scripts: Fig.3B.classification.R, https://github.com/atorang/mouseCMS
* Datasets: organoids.RData
### Metabolic and enterocyte scoring:
* Scripts: Fig.4D.Met.Ent.Scores.R
* Datasets: Synapse.High.RData, msigdb.v7.4.symbols.gmt
### ssGSEA: 
* Scripts: Fig.S4B.ssGSEA.organoids.R, Fig.6B.ssGSEA.CellLines.R
* Datasets: organoids.RData, Synapse.genesets.gmt, msigdb.v7.4.symbols.gmt 
