##PCA

###Change directory based on where the data is deposited

library(biomaRt) 
library(readxl)
library(factoextra)
library(FactoMineR)
library(caret)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("/Users/u6038961/Documents/R/Deseq workflow")

Prostate = read_excel("TL_RNAv1_Prostate_Protein_Coding_MINUS_TL.20.8213D3.xlsx")

Prostate
ncol(Prostate) #30 tls
#drop the columns of hgnc and biotype, will give us issues if not removed. Purpose no longer needed

Prostate <- Prostate[ -c(2:3) ]

Prostate

###genes are row names
Prostate = Prostate %>% remove_rownames %>% column_to_rownames(var="Geneid")

Prostate
ncol(Prostate)

##Remove whole rows with less than 10 total reads
threshold = 10
culsumzero = rownames(Prostate)[rowSums(Prostate)>threshold]

PCA_Cluster_KD_Variance = Prostate[culsumzero,]
PCA_Cluster_KD_Variance[is.na(PCA_Cluster_KD_Variance)] <- 0

##Create a dataframe with the conditions to evaluate to produce the DeSEQ2 object
#coldata = data.frame(ind = colnames(PCA_Cluster_KD_Variance),condition=[(c("")
coldata = read_excel("Data.xlsx")
library(tidyr)
#coldata =coldata %>% drop_na()
coldata
nrow(coldata)
##coldata$condition[c(3:5)]="B"

###ind as rownames
coldata = coldata %>% remove_rownames %>% column_to_rownames(var="ID")
coldata
coldata$De_Novo = as.factor(coldata$De_Novo)

###Plot PCA before DESeq2
mtcars.pca_B <- prcomp(as.matrix(PCA_Cluster_KD_Variance),center=T,scale. = T)

##Get the groups to color from the metadata file
habi = coldata$De_Novo

###PCA for individuals
Pre_DESeq2_PCA = fviz_pca_var(mtcars.pca_B,pointsize = 5,pointshape = 19,
                         #geom = c("point","text"),invisible="quali",
                         geom = c("point"),invisible="quali",
                         habillage =  habi,axes = c(2,3), # Color by the quality of representation
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = F, title = "PCA_Wildcat"    # Avoid text overlapping
)
Pre_DESeq2_PCA



#coldata$De_Novo = as.factor(coldata$De_Novo)
coldata$De_Novo = as.factor(coldata$De_Novo)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = PCA_Cluster_KD_Variance,
                              colData = coldata,
                              design = ~De_Novo)


#PCA_Cluster_KD_Variance
##Prefiltering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

##There are two options VST or rlog, they are similiar but rlog get an extra
##shrinkage element on each sample and each gene, but if it takes a long time then use vst
vsd <- vst(dds, blind=FALSE)
rld <- vst(dds, blind=FALSE)

##Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))

library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# this gives log2(n + 1)
ntd <- normTransform(dds)

library(vsn)



meanSdPlot(assay(rld))

##Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:1000]
df <- as.data.frame(colData(dds))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

##This is from DeSEq2
plotPCA(rld, intgroup=c("De_novo"))


pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


###This is my code to generate the PCA

##Get a dataframe from the rlog normalization
vsd_df = assay(rld)

##apply the PCA function
#mtcars.pca_B <- prcomp(as.matrix(vsd_df), center = TRUE,scale. = TRUE)
mtcars.pca_B <- prcomp(as.matrix(vsd_df))

##Get the groups to color from the metadata file
  habi = coldata$De_Novo

###PCA for individuals
Final_PCA = fviz_pca_var(mtcars.pca_B,pointsize = 5,pointshape = 19,
                         geom = c("point","text"),invisible="quali",
                         #geom = c("point"),invisible="quali",
                         habillage =  habi,axes = c(1,3), # Color by the quality of representation
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = F, title = "PCA_Wildcat"    # Avoid text overlapping
)
Final_PCA

###PCA for Genes
Final_PCA_genes = fviz_pca_ind(mtcars.pca_B,pointsize = 5,pointshape = 19,
                         #geom = c("point","text"),invisible="quali",
                         geom = c("point"),invisible="quali",axes = c(1,2), # Color by the quality of representation
                         #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = F, title = "PCA_Wildcat"    # Avoid text overlapping
)
Final_PCA_genes

dds <- DESeq(dds)
res <- results(dds)
res




names =coldata$De_Novo


dds$Volume_of_Disease
res <- results(dds)
res
#writing off data
write.csv(as.data.frame((res)), file ="Result.csv")

summary(res)

#data =mcols(,use.names=TRUE)
data

