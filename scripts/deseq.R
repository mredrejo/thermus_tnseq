library(DESeq2)
library(ggplot2)
#meethod from https://lashlock.github.io/compbio/R_presentation.html

countData <- read.csv('tnseq_counts_scores.csv', header = TRUE, sep = ";")
countData <- countData[,c(1,2,5,8,11)]
countData[,2:5] <- lapply(countData[,2:5], as.integer)
head(countData)
metaData <- data.frame(names(countData[,2:5]),c("HB27","HB27","Ppol","Ppol"),c("Exp1","Exp2","Exp1","Exp2"))

names(metaData) <- c("Sample","Strain","Experiment")

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Strain, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res <- res[order(res$padj),]
head(res)

#PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Strain")



tnseq.pca <- prcomp(tnseq[,c(4,7,10,13)], 
                   center = TRUE, 
                   scale. = TRUE) 
biplot(tnseq.pca) 
library(ggfortify) 
iris.pca.plot <- autoplot(tnseq.pca, 
                          data = stack(tnseq[,c(4,7,10,13)]), 
                          colour = 'ind') 


#volcano
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-7,7)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="cyan"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


subset(res, padj<.01 & abs(log2FoldChange)>2)


