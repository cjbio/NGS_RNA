# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

plotDispEsts(dds2, main="Dispersion plot")

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
rld <- vst(dds_bio)

head(assay(rld))
hist(assay(rld))

# Colors for plots below
## Ugly:
#(mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
#library(RColorBrewer)
(mycols <- brewer.pal(6, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis

#for colour coding the young and old donors later on (chao's own code)

#subjectnumb <-read.csv("subject_dex_only_blk_1tp_removed.csv", header=TRUE)
agegroup_bio <- coldata_YvO_3_clean_3$agegroup
condition_bio <- coldata_YvO_3_clean_3$condition
treatment_bio <- coldata_YvO_3_clean_3$treatment

ageTreatList <- as.numeric(factor(paste(agegroup_bio, treatment_bio), levels=c("old dex", "young bio", "old bio", "young dex")))

subjectlist <- as.numeric(as.factor(coldata_YvO_3_clean_3$subject))
agegrouplist <- as.numeric(as.factor(agegroup_bio))
conditionlist <- as.numeric((as.factor(condition_bio)))
treatmentlist <- as.numeric((as.factor(treatment_bio)))

lengendCol <- c(0.25, 0.75, 0.5, 1)

## Could do with built-in DESeq2 function:
DESeq2::plotPCA(vst, intgroup="agegroup")
DESeq2::plotPCA(rld, intgroup="condition")
DESeq2::plotPCA(rld, intgroup="agegroup")
#?plotPCA
#install.packages("calibrate")
library("calibrate")
## I like mine better:
rld_pca <- function (rld, intgroup = c("condition","subject", "treatment"), ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA plot", textcx=0.55, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  write.csv(pca$rotation[,"PC1"], "PC1-vst_dex_bio.csv")
  write.csv(pca$rotation[,"PC2"], "PC2-vst_dex_bio.csv")
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac],pch=19, col=hsv(h=ageTreatList/4,s=ageTreatList/4,v=ageTreatList/4, alpha = log(as.numeric(levels(condition_bio)[condition_bio])+1.5,600)), xlab=pc1lab, ylab=pc2lab, main=main, ...)
  #with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=c("old dex","old bio", "young bio", "young dex"), col=hsv(h=lengendCol, s=lengendCol, v=lengendCol), pch=19)
  #    rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #          pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca-adipo.png", 2400, 1800, pointsize=40)
rld_pca(rld, intgroup=c("condition","subject", "treatment"))#, xlim=c(-30, 100))
dev.off()

?hsv()
#to convert the factor back to its original numer use : as.numeric(levels(x)[x])
log(as.numeric(levels(condition)[condition])+1,600)
## MA plot
## Could do with built-in DESeq2 function:
DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
