library(DESeq2)

filtered_nc[rownames(filtered_nc)=="ENSG00000081059",coldata_YvO_3_dex$condition==504]

countdata_YvO_3_dex <- countdata_YvO_3_dex[,!colnames(countdata_YvO_3_dex)=="old7_dex_21d"]
coldata_YvO_3_dex <- coldata_YvO_3_dex[!rownames(coldata_YvO_3_dex)=="old7_dex_21d",]

dds <- DESeqDataSetFromMatrix(countData=countdata_YvO_3_dex, 
                              colData=coldata_YvO_3_dex, 
                              design=~agegroup + condition + agegroup:condition)

#filter normalized count of at least 10 in two or more samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]
filtered_nc <- counts(dds, normalized=TRUE)
# Run the DESeq pipeline
#dds <- estimateDispersions(dds)
#dds <- nbinomLRT(dds, reduced = ~ agegroup + condition, maxit = 1000)
dds <- DESeq(dds, test = "LRT", reduced = ~ agegroup + condition)
dds_clean <- dds[which(mcols(dds)$fullBetaConv)]

resTC <- as.data.frame(results(dds_clean))
summary(resTC)
resTC$symbol <- mcols(dds_clean)$symbol
head(resTC[order(resTC$padj),], 4)
orderedTC <- as.data.frame(resTC[order(resTC$padj),])

sigfilter <- orderedTC$padj<0.05
sigTC <- orderedTC[sigfilter,]
sigTC <- na.omit(sigTC)

resTC[rownames(resTC)=="ENSG00000081059",]

###################################################################
filtered_nc[rownames(filtered_nc)=="ENSG00000081059",coldata_YvO_3_dex$condition==6]

countdata_YvO_3_dex <- countdata_YvO_3_dex[,!colnames(countdata_YvO_3_dex)=="old6_ost_6"]
coldata_YvO_3_dex <- coldata_YvO_3_dex[!rownames(coldata_YvO_3_dex)=="old6_ost_6",]

dds <- DESeqDataSetFromMatrix(countData=countdata_YvO_3_dex, 
                              colData=coldata_YvO_3_dex, 
                              design=~agegroup + condition + agegroup:condition)

#filter normalized count of at least 10 in two or more samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]
filtered_nc <- counts(dds, normalized=TRUE)
# Run the DESeq pipeline
#dds <- estimateDispersions(dds)
#dds <- nbinomLRT(dds, reduced = ~ agegroup + condition, maxit = 1000)
dds <- DESeq(dds, test = "LRT", reduced = ~ agegroup + condition)
dds_clean <- dds[which(mcols(dds)$fullBetaConv)]

resTC <- as.data.frame(results(dds_clean))
summary(resTC)
resTC$symbol <- mcols(dds_clean)$symbol
head(resTC[order(resTC$padj),], 4)
orderedTC <- as.data.frame(resTC[order(resTC$padj),])

sigfilter <- orderedTC$padj<0.05
sigTC <- orderedTC[sigfilter,]
sigTC <- na.omit(sigTC)

resTC[rownames(resTC)=="ENSG00000081059",]
write.csv(sigTC, "reAnalysed2Siggenes.csv")
#####################################################################
rm(list=ls())

table(coldata_YvO_3_dex$condition)
