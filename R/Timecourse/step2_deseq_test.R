
setwd("~/rnaseq/dexTcYoungVsOld_3_final/")

coldata_YvO_3
nrow(coldata_YvO_3_clean_3)
ncol(countdata_YvO_3_clean_3)

#remove bio treatment
coldata_YvO_3_dex <- coldata_YvO_3_clean_3[!coldata_YvO_3_clean_3$treatment=="bio",]
countdata_YvO_3_dex <- countdata_YvO_3_clean_3[,!coldata_YvO_3_clean_3$treatment=="bio"]

nrow(coldata_YvO_3_dex)
ncol(countdata_YvO_3_dex)
agegroup <- coldata_YvO_3_dex$agegroup
condition <- coldata_YvO_3_dex$condition
treatment <- coldata_YvO_3_dex$treatment
# make dataframe
dds <- DESeqDataSetFromMatrix(countData=countdata_YvO_3_dex, colData=coldata_YvO_3_dex, design=~agegroup + condition + agegroup:condition)

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


dds
View(mcols(dds_clean))

resTC <- as.data.frame(results(dds_clean))
summary(resTC)
resTC$symbol <- mcols(dds_clean)$symbol
head(resTC[order(resTC$padj),], 4)
orderedTC <- as.data.frame(resTC[order(resTC$padj),])

sigfilter <- orderedTC$padj<0.05
sigTC <- orderedTC[sigfilter,]
sigTC <- na.omit(sigTC)



#filter nc by significant genes
nc_sig <- as.data.frame(filtered_nc[rownames(filtered_nc) %in% rownames(sigTC),])
ordered_nc_sig <- nc_sig[order(rownames(nc_sig)),]
write.csv(ordered_nc_sig, "sig_genes_clean.csv")

siggeneID <- read.csv("sig_genes_clean_names.csv")
rownames(siggeneID)<- siggeneID[,1]
ordered_siggeneID <- siggeneID[order(rownames(siggeneID)),]

named_nc_sig <- cbind(ordered_siggeneID, ordered_nc_sig[1:63,]) 
rownames(named_nc_sig) <- named_nc_sig[,2]
named_nc_sig<- named_nc_sig[,-c(1,2)]

######################################################################################################################################
###remove all 31177 7d, 21d, 1209 14d, and try again
colnames(countdata_YvO3_c)
countdata_YvO3_c <- countdata_YvO_3[,-c(3:11)]
colnames(countdata_YvO3_c)

coldata_YvO_3_c <- coldata_YvO_3[-c(3:11),]
condition <- coldata_YvO_3_c$condition
subject <- coldata_YvO_3_c$subject
agegroup <- coldata_YvO_3_c$agegroup
coldata_YvO_3_c[rownames(coldata_YvO_3_c)=="old 5 ost 7d",]
coldata_YvO_3_c[rownames(coldata_YvO_3_c)=="old 5 ost 24",]
View(coldata_YvO_3_c)

dds_c <- DESeqDataSetFromMatrix(countData=countdata_YvO3_c, colData=coldata_YvO_3_c, design=~agegroup + condition + agegroup:condition)

#filter normalized count of at least 10 in two or more samples
dds_c <- estimateSizeFactors(dds_c)
nc_c <- counts(dds_c, normalized=TRUE)
filter <- rowSums(nc_c >= 20) >= 2
dds_c <- dds_c[filter,]
filtered_nc_c <- counts(dds_c, normalized=TRUE)
# Run the DESeq pipeline
dds_c <- DESeq(dds_c, test = "LRT", reduced = ~ agegroup + condition)
#dds_c <- DESeq(dds_c)
resTC <- as.data.frame(results(dds_c))
summary(resTC)
resTC$symbol <- mcols(dds_c)$symbol
head(resTC[order(resTC$padj),], 4)
orderedTC <- as.data.frame(resTC[order(resTC$padj),])

sigfilter <- orderedTC$padj<0.05
sigTC <- orderedTC[sigfilter,]
sigTC <- na.omit(sigTC)



#filter nc by significant genes
nc_sig <- as.data.frame(filtered_nc[rownames(filtered_nc) %in% rownames(sigTC),])
ordered_nc_sig <- nc_sig[order(rownames(nc_sig)),]
write.csv(ordered_nc_sig, "sig_genes_clean.csv")

siggeneID <- read.csv("sig_genes_clean_names.csv")
rownames(siggeneID)<- siggeneID[,1]
ordered_siggeneID <- siggeneID[order(rownames(siggeneID)),]

named_nc_sig <- cbind(ordered_siggeneID, ordered_nc_sig[1:63,]) 
rownames(named_nc_sig) <- named_nc_sig[,2]
named_nc_sig<- named_nc_sig[,-c(1,2)]



################################with BIO


#with bio treatment

agegroup_bio <- coldata_YvO_3_clean_3$agegroup
condition_bio <- coldata_YvO_3_clean_3$condition
treatment_bio <- coldata_YvO_3_clean_3$treatment
# make dataframe
dds_bio <- DESeqDataSetFromMatrix(countData=countdata_YvO_3_clean_3, colData=coldata_YvO_3_clean_3, design=~agegroup + condition + agegroup:condition)

#filter normalized count of at least 10 in two or more samples
dds_bio <- estimateSizeFactors(dds_bio)
nc_bio <- counts(dds_bio, normalized=TRUE)
filter_bio <- rowSums(nc >= 10) >= 2
dds_bio <- dds_bio[filter_bio,]
filtered_nc_bio <- counts(dds_bio, normalized=TRUE)


save.image("dds_31177rm_old5MSC0rm.RData")




