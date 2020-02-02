
setwd("~/rnaseq/passageTcYoungVsOld/R")
rm(list=ls())


# get feature files from several place


#first one
countdata1 <- read.table("~/rnaseq/dexTcYoungVsOld_2/R/featurecounts_dedup.txt", header=TRUE, row.names=1) 

# Remove first five columns (chr, start, end, strand, length)
countdata1 <- countdata1[ ,6:ncol(countdata1)]
colnames(countdata1[c(1,35,51,65,86,118,121,129)])

renamecol1 <-read.csv("~/rnaseq/dexTcYoungVsOld_2/R/rename_samples_v2.csv", header=TRUE)
colnames(renamecol1)

samples <- colnames(renamecol1)
samples <- samples[-1]
# rename col

colnames(countdata1) <- samples
# Convert to matrix

countdata1 <- countdata1[,c(1,35,51,65,86,118,121,129)]
colnames(countdata1)



#second one
countdata2 <- read.table("~/rnaseq/passageTcYoungVsOld/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1) 
# Remove first five columns (chr, start, end, strand, length)
countdata2 <- countdata2[ ,6:ncol(countdata2)]
colnames(countdata2)


#third one
countdata3 <- read.table("~/rnaseq/martin_olddonor1_tc/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1) 
# Remove first five columns (chr, start, end, strand, length)
countdata3 <- countdata3[ ,6:ncol(countdata3)]
countdata3 <- countdata3[,c(1,2)]
names3 <- c("old1_0", "old1_0.5")
colnames(countdata3)<- names3

#combine samples
countdata <- cbind(countdata1,countdata2,countdata3)
countdata <- countdata[,1:30]
colnames(countdata)

#young(1,2,3,4,12,18)
#old (5,6,7,8,30)
#0014 passage(10,14,16)
#0015 passage(17,19,20)
#0025 passage(28,29,21,22,23,24,25,26,27)

countdata <- countdata[,c(1,2,3,12,18,5,6,7,8)]
colnames(countdata)

countdata_P5_1 <- countdata
###############################################################################################

# script
# 1. youngVsOld_at_P5_only_import_data_1.R
# 1. youngVsOld_at_P5_only_import_data_2.R



##########################################################################################
setwd("~/rnaseq/youngVsOld_atP5_final")

countdata <- cbind(countdata_P5_1,countdata_P5_2,countdata_P5_3)
colnames(countdata)



countdata_young <- countdata[,c(1,2,3,11,14,17,18,19,21,22,24,27)]
colnames(countdata_young)
countdata_old <- countdata[,c(6,7,8,9,10,20,23,25,26)]
colnames(countdata_old)

countdata_P5_final <- cbind(countdata_young, countdata_old)
colnames(countdata_P5_final)

#names4 <- c("young2_F","young2.1_F","young4_M","young2.2_F","young3_F","old4_F","old3_M","old3.1_M","old2_M")
#colnames(countdata)<-names4
condition <- factor(c(rep("young", 12),rep("old", 9)))
#subject <- factor(as.character(dex_conditions[,7]))
#agegroup <-  factor(as.character(dex_conditions[,10]))
# Assign condition (first four are controls, second four contain the expansion)
#condition2 <- factor(c(rep("0", 4), rep("0.5", 4), rep("1", 4), rep("2", 4),rep("6", 4),rep("12", 4),rep("16", 4),rep("24", 4),rep("336", 4),rep("504", 4)))
#subject <- factor(c(rep(c(1,2,3,4),10)))
#agegroup <-  factor(c(rep(c("young","young","old","old"),10)))
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata_P5_final), condition)


# make dataframe
dds <- DESeqDataSetFromMatrix(countData=countdata_P5_final, colData=coldata, design=~condition)

#filter normalized count of at least 10 in two or more samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 2
dds <- dds[filter,]
filtered_nc <- counts(dds, normalized=TRUE)
# Run the DESeq pipeline
dds <- DESeq(dds)
dds
results(dds)
#check results

resTC <- as.data.frame(results(dds))
head(resTC)
summary(resTC)
resTC$symbol <- mcols(dds)$symbol
head(resTC[order(resTC$padj),], 4)
orderedTC <- as.data.frame(resTC[order(resTC$padj),])

sigfilter <- orderedTC$padj<0.05
sigTC <- orderedTC[sigfilter,]
sigTC <- na.omit(sigTC)

#######################################################
library(biomaRt)
listMarts(host = "uswest.ensembl.org")
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")

IDs <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = rownames(sigTC),  mart = ensembl)

ordered_IDs <- IDs[order(IDs$ensembl_gene_id),]
ordered_sigTC <- sigTC[order(rownames(sigTC)),]
ordered_sigTC$Names <- rownames(ordered_sigTC)

named_sigTC <- cbind(ordered_IDs, ordered_sigTC[rownames(ordered_sigTC) %in% ordered_IDs$ensembl_gene_id,] )

#write.csv(named_sigTC,"sig_youngVsoldatP5_n=21")

log2Up <- named_sigTC[named_sigTC$log2FoldChange> 1,]
log2Down <- named_sigTC[named_sigTC$log2FoldChange< -1,]
