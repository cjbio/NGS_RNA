
setwd("~/rnaseq/ex0080/R/")
rm(list=ls())

countdata <- read.table("~/rnaseq/ex0080/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata)


#write.csv(colnames(countdata), "samplenames.csv")


renamecol <-read.csv("rename_sample.csv", header=TRUE)
renamecol <- t(renamecol)
samples <- renamecol[2,]
# rename col

colnames(countdata) <- samples
# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)



#same for run2
countdata2 <- read.table("~/rnaseq/ex0080_2/alignments/featurecounts_dedup.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata2 <- countdata2[ ,6:ncol(countdata2)]
colnames(countdata2)

#write.csv(colnames(countdata2), "samplenames2.csv")

renamecol2 <-read.csv("rename_sample2.csv", header=TRUE)
renamecol2 <- t(renamecol2)
samples2 <- renamecol2[2,]
# rename col

colnames(countdata2) <- samples2
# Convert to matrix
countdata2 <- as.matrix(countdata2)
colnames(countdata2)


countdata_ex0080 <- cbind(countdata, countdata2)
colnames(countdata_ex0080)

#make DESeq samples sheet
#write.csv(colnames(countdata_f), "samples.csv")




#
dex_conditions <-read.csv("sample_conditions2.csv", header=TRUE)
treatment <- factor(as.character(dex_conditions[,5]))
condition <- factor(as.character(dex_conditions[,6]))
subject <- factor(as.character(dex_conditions[,7]))
agegroup <-  factor(as.character(dex_conditions[,8]))
# Assign condition (first four are controls, second four contain the expansion)
#condition2 <- factor(c(rep("0", 4), rep("0.5", 4), rep("1", 4), rep("2", 4),rep("6", 4),rep("12", 4),rep("16", 4),rep("24", 4),rep("336", 4),rep("504", 4)))
#subject <- factor(c(rep(c(1,2,3,4),10)))
#agegroup <-  factor(c(rep(c("young","young","old","old"),10)))
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata_ex0080 <- data.frame(row.names=colnames(countdata_f), treatment, subject, condition, agegroup)
###################################################################################################################


setwd("~/rnaseq/ex0072/ex0072-95626531/run1/R")


countdata1 <- read.table("~/rnaseq/ex0072/ex0072-95626531/run1/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1)
countdata2 <- read.table("~/rnaseq/ex0072/ex0072-95626531/run2/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata1 <- countdata1[ ,6:ncol(countdata1)]
colnames(countdata1)
countdata2 <- countdata2[ ,6:ncol(countdata2)]
colnames(countdata2)
countdata <- cbind(countdata1,countdata2)
#write.csv(colnames(countdata), "samplenames1.csv")


renamecol <-read.csv("rename_samplenames1.csv", header=TRUE)
renamecol <- t(renamecol)
samples <- renamecol[2,]
# rename col

colnames(countdata) <- samples
# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

test <- countdata[, -grep("dmem", colnames(countdata))]
test2 <- test[, -grep("adi", colnames(test))]
countdata_f1 <- test2[, -grep("amem", colnames(test2))]
colnames(countdata_f1)


#same for run3 and run4
countdata3 <- read.table("~/rnaseq/ex0072/ex0072-95626531/run3/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1)
countdata4 <- read.table("~/rnaseq/ex0072/ex0072-95626531/run4/quantification/featurecounts_dedup.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata3 <- countdata3[ ,6:ncol(countdata3)]
colnames(countdata3)
countdata4 <- countdata4[ ,6:ncol(countdata4)]
colnames(countdata4)
countdata_3_4 <- cbind(countdata3,countdata4)
#write.csv(colnames(countdata_3_4), "samplenames2.csv")

renamecol2 <-read.csv("rename_samplenames2.csv", header=TRUE)
renamecol2 <- t(renamecol2)
samples2 <- renamecol2[2,]
# rename col

colnames(countdata_3_4) <- samples2
# Convert to matrix
countdata_3_4 <- as.matrix(countdata_3_4)
colnames(countdata_3_4)

test3 <- countdata_3_4[, -grep("amem", colnames(countdata_3_4))]
countdata_f2 <- test3[, -grep("adipo", colnames(test3))]
colnames(countdata_f2)

countdata_f <- cbind(countdata_f1, countdata_f2)
colnames(countdata_f)

#remove repeated 31177 MSC 0, and 121217_p2, all F21
countdata_ex0072 <- countdata_f[,-c(1,6,31,32,43,44)]
colnames(countdata_ex0072)

#make DESeq samples sheet
#write.csv(colnames(countdata_clean), "samples.csv")




#
dex_conditions <-read.csv("samples_conditions.csv", header=TRUE)
condition <- factor(as.character(dex_conditions[,6]))
subject <- factor(as.character(dex_conditions[,7]))
agegroup <-  factor(as.character(dex_conditions[,10]))
# Assign condition (first four are controls, second four contain the expansion)
#condition2 <- factor(c(rep("0", 4), rep("0.5", 4), rep("1", 4), rep("2", 4),rep("6", 4),rep("12", 4),rep("16", 4),rep("24", 4),rep("336", 4),rep("504", 4)))
#subject <- factor(c(rep(c(1,2,3,4),10)))
#agegroup <-  factor(c(rep(c("young","young","old","old"),10)))
library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata_ex0072 <- data.frame(row.names=colnames(countdata_ex0072), subject, condition, agegroup)
###################################################################################################################



#Add the youngVSold_2 samples (deseq2_younvsold_LRT_test.r)
setwd("~/rnaseq/dexTcYoungVsOld_2/R")

countdata <- read.table("featurecounts_dedup.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
colnames(countdata)
#write.csv(colnames(countdata), "samplenames.csv")
renamecol <-read.csv("rename_samples_v2.csv", header=TRUE)
colnames(renamecol)

samples <- colnames(renamecol)
samples <- samples[-1]

# rename col

colnames(countdata) <- samples
# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)


#names <- colnames(countdata)
#gsub("^.*?_",".",names)
#names[c(12,6,3,5,8,9,11,1,2,4,7,10)]

#seperate the samples to DMEM and dex
dmemsamples <-read.csv("dmem_samples.csv", header=TRUE)
dmemlist <- as.numeric(dmemsamples[1,])
dmem_countdata <- countdata[,dmemlist]
dex_countdata <- countdata[,-dmemlist]
colnames(dex_countdata)

#remove blk_p7_and timepoint with n=1
blksamples <-read.csv("blk_p7_and_only1timepointsamples.csv", header=TRUE)
blklist <- as.numeric(blksamples[1,])
dex_countdata <- dex_countdata[,-blklist]
colnames(dex_countdata)
#
dex_conditions <-read.csv("dex_only_sample_conditions_2.csv", header=TRUE)
condition <- factor(as.character(dex_conditions[,6]))
subject <- factor(as.character(dex_conditions[,7]))
agegroup <-  factor(as.character(dex_conditions[,10]))
# Assign condition (first four are controls, second four contain the expansion)
#condition2 <- factor(c(rep("0", 4), rep("0.5", 4), rep("1", 4), rep("2", 4),rep("6", 4),rep("12", 4),rep("16", 4),rep("24", 4),rep("336", 4),rep("504", 4)))
#subject <- factor(c(rep(c(1,2,3,4),10)))
#agegroup <-  factor(c(rep(c("young","young","old","old"),10)))
#library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(dex_countdata), subject, condition, agegroup)

#remove the MSC timepoint (tp) that is not at 0h
colnames(dex_countdata[,c(67,69,70,71)])
countdata_YvO_2 <- dex_countdata[,-c(67,69,70,71)]
colnames(countdata_YvO_2)
coldata_YvO_2 <- coldata[-c(67,69,70,71),]




# combine the two sets together!
countdata_YvO_3 <- cbind(countdata_ex0080, countdata_ex0072, countdata_YvO_2)

##reaarange and sub "msc" to "dex" at 0 timepoint
coldata_ex0080 <- coldata_ex0080[,c("subject", "condition", "agegroup", "treatment")]
coldata_ex0080$treatment <- sub("msc", "dex", coldata_ex0080$treatment)

#add treatment catergory for ex0072 and YvO_2
coldata_ex0072$treatment <- "dex"
coldata_YvO_2$treatment <- "dex"

coldata_YvO_3 <- rbind(coldata_ex0080 ,coldata_ex0072, coldata_YvO_2)
condition <- coldata_YvO_3$condition
subject <- coldata_YvO_3$subject
agegroup <- coldata_YvO_3$agegroup

setwd("~/rnaseq/dexTcYoungVsOld_3_final/")


###remove 31177
coldata_YvO_3_clean <- coldata_YvO_3[!coldata_YvO_3$subject=="Young31177",]
countdata_YvO_3_clean <- countdata_YvO_3[,!coldata_YvO_3$subject=="Young31177"]

coldata_YvO_3_clean_2 <- coldata_YvO_3_clean[!(coldata_YvO_3_clean$subject=="Young1309" & coldata_YvO_3_clean$condition=="336"),]
countdata_YvO_3_clean_2 <- countdata_YvO_3_clean[,!(coldata_YvO_3_clean$subject=="Young1309" & coldata_YvO_3_clean$condition=="336")]

#remove old6 0 outliner
coldata_YvO_3_clean_3 <- coldata_YvO_3_clean_2[-57,]
countdata_YvO_3_clean_3  <- countdata_YvO_3_clean_2 [,-57]
