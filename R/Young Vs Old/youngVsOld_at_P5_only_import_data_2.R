
setwd("~/rnaseq/ex0080/R/")


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


countdata_f <- cbind(countdata, countdata2)
colnames(countdata_f)

test5 <- countdata_f[, grep("msc_0", colnames(countdata_f))]
colnames(test5)
countdata_P5_3 <- test5
