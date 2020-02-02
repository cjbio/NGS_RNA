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
countdata_clean <- countdata_f[,-c(1,6,31,32,43,44)]
colnames(countdata_clean)

test4 <- countdata_clean[, -grep("1h", colnames(countdata_clean))]
test4 <- test4[, -grep("14d", colnames(test4))]
test4 <- test4[, -grep("2h", colnames(test4))]
test4 <- test4[, -grep("21d", colnames(test4))]
test4 <- test4[, -grep("dex_6", colnames(test4))]
test4 <- test4[, -grep("dex_24", colnames(test4))]
test4 <- test4[, -grep("3d", colnames(test4))]
test4 <- test4[, -grep("7d", colnames(test4))]
test4 <- test4[, -grep("ost_1", colnames(test4))]
test4 <- test4[, -grep("ost_2", colnames(test4))]
test4 <- test4[, -grep("ost_6", colnames(test4))]
test4 <- test4[, -grep("ost_24", colnames(test4))]
test4 <- test4[, -grep("6h", colnames(test4))]
test4 <- test4[, -grep("ost 1", colnames(test4))]
test4 <- test4[, -grep("ost 2", colnames(test4))]
test4 <- test4[, -grep("ost 6", colnames(test4))]
test4 <- test4[, -grep("1h", colnames(test4))]
test4 <- test4[, -grep("1h", colnames(test4))]
test4 <- test4[, -grep("1h", colnames(test4))]
test4 <- test4[, -grep("1h", colnames(test4))]

colnames(test4)
countdata_P5_2 <- test4

#make DESeq samples sheet
#write.csv(colnames(countdata_clean), "samples.csv")

