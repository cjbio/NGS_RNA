#markers <- read.csv("~/database/markerpanel.csv", header = TRUE)
markers <- read.csv("~/database/wnt target genes human ensemble ID final updated.csv", header = TRUE)
rownames(markers) <- markers[,1]
ordered_markers <- markers
#filter nc by marker panel
#sig_filtered_nc <- as.data.frame(filtered_nc[rownames(filtered_nc) %in% rownames(sigTC),])
nc_markers <- as.data.frame(filtered_nc[rownames(filtered_nc) %in% rownames(markers),])
ordered_nc_markers <- nc_markers[rownames(markers),]

ordered_nc_markers <- na.omit(ordered_nc_markers)

rownames(ordered_markers)
rownames(ordered_nc_markers)

ordered_markers_found <- as.data.frame(ordered_markers[rownames(ordered_markers) %in% rownames(ordered_nc_markers),])


named_nc_markers <- cbind(ordered_markers_found, ordered_nc_markers) 
rownames(named_nc_markers) <- named_nc_markers[,2]
named_nc_markers <- named_nc_markers[,-c(1:2)]

##################boxplot
coldata$condition <- relevel(coldata$condition, ref = "young")
pdf("boxplot_sigDE_wnt_target_genes_thesis.pdf")
for(i in rownames(named_nc_markers)){
  boxcounts <- as.data.frame(t(named_nc_markers))
  boxcounts <- as.data.frame(cbind(boxcounts[,i],coldata$condition))
  boxcounts$V2 <- coldata$condition
  boxplot <- ggplot(boxcounts, aes(x=V2, y=V1)) +
    geom_boxplot(color = c("red", "black")) +
    geom_point(colour = c(rep("red",12), rep("black",9)), cex=2) +
    ggtitle(i) +
    xlab("age group") +
    ylab("normalised counts") +
    theme_classic() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=18))
  plot(boxplot)
}
dev.off()

?geom_point
#####barplot of each gene expression

pdf("barchart sigDE wnt target genes.pdf")
par(mfrow=c(2,2))
for(i in rownames(named_nc_markers)){
  barcounts <- as.data.frame(t(named_nc_markers))
  barcounts <- barcounts[,i]
  barplot(barcounts, names.arg = "age group", main = i, 
          col = c(rep("light green",12), rep("light grey", 9)), ylab = "FPKM")
  legend(1,1, legend = c("young", "old"), pch = 15, col = c("light green", "light grey"))
}
dev.off()


# name base on the age group and timepoint/condition
tp_age <- paste(as.character(coldata$condition), as.character(coldata$agegroup))
colnames(named_nc_markers) <- tp_age

test1 <- named_nc_markers



# calculate the mean of the columns with the same column name
test2 <- as.data.frame(
  sapply(unique(names(test1)), 
         function(col) rowMeans(test1[names(test1) == col]) # calculate row means
  )
)


colnames(test2)



library(qPCRheatmaps)
library(pheatmap)


wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}



ordered <- test2

#fcdexp9<-ordereddexp9_fpkm/ordereddexp9_fpkm[,1]
ordered_fpkm<-ordered+1


fc<-ordered_fpkm/ordered_fpkm[,1]


#pdf("ed_script_heatmaps_dex_p9.pdf")
#starting the heatmaps analysis
#correlation of pairwise comparison of overall gene expression between one time point and another
corr_fpkm<-cor(ordered_fpkm)
corr_fc<-cor(fc)

pheatmap(corr_fpkm^2,cluster_rows=F, cluster_cols=F, main="Correlation of FPKM values between timepoints", cellwidth=20,cellheight=20)
pheatmap(corr_fc^2, cluster_cols=F,cluster_rows=F, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)

set.seed(1234)
pheatmap(ordered_fpkm, cluster_rows=T,cluster_cols=F, show_rownames=T, main="heatmap of FPKM")

pheatmap(log2(ordered_fpkm), cluster_cols=F, show_rownames=T, main="heatmap of log2 FPKM")

pheatmap(log2(fc), cluster_cols=F, show_rownames=T, main="heatmap of log2 FPKM")




qPCRheatmap(log2(ordered_fpkm), limit=50, cluster_rows=F, cluster_cols=F, show_rownames=T, main="Heatmap of log10 FPKM")
qPCRheatmap(fc, limit=50, cluster_rows=T, cluster_cols=F, show_rownames=T,fontsize_row =6, main="Heatmap of fold changes")
pheatmap(fc, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)
dev.off()

###########################again for inhibitors - too lazy to write loop########################
# calculate the mean of the columns with the same column name
test2 <- as.data.frame(
  sapply(unique(names(test1_inh)), 
         function(col) rowMeans(test1_inh[names(test1_inh) == col]) # calculate row means
  )
)


colnames(test2)



library(qPCRheatmaps)
library(pheatmap)


wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}



ordered <- test2

#fcdexp9<-ordereddexp9_fpkm/ordereddexp9_fpkm[,1]
ordered_fpkm<-ordered+1


fc<-ordered_fpkm/ordered_fpkm[,1]


#pdf("ed_script_heatmaps_dex_p9.pdf")
#starting the heatmaps analysis
#correlation of pairwise comparison of overall gene expression between one time point and another
corr_fpkm<-cor(ordered_fpkm)
corr_fc<-cor(fc)

pheatmap(corr_fpkm^2,cluster_rows=F, cluster_cols=F, main="Correlation of FPKM values between timepoints", cellwidth=20,cellheight=20)
pheatmap(corr_fc^2, cluster_cols=F,cluster_rows=F, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)

set.seed(1234)
pheatmap(ordered_fpkm, cluster_rows=T,cluster_cols=F, show_rownames=T, main="heatmap of FPKM")

pheatmap(log2(ordered_fpkm), cluster_cols=F, show_rownames=T, main="heatmap of log2 FPKM")

pheatmap(log2(fc), cluster_cols=F, show_rownames=T, main="heatmap of log2 FPKM")




qPCRheatmap(log2(ordered_fpkm), limit=50, cluster_rows=F, cluster_cols=F, show_rownames=T, main="Heatmap of log10 FPKM")
qPCRheatmap(fc, limit=50, cluster_rows=T, cluster_cols=T, show_rownames=T,fontsize_row =6, main="Heatmap of fold changes")
pheatmap(fc, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)
dev.off()


save.image("dex+BIO_treatment.RData")