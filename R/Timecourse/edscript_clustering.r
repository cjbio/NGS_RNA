library(ggplot2)
library(pheatmap)
library(qPCRheatmaps)

#filter for sig upanddown reg genes YOUNG ONLY

#markers_down <- read.csv("../young_only_downregulated_genes_allTC.csv", header = TRUE)
#markers_up <- read.csv("../young_only_upregulated_genes_allTC.csv", header = TRUE)
#markers <- rbind(markers_up, markers_down)
#markers <- markers[order(markers$Gene.stable.ID),]

#write.csv(sigTC, "sigDEG LRT.csv")

markers <- read.csv("sigDEG namesLRT.csv", header = TRUE)
rownames(markers) <- markers[,1]
ordered_markers <- markers


#filter nc by marker panel

nc_markers <- as.data.frame(filtered_nc[rownames(filtered_nc) %in% rownames(markers),])
ordered_nc_markers <- nc_markers[rownames(markers),]

ordered_nc_markers <- na.omit(ordered_nc_markers)

rownames(ordered_markers)
rownames(ordered_nc_markers)

ordered_markers_found <- as.data.frame(ordered_markers[rownames(ordered_markers) %in% rownames(ordered_nc_markers),])


named_nc_markers <- cbind(ordered_markers_found, ordered_nc_markers) 
rownames(named_nc_markers) <- named_nc_markers[,2]
named_nc_markers <- named_nc_markers[,-c(1:2)]

########relevel

condition <- coldata_YvO_3_dex$condition
order(condition)
timepoints <- unique(as.numeric(levels(condition))[condition])
timepoint_order <- timepoints[order(timepoints)]

condition_c <- factor(coldata_YvO_3_dex$condition, levels = timepoint_order)

# name base on the age group and timepoint/condition
tp_age <- paste(as.character(coldata_YvO_3_dex$condition), as.character(coldata_YvO_3_dex$agegroup))
colnames(named_nc_markers) <- tp_age

test1 <- named_nc_markers



# calculate the mean of the columns with the same column name
test2 <- as.data.frame(
  sapply(unique(names(test1)), 
         function(col) rowMeans(test1[names(test1) == col]) # calculate row means
  )
)


colnames(test2)
old_markers <- test2[,c("0 old","0.5 old","1 old","2 old","4 old","6 old","8 old","12 old","16 old","24 old","48 old","72 old","96 old","120 old","168 old","336 old","504 old")]
young_markers <- test2[, c("0 young","0.5 young","1 young","2 young","4 young","6 young","8 young","12 young","16 young","24 young","48 young","72 young","96 young","120 young","168 young","336 young","504 young")]
#"0 old","0.5 old","1 old","2 old","4 old","6 old","8 old","12 old","16 old","24 old","48 old","72 old","96 old","120 old","168 old","336 old","504 old"
#"0 young","0.5 young","1 young","2 young","4 young","6 young","8 young","12 young","16 young","24 young","48 young","72 young","96 young","120 young","168 young","336 young","504 young"
#write.csv(ordered_markers, "youngvsold_markers_63_clean_deseq2.csv")


wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}



ordered <- cbind(young_markers, old_markers)
#ordered["SP7",]
#ordered["ALP",]
#fcdexp9<-ordereddexp9_fpkm/ordereddexp9_fpkm[,1]
ordered_fpkm<-ordered+1
#ordered_fpkm["SP7",]

colnames((ordered_fpkm))
#write.csv(ordered_fpkm, "sig_ordered_fpkm.csv")
fc1<-ordered_fpkm[,c(1:17)]/ordered_fpkm[,1]

fc2<-ordered_fpkm[,c(18:34)]/ordered_fpkm[,18]
fc <- cbind(fc1, fc2)
#fc["SP7",]
#?pdf()
#png("siggnenes_31177rm_heatmaps_.png", width = 1600, height = 1000)
#starting the heatmaps analysis
#correlation of pairwise comparison of overall gene expression between one time point and another
corr_fpkm<-cor(ordered_fpkm)
corr_fc<-cor(fc)

pheatmap(corr_fpkm^2,cluster_rows=F, cluster_cols=F, main="Correlation of FPKM values between timepoints", cellwidth=20,cellheight=20)
pheatmap(corr_fc[-1,-1]^2, cluster_cols=F,cluster_rows=F, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)

set.seed(1234)
pheatmap(ordered_fpkm, cluster_rows=T,cluster_cols=T, show_rownames=F, main="heatmap of FPKM")
pheatmap(log2(ordered_fpkm), cluster_rows=T,cluster_cols=T, show_rownames=F, main="heatmap of log2(FPKM)")
#pheatmap(log2(ordered_fpkm[,c()]), cluster_cols=F, show_rownames=F, main="heatmap of log2 FPKM")

qPCRheatmap(log10(ordered_fpkm+1), cluster_rows=F, cluster_cols=F, show_rownames=T, main="Heatmap of fold changes")
qPCRheatmap(fc[,z], limit=100, cluster_rows=T, cluster_cols=F, show_rownames=T,fontsize_row =6, main="Heatmap of fold changes")
qPCRheatmap(fc, limit=100, cluster_rows=T, cluster_cols=F, show_rownames=T,fontsize_row =6, main="Heatmap of fold changes")

pheatmap(fc, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)

fc1



pdf("wssplot_clusters.pdf")
wssplot(young_markers, nc=50)
wssplot(old_markers,nc=50)
wssplot(ordered_fpkm,nc=50)
dev.off()


pdf("ed_script_numberofclusters_p9.pdf")

for (i in 2:5)
{
  qPCRheatmap(fc1[!rownames(fc1)=="MTND1P23",], limit=100,cluster_cols=F, cluster_rows=F ,main=i, kmeans_k=i)
}

for (i in 2:5)
{
  qPCRheatmap(fc2[!rownames(fc1)=="MTND1P23",], limit=100,cluster_cols=F,cluster_rows=F, main=i, kmeans_k=i)
}


#pull out the clustering information ##km is a variable from kmeans from within the modified pheatmap engine used to drive qPCRheatmap

clusters<-km$cluster
clusters
##itterate through clusters extracting genes

Members_1<-names(clusters[which(clusters==1)])
Members_2<-names(clusters[which(clusters==2)])
Members_3<-names(clusters[which(clusters==3)])
Members_4<-names(clusters[which(clusters==4)])
Members_5<-names(clusters[which(clusters==5)])


write.table(Members_1, file="young_Members_1", row.names = FALSE, quote= FALSE, col.names=F)
write.table(Members_2, file="young_Members_2", row.names = FALSE, quote= FALSE, col.names=F)
write.table(Members_3, file="young_Members_3", row.names = FALSE, quote= FALSE, col.names=F)
write.table(Members_4, file="young_Members_4", row.names = FALSE, quote= FALSE, col.names=F)
write.table(Members_5, file="young_Members_5", row.names = FALSE, quote= FALSE, col.names=F)



Members_array_1<-fc1[rownames(fc1) %in% Members_1,]
Members_array_2<-fc1[rownames(fc1) %in% Members_2,]
Members_array_3<-fc1[rownames(fc1) %in% Members_3,]
Members_array_4<-fc1[rownames(fc1) %in% Members_4,]
Members_array_5<-fc1[rownames(fc1) %in% Members_5,]






pdf("young_cluster1.pdf")
#cluster1
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(log(timepoint_order,2),km$centers[1,], xaxt="n", main="Cluster 1",type="b",xlab="timepoint (hr)", ylab="Fold change")
axis(1, labels=timepoint_order, at=log(timepoint_order,2))
for (i in 1:length(Members_array_1[,1]))
{
  plot(log(timepoint_order,2),Members_array_1[i,], xaxt="n", main=row.names(Members_array_1[i,]), type="b", xlab="timepoint (hr)", ylab="Fold Change")
  axis(1, labels=timepoint_order, at=log(timepoint_order,2))
}
dev.off()

pdf("young_cluster2.pdf")
#cluster1
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(log(timepoint_order,2),km$centers[2,], xaxt="n", main="Cluster 2",type="b",xlab="timepoint (hr)", ylab="Fold change")
axis(1, labels=timepoint_order, at=log(timepoint_order,2))
for (i in 1:length(Members_array_2[,1]))
{
  plot(log(timepoint_order,2),Members_array_2[i,], xaxt="n", main=row.names(Members_array_2[i,]), type="b", xlab="timepoint (hr)", ylab="Fold Change")
  axis(1, labels=timepoint_order, at=log(timepoint_order,2))
}
dev.off()

pdf("young_cluster3.pdf")
#cluster1
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(log(timepoint_order,2),km$centers[3,], xaxt="n", main="Cluster 3",type="b",xlab="timepoint (hr)", ylab="Fold change")
axis(1, labels=timepoint_order, at=log(timepoint_order,2))
for (i in 1:length(Members_array_3[,1]))
{
  plot(log(timepoint_order,2),Members_array_3[i,], xaxt="n", main=row.names(Members_array_3[i,]), type="b", xlab="timepoint (hr)", ylab="Fold Change")
  axis(1, labels=timepoint_order, at=log(timepoint_order,2))
}
dev.off()

pdf("young_cluster4.pdf")
#cluster1
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(log(timepoint_order,2),km$centers[4,], xaxt="n", main="Cluster 4",type="b",xlab="timepoint (hr)", ylab="Fold change")
axis(1, labels=timepoint_order, at=log(timepoint_order,2))
for (i in 1:length(Members_array_4[,1]))
{
  plot(log(timepoint_order,2),Members_array_4[i,], xaxt="n", main=row.names(Members_array_4[i,]), type="b", xlab="timepoint (hr)", ylab="Fold Change")
  axis(1, labels=timepoint_order, at=log(timepoint_order,2))
}
dev.off()


pdf("young_cluster5.pdf")
#cluster1
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(log(timepoint_order,2),km$centers[5,], xaxt="n", main="Cluster 5",type="b",xlab="timepoint (hr)", ylab="Fold change")
axis(1, labels=timepoint_order, at=log(timepoint_order,2))
for (i in 1:length(Members_array_5[,1]))
{
  plot(log(timepoint_order,2),Members_array_5[i,], xaxt="n", main=row.names(Members_array_5[i,]), type="b", xlab="timepoint (hr)", ylab="Fold Change")
  axis(1, labels=timepoint_order, at=log(timepoint_order,2))
}
dev.off()

