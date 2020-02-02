
markers <- read.csv("~/rnaseq/ex0080/R/markerpanel.csv", header = FALSE)
rownames(markers) <- markers[,1]
ordered_markers <- markers
#filter nc by marker panel
nc_markers <- as.data.frame(filtered_nc_bio[rownames(filtered_nc_bio) %in% rownames(markers),])
ordered_nc_markers <- nc_markers[rownames(markers),]

ordered_nc_markers <- na.omit(ordered_nc_markers)

ordered_markers <- as.data.frame(ordered_markers[rownames(ordered_markers) %in% rownames(ordered_nc_markers),])
rownames(ordered_markers)
rownames(ordered_nc_markers)


named_nc_markers <- cbind(ordered_markers, ordered_nc_markers) 
rownames(named_nc_markers) <- named_nc_markers[,2]
named_nc_markers<- named_nc_markers[,-c(1,2)]

######################## plot the young MSCs by subject
boxcounts <- as.data.frame(t(named_nc_markers))
boxcounts <- as.data.frame(cbind(rownames(boxcounts),boxcounts[,16],condition_c))
boxcounts$Vcon1 <- condition_c
boxcounts$Vcon2 <- condition_n
boxcounts$AgeGroup <- coldata_YvO_3_clean_3$agegroup
boxcounts$treatment <- coldata_YvO_3_clean_3$treatment

boxcounts$V4 <- paste(coldata_YvO_3_clean_3$subject, boxcounts$AgeGroup, boxcounts$treatment)
boxcounts$V5 <- paste(boxcounts$AgeGroup, boxcounts$treatment)
boxcounts <- boxcounts[order(boxcounts$AgeGroup),]
boxcounts <- boxcounts[order(boxcounts$treatment),]
boxcounts <- boxcounts[order(boxcounts$Vcon2),]
####plot of markers

condition <- coldata_YvO_3_clean_3$condition
order(condition)
timepoints <- unique(as.numeric(levels(condition))[condition])
timepoint_order <- timepoints[order(timepoints)]

condition_c <- factor(coldata_YvO_3_clean_3$condition, levels = timepoint_order)
condition_n <- as.numeric(levels(condition_c)[condition_c])
pdf("boxplot_BIO+DEX_negative regulation of osteoclast differentiation and bone resorption_genes.pdf")
for(i in rownames(named_nc_markers)){
  boxcounts <- as.data.frame(t(named_nc_markers))
  boxcounts <- as.data.frame(cbind(boxcounts[,7],condition_c))
  boxcounts$V2 <- condition_c
  boxcounts$V3 <- condition_n
  boxcounts$AgeGroup <- coldata_YvO_3_clean_3$agegroup
  boxcounts$treatment <- coldata_YvO_3_clean_3$treatment
  boxcounts$V4 <- paste(boxcounts$V2, boxcounts$AgeGroup, boxcounts$treatment)
  boxcounts$V5 <- paste(boxcounts$AgeGroup, boxcounts$treatment)
  boxcounts <- boxcounts[order(boxcounts$AgeGroup),]
  boxcounts <- boxcounts[order(boxcounts$treatment),]
  colourByGroup <- sub("old dex", "grey", boxcounts$V5)
  colourByGroup <- sub("young dex", "green", colourByGroup)
  colourByGroup <- sub("old bio", "red", colourByGroup)
  colourByGroup <- sub("young bio", "blue", colourByGroup)
  print(boxplot <- ggplot(boxcounts, aes(x=log2(V3+1), y=V1, colour=V5)) +
          scale_colour_manual(values=c("red", "light grey", "blue", "light green")) +
          geom_smooth(method = lm, formula = y ~ splines::bs(x,3),fill=NA) +
          geom_point(color = colourByGroup) +
          ggtitle(i) +
          xlab("log2 time points (hr)") +
          ylab("log2 FPKM"))
  #labs(fill = "young")
}

dev.off()


# name base on the age group and timepoint/condition
tp_age <- paste(as.character(coldata_YvO_3_clean_3$treatment),
                as.character(coldata_YvO_3_clean_3$condition), 
                as.character(coldata_YvO_3_clean_3$agegroup))
colnames(named_nc_markers) <- tp_age
test1 <- named_nc_markers

# calculate the mean of the columns with the same column name
test2 <- as.data.frame(
  sapply(unique(names(test1)), 
         function(col) rowMeans(test1[names(test1) == col]) # calculate row means
  )
)

# check if the cal is correct
asd <- test1[1,]
asd[colnames(asd)=="dex 0 young"]
359.601 + 277.0566 + 433.0526
3.565700e+02*3


colnames(test2)
young_markers <- test2[,c("dex 0 young", "dex 2 young", "dex 6 young", "dex 24 young", "dex 168 young", "dex 504 young")]
old_markers <- test2[,c("dex 0 old", "dex 2 old", "dex 6 old", "dex 24 old", "dex 168 old", "dex 504 old")]

bio_young_markers <- test2[,c("dex 0 young", "bio 2 young", "bio 6 young", "bio 24 young", "bio 168 young", "bio 504 young")]
bio_old_markers <- test2[,c("dex 0 old", "bio 2 old", "bio 6 old", "bio 24 old", "bio 168 old", "bio 504 old")]

#write.csv(ordered_markers, "youngvsold_markers_deseq2.csv")


#using ed's R script

library(qPCRheatmaps)
library(pheatmap)


wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}



ordered <- cbind(young_markers,old_markers,bio_young_markers,bio_old_markers)
ordered["RUNX2",]
ordered["ALP",]
ordered["GSK3B",]
ordered["NGFR",]
#fcdexp9<-ordereddexp9_fpkm/ordereddexp9_fpkm[,1]

plot(range(as.numeric(levels(condition))), range(log2(ordered["ALP",])), type='n',xaxt='n')
lines(c(0,2,6,24,168,504), log2(young_markers["ALP",]) , type='b', col='green', pch=1)
lines(c(0,2,6,24,168,504), log2(old_markers["ALP",]) , type='b', col='red', pch=1)
lines(c(0,2,6,24,168,504), log2(bio_young_markers["ALP",]) , type='b', col='green', pch=4)
lines(c(0,2,6,24,168,504), log2(bio_old_markers["ALP",]) , type='b', col='red', pch=4)
axis(1, labels = c(0,2,6,24,168,504), at=c(0,2,6,24,168,504))
#mtext(rownames(ordered)[1])

pdf("marker_gene_exp_2.pdf")
for (i in 1:length(rownames(ordered))) {
  plot(range(as.numeric(levels(condition))), range(log2(ordered_fpkm[i,])), type='n',xaxt='n')
  lines(c(0,2,6,24,168,504), log2(young_markers[i,]) , type='b', col='green', pch=1)
  lines(c(0,2,6,24,168,504), log2(old_markers[i,]) , type='b', col='red', pch=1)
  lines(c(0,2,6,24,168,504), log2(bio_young_markers[i,]) , type='b', col='green', pch=4)
  lines(c(0,2,6,24,168,504), log2(bio_old_markers[i,]) , type='b', col='red', pch=4)
  axis(1, labels = c(0,2,6,24,168,504), at=c(0,2,6,24,168,504))
  mtext(rownames(ordered)[i])
  mtext("fpkm", side=2, line=+2)
  mtext("timepoints (hr)", side=1, line=+2)
  par(xpd=TRUE)
  legend("topright", inset=c(0,-0.25), c("young donor", "old donor", "Bio treated young", "bio treated old"), col=c("green","green", "red", "red"), pch=c(1,1,4,4))
}
dev.off()




ordered_fpkm<-ordered+1
ordered_fpkm["SP7",]
colnames(ordered_fpkm)
fc1<-ordered_fpkm[,c(1:6)]/ordered_fpkm[,1]

fc2<-ordered_fpkm[,c(7:12)]/ordered_fpkm[,7]
fc3<-ordered_fpkm[,c(13:18)]/ordered_fpkm[,13]
fc4<-ordered_fpkm[,c(19:24)]/ordered_fpkm[,19]
fc <- cbind(fc1, fc2, fc3, fc4)
fc["SP7",]

#pdf("ed_script_heatmaps_bio_dex.pdf")
#starting the heatmaps analysis
#correlation of pairwise comparison of overall gene expression between one time point and another
corr_fpkm<-cor(ordered_fpkm)
corr_fc<-cor(fc)

pheatmap(corr_fpkm^2,cluster_rows=F, cluster_cols=F, main="Correlation of FPKM values between timepoints", cellwidth=20,cellheight=20)
pheatmap(corr_fc[-1,-1]^2, cluster_cols=F,cluster_rows=F, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)

set.seed(1234)
pheatmap(ordered_fpkm, cluster_rows=F,cluster_cols=F, show_rownames=F, main="heatmap of FPKM")

pheatmap(log2(ordered_fpkm), cluster_cols=F, show_rownames=F, main="heatmap of log2 FPKM")

qPCRheatmap(log2(ordered_fpkm), cluster_rows=F, cluster_cols=F, show_rownames=T, main="Heatmap of fold changes")
qPCRheatmap(fc, limit=100, cluster_rows=F, cluster_cols=F, show_rownames=T,fontsize_row =6, main="Heatmap of fold changes")
pheatmap(fc, main="correlation of fold change relative to baseline", cellwidth=20,cellheight=20)
dev.off()

pdf("wssplot_clusters.pdf")
wssplot(ordered_fpkm[,c(1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20)], nc=50)
wssplot(fcdexp9, nc=50)
dev.off()

