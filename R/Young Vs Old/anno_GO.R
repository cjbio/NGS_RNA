
library("clusterProfiler")



anno_up <- log2Up$hgnc_symbol[!log2Up$hgnc_symbol==""]
anno_down <- log2Down$hgnc_symbol[!log2Down$hgnc_symbol==""]
go_up <- enrichGO(log2Up$ensembl_gene_id,keytype = 'ENSEMBL', OrgDb = "org.Hs.eg.db", ont ="BP")
go_down <- enrichGO(log2Down$ensembl_gene_id,keytype = 'ENSEMBL', OrgDb = "org.Hs.eg.db", ont ="BP")
write.csv(anno_go, "anno_consensus_peak.csv")
View(go_up)
View(go_down)

up_skel_genes <- unlist(strsplit(go_up[1,]$geneID, split="/"))
down_skel_genes <- unlist(strsplit(go_down[1,]$geneID, split="/"))

library(annotate)
unlist(lookUp(up_skel_genes, "org.Hs.eg.db", "SYMBOL"))
unlist(lookUp(down_skel_genes, "org.Hs.eg.db", "SYMBOL"))

View(anno_go)
##########################################################
anno_up_2 <- as.data.frame(anno_up_2)[anno_up_2$`y$Fold`>1,]
anno_down_2 <- anno_down[anno_down$`y$Fold`< -1,]

genes_up_2 <- anno_up_2$SYMBOL
genes_down_2 <- anno_down_2$SYMBOL
genes_up_2 <- genes_up_2[order(genes_up_2)]
genes_down_2 <- genes_down_2[order(genes_down_2)]
genes_up_2 <- unique(genes_up_2)
genes_down_2 <- unique(genes_down_2)

write.csv(anno_up_2,"genes_up_log2fold>1_MSC0_only.csv")
write.csv(anno_down_2,"genes_down_log2fold>1_MSC0_only.csv")

go_up_2 <- enrichGO(anno_up_2$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")
go_down_2 <- enrichGO(anno_down_2$geneId, OrgDb = "org.Hs.eg.db", ont ="BP")


View(go_up_2)
View(go_up_2)
up_skel_genes_2 <- unlist(strsplit(go_up_2[1,]$geneID, split="/"))
down_skel_genes_2 <- unlist(strsplit(go_down_2[1,]$geneID, split="/"))

#library(annotate)
unlist(lookUp(up_skel_genes_2, "org.Hs.eg.db", "SYMBOL"))
unlist(lookUp(down_skel_genes_2, "org.Hs.eg.db", "SYMBOL"))

save.image("complete_analysis.RData")
load("complete_analysis_Diffbind_chIPseeker_GO.RData")

