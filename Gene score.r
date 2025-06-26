library(MuDataSeurat)
my <- ReadH5MU("Myedata_05_annotation_all.h5mu")

#inflammatory gene
INFLAMMATORY = readLines("HALLMARK_INFLAMMATORY_RESPONSE.v2024.1.Hs.gmt")
res <- strsplit(INFLAMMATORY, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
res <- lapply(res, "[", -c(1:2))
inflammatory.list <- list(res$HALLMARK_INFLAMMATORY_RESPONSE)
#antiviral gene
innate_antiviral_response_gene <- read_excel("GO_term_summary_20241212_224248.xlsx")
antiviral.list <- list(unique(toupper(innate_antiviral_response_gene$Symbol)))
#cytokine gene
cytokine <- read.table("cytokine_gene.txt",sep="\t")
cytokine.list <- list(unique(cytokine$V1))
day_inflammatory <- AddModuleScore(mono,assay="rna",features=inflammatory.list,nbin = 24,min.cells=1,ctrl = 100,name = "inflammatory_score")

timepoints <- list("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
subset <-list()
day_inflammatory_aggregate<-list()
day_antiviral_aggregate <- list()
day_cytokine_aggregate <- list()
for (i in 1:length(timepoint)){
print(timepoints[[i]])
subset[[i]]<- subset(mono, timepoint==timepoints[[i]])
print(table(subset[[i]]$timepoint))
day_inflammatory_s <- AddModuleScore(subset[[i]],assay="rna",features=inflammatory.list,nbin = 24,min.cells=1,ctrl = 100,name = "inflammatory_score")
day_antiviral_s <- AddModuleScore(subset[[i]],assay="rna",features=antiviral.list,nbin = 24,min.cells=1,ctrl = 100,name = "antiviral_score")
day_cytokine_s <- AddModuleScore(subset[[i]],assay="rna",features=cytokine.list,nbin = 24,min.cells=1,ctrl = 100,name = "cytokine_score")

day_inflammatory_aggregate[[i]]<-aggregate(day_inflammatory_s$inflammatory_score1, by=list(type=day_inflammatory_s$celltype),mean)
day_antiviral_aggregate[[i]]<-aggregate(day_antiviral_s$antiviral_score1, by=list(type=day_antiviral_s$celltype),mean)
day_cytokine_aggregate[[i]]<-aggregate(day_cytokine_s$cytokine_score1, by=list(type=day_cytokine_s$celltype),mean)
print(paste0(timepoint[[i]], " is down"))
}
names(day_inflammatory_aggregate) <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
names(day_antiviral_aggregate) <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
names(day_cytokine_aggregate) <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")

df_inflammatory <- data.table::rbindlist(day_inflammatory_aggregate,use.names=TRUE, fill=TRUE, idcol="timepoint") 
df_antiviral <- data.table::rbindlist(day_antiviral_aggregate,use.names=TRUE, fill=TRUE, idcol="timepoint")
df_cytokine <- data.table::rbindlist(day_cytokine_aggregate,use.names=TRUE, fill=TRUE, idcol="timepoint")

df1 <- spread(df_inflammatory, key = "type",value = "x", fill = 0)
library("dplyr")
timepoints <- c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
df1.1<- df1[match(timepoints, df1$time), ]   
rownames(df1.1) <- df1.1$time
df1.1 <- df1.1[,-1]%>%as.matrix
annotation_row = data.frame(timeClass = factor(rep(c("first_dose", "second_dose"), each=5 )))
rownames(annotation_row) = c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
ann_colors = list(timeClass = c(first_dose = "#add9e6", second_dose = "#962828") )

pheatmap(df1.1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
         cellwidth = 15,cellheight = 15,annotation_legend=T, main ="Inflammatory score", gaps_row = c(5),annotation_row = annotation_row, annotation_colors = ann_colors)

pheatmap(df1.1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("#552e88", "white", "#b35c29"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 15,cellheight = 15,annotation_legend=T, main ="Inflammatory score", gaps_row = c(5),annotation_row = annotation_row, annotation_colors = ann_colors)
ggsave("dose_inflammatory.score.pdf",width = 5, height = 7)

df2 <- spread(df_antiviral, key = "type",value = "x", fill = 0)
df2.1<- df2[match(timepoints, df2$timepoint), ]   
rownames(df2.1) <- df2.1$timepoint
df2.1 <- df2.1[,-1]%>%as.matrix
pheatmap(df2.1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 15,cellheight = 15,annotation_legend=T, main ="Antiviral score", gaps_row = c(5),annotation_row = annotation_row, annotation_colors = ann_colors)

df3 <- spread(df_cytokine, key = "type",value = "x", fill = 0)
df3.1<- df3[match(timepoints, df3$timepoint), ]   
rownames(df3.1) <- df3.1$timepoint
df3.1 <- df3.1[,-1]%>%as.matrix
pheatmap(df3.1,border="white", cluster_cols = F, cluster_rows = F,color =  colorRampPalette(c("skyblue3", "white", "lightcoral"))(100),fontsize_row = 14,fontsize_col = 14,
cellwidth = 15,cellheight = 15,annotation_legend=T, main ="Cytokine score", gaps_row = c(5),annotation_row = annotation_row, annotation_colors = ann_colors)