first_dose_time <- c("day1","day3","day6","day14","day30","day31","day33","day36","day44")
DA_limma_voom <- list()
for (i in 1:length(first_dose_time)){
  print(first_dose_time[i])
  peaks <- paste0("peak.matrix.",first_dose_time[i])
peaks <- mono.peak.matrix[,mono.peak.matrix$timepoint%in%c("day0",first_dose_time[i])]@assays@data$PeakMatrix
rownames(peaks) <- peak.matrix.name

meta=colData(mono.peak.matrix[,mono.peak.matrix$timepoint%in%c("day0",first_dose_time[i])])
meta$Sample <- as.factor(meta$Sample)
meta$timepoint <- factor(meta$timepoint,levels = c(first_dose_time[i],"day0"))

DA_limma_voom[[i]] = run_de(peaks, meta =meta ,replicate_col = "Sample", cell_type_col = "celltype", label_col = "timepoint",de_family = 'pseudobulk',input_type = "scATAC",de_method = "limma",de_type = "voom")
 print(paste0(first_dose_time[i], " was done!"))
}

DA_filter <- list()
for (i in 1:length(first_dose_time)){
DA_filter[[i]] <- DA_limma_voom[[i]][DA_limma_voom[[i]]$p_val_adj<= 0.1 & abs(DA_limma_voom[[i]]$avg_logFC) >= 1,]
DA_filter[[i]]$times <- first_dose_time[[i]]
colnames(DA_filter[[i]])[ 4 ]  <- 'dayn.pct'
colnames(DA_filter[[i]])[ 6 ] <- 'dayn.exp'
}

DA_filter.df <- do.call(rbind, DA_filter)
tca_df <- DA_filter.df[,c(2,3,13)]
tca_df_w <- dcast(tca_df,gene ~ times,value.var = "avg_logFC",na.rm = TRUE,fun.aggregate = mean)
tca_df_w <- cbind(day0=0,tca_df_w)
rownames(tca_df_w) <- tca_df_w$gene
 tca_df_w <- tca_df_w[,-2]
 tca_df_w[is.na(tca_df_w)] <- 0
 b <- c(1,2,4,10,3,5,6,7,8,9)
tca_df_w<- tca_df_w[b]
 tca_df_nonna <- as.matrix(tca_df_w)
tca <- timeclust(tca_df_nonna, algo = "cm", k = 5, standardize = TRUE)
  p <- timeclustplot(tca, value = "z-score(PRKM)", cols = 3)
ann_colors = c("day0" = pal_npg("nrc")(10)[1], 
               "day1" = pal_npg("nrc")(10)[2],
               "day3" = pal_npg("nrc")(10)[3], 
               "day6" = pal_npg("nrc")(10)[4],
               "day14" = pal_npg("nrc")(10)[5],
               "day30" = pal_npg("nrc")(10)[6], 
               "day31" = pal_npg("nrc")(10)[7], 
               "day33" = pal_npg("nrc")(10)[8], 
               "day36" = pal_npg("nrc")(10)[9], 
               "day44" =pal_npg("nrc")(10)[10] )
for (i in 1:5){
  
  cluster.peaks <- tca@membership[tca@cluster == i,]
  cluster.melt <- melt(cluster.peaks)
  
  ggplot(cluster.melt, aes(x = Var2, y= value, color = Var2)) +
    geom_line(aes(group = Var1), alpha = 0.1, color = "darkgray") + 
    geom_violin() + geom_sina(size = 1) + scale_color_manual(values = ann_colors) +
    theme_classic() +
    xlab("") + ylab("z-score(PRKM)") + theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45,size = rel(1.2)))
  ggsave(paste0("TC_membership_c",i,"_violinplot.svg"), width = 4, height = 3)
}
cluster.df <- as.data.frame(tca@cluster)
cluster.df$peak <- rownames(cluster.df)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cluster.df <- separate(data = cluster.df, col = peak, into = c("chr", "site"), sep = ":")
cluster.df <- separate(data = cluster.df, col = site, into = c("start", "end"), sep = "-")

c_df.bed<- list()
for (i in 1:5){
n=0+i
c_df=cluster.df[cluster.df==n,]
c_df=c_df[,-1]
rownames(c_df) <- NULL
c_df$start <- as.numeric(c_df$start)
c_df$end <- as.numeric(c_df$end)
c_df.bed[[i]] <- GRanges(seqnames = c_df$chr, ranges = IRanges(start = c_df$start, end = c_df$end))
}


peaks_ann <- c_df.bed

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks_ann, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks_ann, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci \n relative to TSS")
gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)




#####non tca
DA_filter.df <- do.call(rbind, DA_filter)
tca_df <- DA_filter.df[,c(2,3,13)]
tca_df_w <- dcast(tca_df,gene ~ times,value.var = "avg_logFC",na.rm = TRUE,fun.aggregate = mean)
tca_df_w <- cbind(day0=0,tca_df_w)
rownames(tca_df_w) <- tca_df_w$gene
 tca_df_w <- tca_df_w[,-2]
 tca_df_w[is.na(tca_df_w)] <- 0
 b <- c(1,2,4,10,3,5,6,7,8,9)
tca_df_w<- tca_df_w[b]
 tca_df_nonna <- as.matrix(tca_df_w)
tca_1 <- timeclust(tca_df_nonna, algo = "cm", k = 5, standardize = FALSE)

p <- timeclustplot(tca, value = "z-score(PRKM)", cols = 3)
for (i in 1:5){
  cluster.peaks <- tca_1@data[tca@cluster == i,]
colnames(cluster.peaks) <-c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")
cluster.melt <- melt(cluster.peaks)

  ggplot(cluster.melt, aes(x = Var2, y= value, color = Var2)) +
    geom_line(aes(group = Var1), alpha = 0.1, color = "darkgray") + 
    geom_violin() + geom_sina(size = 1) + scale_color_manual(values = ann_colors) +
    theme_classic() +
    xlab("") + ylab("z-score(PRKM)") + theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45,size = rel(1.2)))
  ggsave(paste0("DA_non_c",i,"_violinplot.pdf"), width = 4, height = 3)
}
save(tca_1,tca_df_nonna,file = "DA_non.RData")

#####p2g
proj <- addIterativeLSI(ArchRProj = motif_archr, clusterParams = list(resolution = 0.2, sampleCells = 10000,n.start = 10),saveIterations = FALSE,useMatrix = "PeakMatrix", depthCol = "nFrags",name = "LSI_ATAC")
proj <- addIterativeLSI(ArchRProj = proj, clusterParams = list(resolution = 0.2, sampleCells = 10000,n.start = 10),saveIterations = FALSE,useMatrix = "GeneExpressionMatrix", depthCol = "Gex_nUMI",varFeatures = 2500,firstSelection = "variable",binarize = FALSE,name = "LSI_RNA")
#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
proj <- addHarmony(ArchRProj = proj,reducedDims = "LSI_Combined",name = "Harmony",groupBy = "Sample")
proj_m_p2g<- addCoAccessibility(ArchRProj = proj,reducedDims = "Harmony")
cA <- getCoAccessibility(ArchRProj = proj_m_p2g,corCutOff = 0.5,resolution = 1,returnLoops = FALSE)

proj_m_p <- addPeak2GeneLinks(ArchRProj = proj_m_p2g,reducedDims = "Harmony",useMatrix = "GeneExpressionMatrix")

knnIteration = 500

p2g_new <- addPeak2GeneLinks(ArchRProj = motif_archr,reducedDims = "Harmony",useMatrix = "GeneExpressionMatrix",
knnIteration = 200)
p2g <- getPeak2GeneLinks(ArchRProj = proj_m_p,corCutOff = 0.6,resolution = 1,returnLoops = FALSE)
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2g$idxATAC]
p2g.df <- as.data.frame(p2g)

c_p2g_df <- list()
for (i in 1:5){
n=0+i
sub_cluster.df=cluster.df[cluster.df==n,]
c_gene <- intersect(unique(rownames(sub_cluster.df)),p2g.df$peakName)
c_p2g_df[[i]] <- p2g.df[p2g.df$peakName%in%c_gene,]
}
c_p2g_df[[1]]$Cluster <- "c1"
c_p2g_df[[2]]$Cluster <- "c2"
c_p2g_df[[3]]$Cluster <- "c3"
c_p2g_df[[4]]$Cluster <- "c4"
c_p2g_df[[5]]$Cluster <- "c5"
c_p2g_df_combine <- do.call(rbind,c_p2g_df)
c_p2g_df_combine %>% arrange (Cluster, Correlation) %>%
    group_by (Cluster) %>% 
    mutate (rank = rank(Correlation))

 c_p2g_gene <- list(c_p2g_df[[1]]$geneName,c_p2g_df[[2]]$geneName,c_p2g_df[[3]]$geneName,c_p2g_df[[4]]$geneName,c_p2g_df[[5]]$geneName)

 p2g_GO <-lapply(unique(c_p2g_gene), bitr,fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb="org.Hs.eg.db")

c3_go <- enrichGO(gene =p2g_GO[[3]]$ENTREZID, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)
save(DA_limma_voom,DA_filter,tca ,proj_m_p,"DA.RData")

c1_go@result$Cluster <- "c1"
c2_go@result$Cluster <- "c2"
c3_go@result$Cluster <- "c3"
c4_go@result$Cluster <- "c4"
c5_go@result$Cluster <- "c5"
co_pathway <- list(c1_go@result[1:7,],c2_go@result[1:7,],c3_go@result[1:7,],c4_go@result[1:7,],c5_go@result[1:7,])
co_pathway_df <-  do.call(rbind,co_pathway)


library(forcats)
co_pathway_df$Description <- factor(co_pathway_df$Description,levels = unique(co_pathway_df$Description))
 
"#0072B5FF", "white","#DD9D93","#BC3C29FF"


fig2_1 <- ggplot(co_pathway_df, aes(Cluster, Description)) +
  geom_point(aes(fill=-log10(p.adjust),size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="white",high="#BC3C29FF")+
  labs(x=NULL,y=NULL)+
  coord_flip()
 ggsave("tca_p2g_enrich.dotplt.pdf", width = 9, height = 5)
save(c_p2g_df,co_pathway,file = "DA_p2g_GO.Rdata")
saveRDS(fig2_1 ,file="DA_cluster_GO.rds")
counts <- mono.bed.oder[,-c(1:3)]
rownames(counts) <- peak.matrix.name
counts.mean <- t(aggregate.Matrix(t(counts), groupings = contrasts, fun = "mean")  )


p <- plotBrowserTrack(
    ArchRProj = proj_m_p, 
    groupBy = "timepoint", 
    geneSymbol = "IL1B", 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(proj_m_p),log2Norm=TRUE,pal = pal_npg("nrc")(10))