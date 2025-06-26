library(CellChat)
library(patchwork)
all_seu <- readRDS("rna_annotation_layer3.rds")
all_seu$timepoint <- as.character(all_seu$sample)
all_seu$timepoint[str_which(all_seu$timepoint,'-10') ] <- 'day44'
all_seu$timepoint[str_which(all_seu$timepoint,'-9') ] <- 'day36'
all_seu$timepoint[str_which(all_seu$timepoint,'-8') ] <- 'day33'
all_seu$timepoint[str_which(all_seu$timepoint,'-7') ] <- 'day31'
all_seu$timepoint[str_which(all_seu$timepoint,'-6') ] <- 'day30'
all_seu$timepoint[str_which(all_seu$timepoint,'-5') ] <- 'day14'
all_seu$timepoint[str_which(all_seu$timepoint,'-4') ] <- 'day6'
all_seu$timepoint[str_which(all_seu$timepoint,'-3') ] <- 'day3'
all_seu$timepoint[str_which(all_seu$timepoint,'-2') ] <- 'day1'
all_seu$timepoint[str_which(all_seu$timepoint,'-1') ] <- 'day0'


 library(CellChat)
  CellChatDB <- CellChatDB.human
  str(CellChatDB)
  colnames(CellChatDB$interaction)
  CellChatDB$interaction[1:4,1:4]
  head(CellChatDB$complex)
  head(CellChatDB$cofactor)
  head(CellChatDB$geneInfo)
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- CellChatDB
pbmc_cellchat_proj <- list()
for (timepoint in unique(all_seu$timepoint)) {
  pbmc_subset <- all_seu[,all_seu$timepoint == timepoint]
  pbmc_cellchat <- createCellChat(pbmc_subset[["rna"]]$data)
  pbmc_meta <- data.frame(celltype = pbmc_subset$celltype_2, row.names =  Cells(pbmc_subset))
  
  pbmc_cellchat <- addMeta(pbmc_cellchat, meta = pbmc_meta, meta.name = "celltype_2")
  pbmc_cellchat <- setIdent(pbmc_cellchat, ident.use = "celltype_2") # set "labels" as default cell identity
  groupSize <- as.numeric(table(pbmc_cellchat@idents)) # number of cells in each cell group
  #############################
 
  pbmc_cellchat@DB <- CellChatDB.use
  pbmc_cellchat <- subsetData(pbmc_cellchat)
  pbmc_cellchat <- identifyOverExpressedGenes(pbmc_cellchat) 
  pbmc_cellchat <- identifyOverExpressedInteractions(pbmc_cellchat)
  #pbmc_cellchat <- smoothData(pbmc_cellchat, adj = PPI.human)
  #############################
  pbmc_cellchat <- computeCommunProb(pbmc_cellchat)###
  pbmc_cellchat <- filterCommunication(pbmc_cellchat, min.cells =10)##
  
  #######
  pbmc_cellchat <- computeCommunProbPathway(pbmc_cellchat)
  pbmc_cellchat_proj[[timepoint]] <- aggregateNet(pbmc_cellchat)
  
  # Compute the network centrality scores
  c <- netAnalysis_computeCentrality(pbmc_cellchat_proj[[timepoint]], slot.name = "netP") 
  saveRDS(pbmc_cellchat_proj[[timepoint]],paste0(timepoint,'_Cellchat.rds'))
  print(paste0(timepoint,'is complete'))
}
for (timepoint in c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")) {
pbmc_cellchat_proj[[timepoint]] <- readRDS(paste0(timepoint,'_Cellchat.rds'))
}
cellchat <- mergeCellChat(pbmc_cellchat_proj, add.names = names(pbmc_cellchat_proj))
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = pal_npg("nrc")(10),group = c(1,2,3,4,5,6,7,8,9,10))
gg2 <- compareInteractions(cellchat, show.legend = F, color.use = pal_npg("nrc")(10), group = c(1,2,3,4,5,6,7,8,9,10),measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

library(uwot)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = 'uwot')
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
rankSimilarity(cellchat, type = "functional")
cellchat@meta$celltype_2 <- factor(cellchat@meta$celltype_2 ,levels = c( "CD14+ Mono","CD16+ Mono","DC","Naïve B","Memory B","CD56bright_CD16low NK","CD56dim_CD16high NK","NKT","CD4+ T","CD8+ T"))
pairLR.use <- as.data.frame(c("CD55_ADGRE5",
                              "PTPRC_CD22",
                              "C3_ITGAM_ITGB2",
                              "PECAM1_CD38",
                              "HLA-DPA1_CD4",
                              "CD99_CD99L2",
                              'HLA-A_CD8A',
                              "ICAM1 − ITGAL"))


colnames(pairLR.use) <- 'interaction_name'
p <-netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,5:10),comparison = c(1,2,3,4,5,6,7,8,9,10), angle.x = 45,color.text = pal_npg("nrc")(10),sort.by.target = TRUE,color.heatmap = c("Spectral", "viridis"))

p2 <- annoSegment(object = p,
            annoPos = 'top', 
            aesGroup = T,
            aesGroName = 'source',
            segWidth = 0.4,
            pCol=c("#EB746A",'#26A5DF'),
            addText=T,
            textSize = 10,
            textCol = c("black","black"))

,pairLR.use = pairLR.use