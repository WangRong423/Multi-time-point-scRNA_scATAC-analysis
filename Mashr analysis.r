#################first dose mashr 
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(limma)
library(edgeR)
library(mashr)
library(corrplot)


mono$participant <- substr(mono$sample,1,2)
counts <- mono@assays$rna@counts
mashr_celltype<- gsub("_",".",mono@meta.data$celltype)
metadata<-mono@meta.data
metadata$cluster_type <- mashr_celltype
metadata<-metadata[,c("sample","participant","timepoint","cluster_type")]
metadata$sample <- gsub("-","",metadata$sample)
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$participant <- as.factor(metadata$participant)
metadata$cluster_type <- as.factor(metadata$cluster_type)
metadata$sample <- as.factor(metadata$sample)
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

timepoint <- purrr::set_names(levels(sce$timepoint))
nt <- length(timepoint)
cluster_type <- purrr::set_names(levels(sce$cluster_type))
na<- length(cluster_type)
participant <- purrr::set_names(levels(sce$participant))
np<- length(participant)
sample <- purrr::set_names(levels(sce$sample))
nr<- length(sample)
table(sce$sample)				
nr_cells <- as.numeric(table(sce$sample))
mr <- match(sample, sce$sample)
eia <- data.frame(colData(sce)[mr, ], nr_cells, row.names = NULL) 
groups <- colData(sce)[, c("cluster_type","sample")]
pb1 <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
class(pb1)
dim(pb1)
pb1[1:6, 1:6]
splitf <- sapply(stringr::str_split(rownames(pb1), pattern = "_",  n = 2), `[`, 1)
library(trqwe)
pb2 <- split.data.frame(pb1, factor(splitf)) %>%lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
str(pb2)
options(width = 100)
table(sce$cluster_type, sce$sample)

get_sample_ids <- function(x){
        pb2[[x]] %>%colnames()
}
de_samples <- map(1:length(cluster_type ), get_sample_ids) %>%unlist()

samples_list <- map(1:length(cluster_type), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb2)[x], each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(cluster_type), get_cluster_ids) %>%unlist()

gg_df <- data.frame(cluster_type = de_cluster_ids,sample = de_samples)
gg_df <- left_join(gg_df, eia[, c("sample","timepoint", "participant")]) 
metadata <- gg_df %>%dplyr::select(cluster_type, sample,timepoint, participant)         

metadata$cluster_type <- as.factor(metadata$cluster_type)
clusters <- levels(metadata$cluster_type)
save(pb2,metadata, clusters, file = "2024.12.17.psudocelltype.momo.matrix.RData")

est_pi.data <- list()
for( i in 1:length(clusters)){
print(clusters[i])

cluster_metadata <- metadata[which(metadata$cluster_type == clusters[i]), ]
print(head(cluster_metadata))
rownames(cluster_metadata) <- cluster_metadata$sample
counts <- pb2[[clusters[i]]]
cluster_counts <- as.data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))

dim(cluster_counts)

condition <- factor(cluster_metadata$timepoint, levels = c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44")) 
design <- model.matrix(~0 + condition ,cluster_metadata)

y <- DGEList(cluster_counts)
y <- calcNormFactors(y)

v = voom(y, design)

timepoints.gene <- lmFit( v , design)
timepoints.Bayes.gene <- eBayes(timepoints.gene)

timepoints.condition.gene.mean <-timepoints.Bayes.gene$coefficients[, 1:10]
timepoints.condition.gene.se <- ( timepoints.Bayes.gene$stdev.unscaled * sqrt(timepoints.Bayes.gene$s2.post) )[,1:10]

colnames(timepoints.condition.gene.mean) <- colnames(timepoints.condition.gene.se) <- gsub("condition","", colnames(timepoints.condition.gene.se))
#machr
data.gene = mash_set_data(timepoints.condition.gene.mean, timepoints.condition.gene.se)
data.L.gene = mash_update_data(data.gene, ref = 1)
U.c.gene = cov_canonical(data.L.gene)
print(head(names(U.c.gene)))
mashcontrast.model.gene = mash(data.L.gene, U.c.gene  )

pdf(file=paste0('timepoints_',clusters[[i]],'.pairwise_sharing.pdf'))
corrplot(get_pairwise_sharing(mashcontrast.model.gene, factor=0.5, lfsr_thresh = 0.1) , method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude\n(< Factor of 2)', mar = c(4,0,4,0))
dev.off()
print(head( get_pm(mashcontrast.model.gene) ))

est_pi <- data.frame( Type = factor( names(get_estimated_pi(mashcontrast.model.gene) ), levels = c(names(get_estimated_pi(mashcontrast.model.gene) )) ), 
                      estimates =  get_estimated_pi(mashcontrast.model.gene)  )
est_pi.data[[i]] <- est_pi
pdf(file=paste0('/public/users/wangrong/home/wangrong/results/new_revision/mono.timepoints_',clusters[[i]],'.est_pi.pdf'))
ggplot(est_pi, aes( x = Type, y = estimates) )+
	geom_bar(stat = "identity")+
	ylab(expression( pi) )+
	theme_classic()+
	theme(text = element_text(color = "black", size = 18), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = -45, vjust = 0), plot.margin = margin(r = 10, b = 10) )+
	scale_y_continuous(expand = c(0,0))
dev.off()
                                                                                                                                                                              est_pi.data[[i]] <- est_pi
sig.DEG_ID <- get_significant_results(mashcontrast.model.gene,thresh = 0.1) 

mashResult.gene <- get_pm(mashcontrast.model.gene)
colnames(mashResult.gene) <- gsub("-day0",".log2FC",colnames(mashResult.gene))
mashResult.gene <- cbind(mashResult.gene, get_lfsr( mashcontrast.model.gene ) )
colnames(mashResult.gene) <- gsub("-day0",".lfsr",colnames(mashResult.gene))
sig.mashResults.gene<- mashResult.gene[sig.DEG_ID,]
sig.mashResults.gene<-as.data.frame(sig.mashResults.gene)
sig.mashResults.gene$total.lfsr<-rowSums(sig.mashResults.gene[, c(10,11,12,13,14,15,16,17,18)])
sig.mashResults.gene<-sig.mashResults.gene[order(sig.mashResults.gene$total.lfsr),]
sig.mashResults.gene<-sig.mashResults.gene[1:25,]

write.csv(sig.mashResults.gene,file=paste0('mono_timepoints_',clusters[[i]],'.sig.mashResults.gene.top25.csv'))


sig.mashResults.gene_melt <- data.frame( geneName = rep( rownames(sig.mashResults.gene), 9 ),
                                                    posterior.log2FC = c(sig.mashResults.gene$day1.log2FC,
                                                                        sig.mashResults.gene$day3.log2FC,
                                                                        sig.mashResults.gene$day6.log2FC,
                                                                        sig.mashResults.gene$day14.log2FC,
																		sig.mashResults.gene$day30.log2FC,
                                                                        sig.mashResults.gene$day31.log2FC,
                                                                        sig.mashResults.gene$day33.log2FC,
                                                                        sig.mashResults.gene$day36.log2FC,
																		sig.mashResults.gene$day44.log2FC),
                                                    lfsr = c( sig.mashResults.gene$day1.lfsr,
                                                              sig.mashResults.gene$day3.lfsr,
                                                              sig.mashResults.gene$day6.lfsr,
                                                              sig.mashResults.gene$day14.lfsr,
															  sig.mashResults.gene$day30.lfsr,
                                                              sig.mashResults.gene$day31.lfsr,
                                                              sig.mashResults.gene$day33.lfsr,
                                                              sig.mashResults.gene$day36.lfsr,
															  sig.mashResults.gene$day44.lfsr),
                                                    condition = c( rep("day1", nrow(sig.mashResults.gene) ),
                                                                   rep("day3", nrow(sig.mashResults.gene) ),
                                                                   rep("day6", nrow(sig.mashResults.gene) ),
                                                                   rep("day14", nrow(sig.mashResults.gene) ), 
                                                                   rep("day30", nrow(sig.mashResults.gene) ),
                                                                   rep("day31", nrow(sig.mashResults.gene) ),
                                                                   rep("day33", nrow(sig.mashResults.gene) ),
                                                                   rep("day36", nrow(sig.mashResults.gene) ),
																   rep("day44", nrow(sig.mashResults.gene) )))

sig.mashResults.gene_melt$condition <- factor(sig.mashResults.gene_melt$condition , levels = c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44") ) 
pdf(file=paste0('timepoints_',clusters[[i]],'.sig.mashResults.gene.pdf'))
ggplot(sig.mashResults.gene_melt, aes( x = condition, y = geneName, size = -log(lfsr), color = posterior.log2FC ) )+geom_point()+scale_color_gradientn( colors=c("cyan","red" ))+scale_size_continuous(breaks= -log(c(0.4,0.1,0.05,0.01,0.001)),labels=c(0.4,0.1,0.05,0.01,0.001), range = c(0.5,8),  name = "lfsr" )+geom_point(data =sig.mashResults.gene_melt[sig.mashResults.gene_melt$lfsr<0.1,], aes( x = condition, y = geneName, size = -log(lfsr) ), shape = 1, color = "black" )+ggtitle("Significant in at least one conditions")+theme_bw()+theme(
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size = 1),
                   plot.title = element_text(size = 16, hjust = 0.5),
                   axis.title=element_blank(),
                   legend.title=element_text(size=15 ),legend.text = element_text(size = 20),
                   axis.text.x = element_text(size = 15,angle = -20, vjust = 0.3 , color = "black") ,axis.text.y = element_text(size = 15, color = "black" ) )
dev.off()
print(paste0(clusters[i],'is done!'))
}

pi.data <- reduce(est_pi.data, inner_join, by="Type")
pi.names <- c("Type",clusters)
colnames(pi.data) <- pi.names
pi.data.subset <- pi.data %>%dplyr::select(Type,"CD14+ Mono","CD16+ Mono",CD56dim.CD16hi,Th2)
pi.data$Type <- factor(pi.data$Type,levels=c("day1-day0","day3-day0","day6-day0","day14-day0","day30-day0","day31-day0","day33-day0","day36-day0","day44-day0","null","identity","equal_effects","simple_het_1","simple_het_2","simple_het_3"))
pi.data.sub.melt <- melt(pi.data, variable.names="Type")
ggplot( pi.data.sub.melt, aes( x = Type, 
                               y = value,
                               fill=variable))+
	geom_bar(stat = "identity",position = 'dodge',width = 0.8,color="black",size=0.1)+
	ylab(expression( pi) )+
	theme_classic()+
	theme(text = element_text(color = "black", size = 18), axis.text = element_text(color = "black"), 
	      axis.text.x = element_text(angle = -45, vjust = 0), 
	      plot.margin = margin(r = 10, b = 10))+
  scale_fill_manual(name='celltype', values=pal_nejm("default")(8))+ 
  expand_limits(y=c(0,1))
  

 #####BTM
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
geneset <- read.gmt(file.path("/public/users/wangrong/database/wangrong/Results/0712_ATAC+RNA/downstream_analysis/BTM_for_GSEA_20131008.gmt"))
topgene <- list()
for( i in 1:length(clusters)){
  print(clusters[i])
  topgene[[i]] <- read.table(file=paste0('/public/users/wangrong/home/wangrong/results/new_revision/mono_timepoints_',clusters[[i]],'.sig.mashResults.gene.top100.csv'),sep = ",",header=TRUE)
topgene[[i]]$celltype <- clusters[i]
}
#marsh_gene <- do.call(rbind, topgene)

#####
BTM_list <- split(geneset$gene, geneset$term)
times <- c("day1","day3","day6","day14","day30","day31","day33","day36","day44")
day1_rank <- topgene[[1]] %>% arrange(desc(day1.log2FC)) %>% dplyr::select(X,day1.log2FC,day1.lfsr,celltype)
FCgenelist_day1 <- day1_rank$day1.log2FC 
names(FCgenelist_day1) <- day1_rank$X
fgseaRes_day1<- fgsea(BTM_list, stats = FCgenelist_day1,minSize = 1)
fgseaRes_day1 <- fgseaRes_day1[fgseaRes_day1$pval<0.05,]
fgseaRes_day1$times <- times[1]

day3_rank <- topgene[[1]] %>% arrange(desc(day3.log2FC)) %>% dplyr::select(X,day3.log2FC,day3.lfsr,celltype)
FCgenelist_day3 <- day3_rank$day3.log2FC 
names(FCgenelist_day3) <- day3_rank$X
fgseaRes_day3<- fgsea(BTM_list, stats = FCgenelist_day3,minSize = 1)
fgseaRes_day3 <- fgseaRes_day3[fgseaRes_day3$pval<0.05,]
fgseaRes_day3$times <- times[2]

day6_rank <- topgene[[1]] %>% arrange(desc(day6.log2FC)) %>% dplyr::select(X,day6.log2FC,day6.lfsr,celltype)
FCgenelist_day6 <- day6_rank$day6.log2FC 
names(FCgenelist_day6) <- day6_rank$X
fgseaRes_day6<- fgsea(BTM_list, stats = FCgenelist_day6,minSize = 1)
fgseaRes_day6 <- fgseaRes_day6[fgseaRes_day6$pval<0.05,]
fgseaRes_day6$times <- times[3]

day14_rank <- topgene[[1]] %>% arrange(desc(day14.log2FC)) %>% dplyr::select(X,day14.log2FC,day14.lfsr,celltype)
FCgenelist_day14 <- day14_rank$day14.log2FC 
names(FCgenelist_day14) <- day14_rank$X
fgseaRes_day14<- fgsea(BTM_list, stats = FCgenelist_day14,minSize = 1)
fgseaRes_day14 <- fgseaRes_day14[fgseaRes_day14$pval<0.05,]
fgseaRes_day14$times <- times[4]

day30_rank <- topgene[[1]] %>% arrange(desc(day30.log2FC)) %>% dplyr::select(X,day30.log2FC,day30.lfsr,celltype)
FCgenelist_day30 <- day30_rank$day30.log2FC 
names(FCgenelist_day30) <- day30_rank$X
fgseaRes_day30<- fgsea(BTM_list, stats = FCgenelist_day30,minSize = 1)
fgseaRes_day30 <- fgseaRes_day30[fgseaRes_day30$pval<0.05,]
fgseaRes_day30$times <- times[5]

day31_rank <- topgene[[1]] %>% arrange(desc(day31.log2FC)) %>% dplyr::select(X,day31.log2FC,day31.lfsr,celltype)
FCgenelist_day31 <- day31_rank$day31.log2FC 
names(FCgenelist_day31) <- day31_rank$X
fgseaRes_day31<- fgsea(BTM_list, stats = FCgenelist_day31,minSize = 1)
fgseaRes_day31 <- fgseaRes_day31[fgseaRes_day31$pval<0.05,]
fgseaRes_day31$times <- times[6]

day33_rank <- topgene[[1]] %>% arrange(desc(day33.log2FC)) %>% dplyr::select(X,day33.log2FC,day33.lfsr,celltype)
FCgenelist_day33 <- day33_rank$day33.log2FC 
names(FCgenelist_day33) <- day33_rank$X
fgseaRes_day33<- fgsea(BTM_list, stats = FCgenelist_day33,minSize = 1)
fgseaRes_day33 <- fgseaRes_day33[fgseaRes_day33$pval<0.05,]
fgseaRes_day33$times <- times[7]

day36_rank <- topgene[[1]] %>% arrange(desc(day36.log2FC)) %>% dplyr::select(X,day36.log2FC,day36.lfsr,celltype)
FCgenelist_day36 <- day36_rank$day36.log2FC 
names(FCgenelist_day36) <- day36_rank$X
fgseaRes_day36<- fgsea(BTM_list, stats = FCgenelist_day36,minSize = 1)
fgseaRes_day36 <- fgseaRes_day36[fgseaRes_day36$pval<0.05,]
fgseaRes_day36$times <- times[8]

day44_rank <- topgene[[1]] %>% arrange(desc(day44.log2FC)) %>% dplyr::select(X,day44.log2FC,day44.lfsr,celltype)
FCgenelist_day44 <- day44_rank$day44.log2FC 
names(FCgenelist_day44) <- day44_rank$X
fgseaRes_day44<- fgsea(BTM_list, stats = FCgenelist_day44,minSize = 1)
fgseaRes_day44 <- fgseaRes_day44[fgseaRes_day44$pval<0.05,]
fgseaRes_day44$times <- times[9]

fgseaRes_c1 <- rbind(fgseaRes_day1,fgseaRes_day3,fgseaRes_day6,fgseaRes_day14,fgseaRes_day30,fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c1 <- fgseaRes_c1 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c1$celltype <- clusters[1]
fgseaResTidy_c1$times <- factor(fgseaResTidy_c1$times,levels=c("day1","day3","day6","day14","day30","day31","day33","day36","day44"))
  
 fgseaResTidy_celltype <- rbind(fgseaResTidy_c1,fgseaResTidy_c2,fgseaResTidy_c3,fgseaResTidy_c4,fgseaResTidy_c5,fgseaResTidy_c6,fgseaResTidy_c7) 
 fgseaResTidy_celltype$times <- factor(fgseaResTidy_celltype$times,levels=c("day1","day3","day6","day14","day30","day31","day33","day36","day44"))
 fgseaResTidy_celltype$celltype_times <- paste(fgseaResTidy_celltype$celltype,fgseaResTidy_celltype$times,sep = "_")
fgseaResTidy_celltype_filter <-  fgseaResTidy_celltype %>% filter(!grepl('TBA', pathway))

ggplot(fgseaResTidy_celltype_filter, aes(x = celltype, y=pathway)) +
  facet_wrap(times~.)+
     geom_point(aes(size = -log10(pval), color = NES), shape =16,show.legend = TRUE)+
  scale_color_gradient2(low='#6699CC',mid= "white",high='#CC3333')+
    theme_bw()+theme_minimal()+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
save(fgseaResTidy_celltype_filter,file="/public/users/wangrong/home/wangrong/results/new_revision/BTM.fgsea.Rdata")


####################################second dose
mono_second<- subset(mono_v2,subset=timepoint %in% c("day30","day31","day33","day36","day44"))
counts <- mono_second@assays$rna@counts
mashr_celltype<- gsub("_",".",mono_second@meta.data$celltype)
metadata<-mono_second@meta.data
metadata$cluster_type <- mashr_celltype
metadata<-metadata[,c("sample","particapants","timepoint","cluster_type")]
metadata$sample <- gsub("-","",metadata$sample)
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$particapant <- as.factor(metadata$particapants)
metadata$cluster_type <- as.factor(metadata$cluster_type)
metadata$sample <- as.factor(metadata$sample)
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

timepoint <- purrr::set_names(levels(sce$timepoint))
nt <- length(timepoint)
cluster_type <- purrr::set_names(levels(sce$cluster_type))
na<- length(cluster_type)
participant <- c("p1","p2","p3","p4")
np<- length(participant)
sample <- purrr::set_names(levels(sce$sample))
nr<- length(sample)
table(sce$sample)				
nr_cells <- as.numeric(table(sce$sample))
mr <- match(sample, sce$sample)
eia <- data.frame(colData(sce)[mr, ], nr_cells, row.names = NULL) 
groups <- colData(sce)[, c("cluster_type","sample")]
pb1 <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
class(pb1)
dim(pb1)
pb1[1:6, 1:6]
splitf <- sapply(stringr::str_split(rownames(pb1), pattern = "_",  n = 2), `[`, 1)
library(trqwe)
pb2 <- split.data.frame(pb1, factor(splitf)) %>%lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
str(pb2)
options(width = 100)
table(sce$cluster_type, sce$sample)

get_sample_ids <- function(x){
        pb2[[x]] %>%colnames()
}
de_samples <- map(1:length(cluster_type ), get_sample_ids) %>%unlist()

samples_list <- map(1:length(cluster_type), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb2)[x], each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(cluster_type), get_cluster_ids) %>%unlist()

gg_df <- data.frame(cluster_type = de_cluster_ids,sample = de_samples)
gg_df <- left_join(gg_df, eia[, c("sample","timepoint", "participant")]) 
metadata <- gg_df %>%dplyr::select(cluster_type, sample,timepoint, participant)         

metadata$cluster_type <- as.factor(metadata$cluster_type)
clusters <- levels(metadata$cluster_type)


est_pi.data <- list()
for( i in 1:length(clusters)){
print(clusters[i])

cluster_metadata <- metadata[which(metadata$cluster_type == clusters[i]), ]
print(head(cluster_metadata))
rownames(cluster_metadata) <- cluster_metadata$sample
counts <- pb2[[clusters[i]]]
cluster_counts <- as.data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))

dim(cluster_counts)

condition <- factor(cluster_metadata$timepoint, levels = c("day30","day31","day33","day36","day44")) 
design <- model.matrix(~0 + condition ,cluster_metadata)

y <- DGEList(cluster_counts)
y <- calcNormFactors(y)

v = voom(y, design)

timepoints.gene <- lmFit( v , design)
timepoints.Bayes.gene <- eBayes(timepoints.gene)

timepoints.condition.gene.mean <-timepoints.Bayes.gene$coefficients[, 1:5]
timepoints.condition.gene.se <- ( timepoints.Bayes.gene$stdev.unscaled * sqrt(timepoints.Bayes.gene$s2.post) )[,1:5]

colnames(timepoints.condition.gene.mean) <- colnames(timepoints.condition.gene.se) <- gsub("condition","", colnames(timepoints.condition.gene.se))
#machr
data.gene = mash_set_data(timepoints.condition.gene.mean, timepoints.condition.gene.se)
data.L.gene = mash_update_data(data.gene, ref = 1)
U.c.gene = cov_canonical(data.L.gene)
print(head(names(U.c.gene)))
mashcontrast.model.gene = mash(data.L.gene, U.c.gene  )

pdf(file=paste0('/public/users/wangrong/home/wangrong/results/new_revision/timepoints_',clusters[[i]],'.momo_second.pairwise_sharing1.pdf'))
corrplot(get_pairwise_sharing(mashcontrast.model.gene, factor=0.5, lfsr_thresh = 0.1) , method='color', cl.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude\n(< Factor of 2)', mar = c(4,0,4,0))
dev.off()
print(head( get_pm(mashcontrast.model.gene) ))

est_pi <- data.frame( Type = factor( names(get_estimated_pi(mashcontrast.model.gene) ), levels = c(names(get_estimated_pi(mashcontrast.model.gene) )) ), 
                      estimates =  get_estimated_pi(mashcontrast.model.gene)  )
est_pi.data[[i]] <- est_pi
pdf(file=paste0('/public/users/wangrong/home/wangrong/results/new_revision/mono.timepoints_',clusters[[i]],'.momo_second.est_pi1.pdf'))
ggplot(est_pi, aes( x = Type, y = estimates) )+
	geom_bar(stat = "identity")+
	ylab(expression( pi) )+
	theme_classic()+
	theme(text = element_text(color = "black", size = 18), axis.text = element_text(color = "black"), axis.text.x = element_text(angle = -45, vjust = 0), plot.margin = margin(r = 10, b = 10) )+
	scale_y_continuous(expand = c(0,0))
dev.off()
                                                                                                                                                                              est_pi.data[[i]] <- est_pi
sig.DEG_ID <- get_significant_results(mashcontrast.model.gene,thresh = 0.1) 

mashResult.gene <- get_pm(mashcontrast.model.gene)
colnames(mashResult.gene) <- gsub("-day30",".log2FC",colnames(mashResult.gene))
mashResult.gene <- cbind(mashResult.gene, get_lfsr( mashcontrast.model.gene ) )
colnames(mashResult.gene) <- gsub("-day30",".lfsr",colnames(mashResult.gene))
sig.mashResults.gene<- mashResult.gene[sig.DEG_ID,]
sig.mashResults.gene<-as.data.frame(sig.mashResults.gene)
sig.mashResults.gene$total.lfsr<-rowSums(sig.mashResults.gene[, c(5,6,7,8)])
sig.mashResults.gene<-sig.mashResults.gene[order(sig.mashResults.gene$total.lfsr),]
sig.mashResults.gene<-sig.mashResults.gene[1:100,]

write.csv(sig.mashResults.gene,file=paste0('mono_timepoints_',clusters[[i]],'.mono_second.sig.mashResults.gene1.top100.csv'))
print(paste0(clusters[i],'is done!'))
}

pi.data <- reduce(est_pi.data, inner_join, by="Type")
pi.names <- c("Type",clusters)
colnames(pi.data) <- pi.names

pi.data$Type <- factor(pi.data$Type,levels=c("day31-day30","day33-day30","day36-day30","day44-day30","null","identity","equal_effects","simple_het_1","simple_het_2","simple_het_3"))
pi.data.sub.melt <- melt(pi.data, variable.names="Type")
ggplot( pi.data.sub.melt, aes( x = Type, 
                               y = value,
                               fill=variable))+
	geom_bar(stat = "identity",position = 'dodge',width = 0.8,color="black",size=0.1)+
	ylab(expression( pi) )+
	theme_classic()+
	theme(text = element_text(color = "black", size = 18), axis.text = element_text(color = "black"), 
	      axis.text.x = element_text(angle = -45, vjust = 0), 
	      plot.margin = margin(r = 10, b = 10))+
  scale_fill_manual(name='celltype', values=pal_nejm("default")(8))+ 
  expand_limits(y=c(0,1))


#####BTM
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
geneset <- read.gmt(file.path("BTM_for_GSEA_20131008.gmt"))
topgene <- list()
for( i in 1:length(clusters)){
  print(clusters[i])
  topgene[[i]] <- read.table(file=paste0('mono_timepoints_',clusters[[i]],'.mono_second.sig.mashResults.gene1.top200.csv'),sep = ",",header=TRUE)
topgene[[i]]$celltype <- clusters[i]
}
#marsh_gene <- do.call(rbind, topgene)

#####
day31_rank <- topgene[[1]] %>% arrange(desc(day31.log2FC)) %>% dplyr::select(X,day31.log2FC,day31.lfsr,celltype)
FCgenelist_day31 <- day31_rank$day31.log2FC 
names(FCgenelist_day31) <- day31_rank$X
fgseaRes_day31<- fgsea(BTM_list, stats = FCgenelist_day31,minSize = 1)
fgseaRes_day31 <- fgseaRes_day31[fgseaRes_day31$pval<0.1,]
fgseaRes_day31$times <- times[6]

day33_rank <- topgene[[1]] %>% arrange(desc(day33.log2FC)) %>% dplyr::select(X,day33.log2FC,day33.lfsr,celltype)
FCgenelist_day33 <- day33_rank$day33.log2FC 
names(FCgenelist_day33) <- day33_rank$X
fgseaRes_day33<- fgsea(BTM_list, stats = FCgenelist_day33,minSize = 1)
fgseaRes_day33 <- fgseaRes_day33[fgseaRes_day33$pval<0.1,]
fgseaRes_day33$times <- times[7]

day36_rank <- topgene[[1]] %>% arrange(desc(day36.log2FC)) %>% dplyr::select(X,day36.log2FC,day36.lfsr,celltype)
FCgenelist_day36 <- day36_rank$day36.log2FC 
names(FCgenelist_day36) <- day36_rank$X
fgseaRes_day36<- fgsea(BTM_list, stats = FCgenelist_day36,minSize = 1)
fgseaRes_day36 <- fgseaRes_day36[fgseaRes_day36$pval<0.1,]
fgseaRes_day36$times <- times[8]

day44_rank <- topgene[[1]] %>% arrange(desc(day44.log2FC)) %>% dplyr::select(X,day44.log2FC,day44.lfsr,celltype)
FCgenelist_day44 <- day44_rank$day44.log2FC 
names(FCgenelist_day44) <- day44_rank$X
fgseaRes_day44<- fgsea(BTM_list, stats = FCgenelist_day44,minSize = 1)
fgseaRes_day44 <- fgseaRes_day44[fgseaRes_day44$pval<0.1,]
fgseaRes_day44$times <- times[9]

fgseaRes_c1 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c1 <- fgseaRes_c1 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c1$celltype <- clusters[1]
fgseaResTidy_c1$times <- factor(fgseaResTidy_c1$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c2 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c2 <- fgseaRes_c2 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c2$celltype <- clusters[2]
fgseaResTidy_c2$times <- factor(fgseaResTidy_c2$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c3 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c3 <- fgseaRes_c3 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c3$celltype <- clusters[3]
fgseaResTidy_c3$times <- factor(fgseaResTidy_c3$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c4 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c4 <- fgseaRes_c4 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c4$celltype <- clusters[4]
fgseaResTidy_c4$times <- factor(fgseaResTidy_c4$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c5 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c5 <- fgseaRes_c5 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c5$celltype <- clusters[5]
fgseaResTidy_c5$times <- factor(fgseaResTidy_c5$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c6 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c6 <- fgseaRes_c6 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c6$celltype <- clusters[6]
fgseaResTidy_c6$times <- factor(fgseaResTidy_c6$times,levels=c("day31","day33","day36","day44"))
fgseaRes_c7 <- rbind(fgseaRes_day31,fgseaRes_day33,fgseaRes_day36,fgseaRes_day44)
fgseaResTidy_c7 <- fgseaRes_c7 %>%as_tibble() %>%arrange(desc(NES))
fgseaResTidy_c7$celltype <- clusters[7]
fgseaResTidy_c7$times <- factor(fgseaResTidy_c7$times,levels=c("day31","day33","day36","day44"))
  
fgseaResTidy_celltype <- rbind(fgseaResTidy_c1,fgseaResTidy_c2,fgseaResTidy_c3,fgseaResTidy_c4,fgseaResTidy_c5,fgseaResTidy_c6,fgseaResTidy_c7) 
 fgseaResTidy_celltype$times <- factor(fgseaResTidy_celltype$times,levels=c("day31","day33","day36","day44"))
 fgseaResTidy_celltype$celltype_times <- paste(fgseaResTidy_celltype$celltype,fgseaResTidy_celltype$times,sep = "_")
fgseaResTidy_celltype_filter <-  fgseaResTidy_celltype %>% filter(!grepl('TBA', pathway))
 
ggplot(fgseaResTidy_celltype_filter, aes(x = celltype, y=pathway)) +
     geom_point(aes(size = -log10(pval), color = NES), shape =16,show.legend = TRUE)+
  scale_color_gradient2(low='#6699CC',mid= "white",high='#CC3333')+
    theme_bw()+theme_minimal()+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
