# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA,lib.loc = "/public/users/wangrong/.conda/envs/hdWGCNA/lib/R/library/")
seurat_obj1 <- SeuratObject::UpdateSeuratObject(mono)
p <- DimPlot(seurat_obj, group.by='celltype', label=TRUE) +
   umap_theme() + ggtitle('Zhou et al Control Cortex') + NoLegend()

p
seurat_obj1 <- SetupForWGCNA(
  seurat_obj1,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj1 <- MetacellsByGroups(
  seurat_obj = seurat_obj1,
  min_cells = 45,
  group.by = c("celltype", "timepoint"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'PCA_HARMONY', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj1 <- NormalizeMetacells(seurat_obj1)

seurat_obj1 <- SetDatExpr(
  seurat_obj1,
  group_name = c("Mono_CD14_c01","Mono_CD14_c02_LUCAT1","Mono_CD14_c03_ATF3","Mono_CD14_c04_JUN","Mono_CD14_c05_FCN1","Mono_CD16_c01_PILRA","Mono_CD16_c02_CCL4"),
  group.by=c("celltype"),
  assay = 'rna', # using RNA assay
  layer = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
seurat_obj1 <- ConstructNetwork(
     seurat_obj1,
     tom_name = 'mono')

PlotDendrogram(seurat_obj1, main='mono hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)
seurat_obj1 <- ScaleData(seurat_obj1, features=VariableFeatures(seurat_obj1))
# compute all MEs in the full single-cell dataset
seurat_obj1 <- ModuleEigengenes(
 seurat_obj1,
 group.by.vars="sample")
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj1)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj1)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

seurat_obj1$timepoint<-  factor(seurat_obj1$timepoint
,levels=c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44"))

cur_traits <- c('timepoint', 'nCount_rna', 'nFeature_rna', 'total_counts_mt')
seurat_obj1$celltype_time <- paste(seurat_obj1$celltype,"-",seurat_obj1$timepoint)
seurat_obj1 <- ModuleTraitCorrelation(
  seurat_obj1,
  traits = cur_traits,
  group.by='timepoint'
)
mt_cor <- GetModuleTraitCorrelation(seurat_obj1)

names(mt_cor)

PlotModuleTraitCorrelation(
  seurat_obj1,
  label = 'fdr', # add p-val label in each cell of the heatmap
  label_symbol = 'stars', # labels as 'stars' or as 'numeric'
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#fc9272',
  mid_color = '#ffffbf',
  low_color = '#9ecae1',
  plot_max = 0.2,
  combine=T )
save(seurat_obj1,file="/public/users/wangrong/home/wangrong/results/new_revision/ConstructNetwork_wgcna.RData")

seurat_obj1 <- ModuleConnectivity(
  seurat_obj1,
  group.by = 'timepoint', 
  corFnc = "bicor", # to obtain Pearson correlation
  harmonized = TRUE,
  assay = NULL,
  slot = "data", # default to normalized 'data' slot
  group_name = 'Mono_CD14_c02_LUCAT1' )


seurat_obj1 <- ResetModuleNames(
  seurat_obj1,
  new_name = "mono_M")

p <- PlotKMEs(seurat_obj1, 
              ncol=5,
              n_hubs = 10, # number of hub genes to display
              text_size = 2,
              plot_widths = c(3, 2) # the relative width between the kME rank plot and the hub gene text
              )
p

hub_df <- GetHubGenes(seurat_obj1, n_hubs = 10)

head(hub_df)

seurat_obj1 <- ModuleExprScore(
  seurat_obj1,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
seurat_obj1 <- ModuleExprScore(
  seurat_obj1,
  n_genes = 25,
  method='UCell'
)

ModuleCorrelogram(seurat_obj1,
                  exclude_grey = TRUE, 
                  features = "hMEs" )

MEs <- GetMEs(seurat_obj1, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj1@meta.data <- cbind(seurat_obj1@meta.data, MEs)

p <- DotPlot(seurat_obj1, features=mods, group.by = 'timepoint')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p
ModuleNetworkPlot(
  seurat_obj1,
  mods = "all", # all modules are plotted.
  outdir = "/public/users/wangrong/home/wangrong/results/new_revision/", # The directory where the plots will be stored.
  plot_size = c(6, 6),
  label_center = FALSE,
  edge.alpha = 0.25,
  vertex.label.cex = 1, 
  vertex.size = 6 
)

HubGeneNetworkPlot(
  seurat_obj1,
  mods = "all", # all modules are plotted.
  n_hubs = 3, 
  n_other=6,
  edge_prop = 0.75,
)
save(seurat_obj1,"/public/users/wangrong/home/wangrong/results/new_revision/WGCNA.mono.celltype.RData")

# define the enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
seurat_obj1 <- RunEnrichr(
  seurat_obj1,
  dbs=dbs,
  max_genes = 100 # use max_genes = Inf to choose all genes
)
EnrichrDotPlot(
  seurat_obj1,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
  n_terms=2, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = FALSE # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))
