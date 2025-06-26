##############################################
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(Seurat)
library(tidyverse)
library(readr)
library(qs)
library(BiocParallel)
library(ComplexHeatmap)
register(MulticoreParam(workers = 8, progressbar = TRUE)) 
library(ComplexHeatmap)
library(grid)
data <- mono@meta.data

###################################################################################################################
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$celltype,
                          colname.patient = "orig.ident",
                          loc = cellInfo.tb$timepoint,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb$adj.p.value <- p.adjust(count.dist.melt.ext.tb$p.value,"BH")
  return(count.dist.melt.ext.tb)
}

OR.time.list <- do.tissueDist(cellInfo.tb=mono@meta.data,
                out.prefix="/public/users/wangrong/home/wangrong/results/new_revision/ORplot",
                pdf.width=4,pdf.height=6,verbose=1,meta.cluster = mono$celltype,loc = mono$timepoint)

#colnames(OR.time.list$OR.dist.mtx)<- factor(colnames(OR.time.list$OR.dist.mtx),levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44"))
#df %>% arrange(factor(OR_plot, levels = LETTERS[c(3, 1, 2)]))

bk <- c(seq(0,1.5,by=0.01),seq(1.6,3,by=0.01))
pheatmap(b,border="white", cluster_cols = F, cluster_rows = F,color = 
          c(colorRampPalette(c("skyblue3", "white"))(length(bk)/2),colorRampPalette(c( "white", "lightcoral"))(length(bk)/2)),fontsize_row = 10,fontsize_col = 10,cellwidth = 20,cellheight = 10,breaks=bk,legend_breaks=seq(0,3,0.5))

b$timepoints <- rownames(b)
c <- melt(b,variable.name = "timepoints")
names(c)
#abc$timepoints <- factor(abc$timepoints,levels=c("day0","day1","day3","day6","day14","day30","day31","day33","day36","day44"))
ggplot(c,aes(x=timepoints,y=OR,colour=celltype,group=celltype,fill=celltype,shape=timepoints)) +
     facet_wrap(celltype~.)+
     scale_shape_manual(values = c(20,20,20,20,20,20,20,20,20,20))+
     geom_line(size =0.8)+
     geom_point(size=3)+
     scale_color_manual(values= pal_nejm("default")(8))+
     theme(axis.line = element_line(colour = "black"))+
     labs( x = 'Cell types', y = 'OR')+
     theme_bw()+
     geom_hline(yintercept= 1.5 , linetype='dashed', color='red')+
  theme(axis.title = element_text(size=20),panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = -90,size = rel(1.2)))+
  theme(axis.text.y= element_text(size = rel(1.2)))+
  theme(strip.background = element_rect(color="black", size=1.5))