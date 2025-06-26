phenotype <- c("BNT162b2_Lym_Magri","BNT162b2_Nab_Magri_F2","BNT162b2_Nab_Magri_F2_4","BNT162b2_Nab_Magri_F4","BNT162b2_fatigue_Magri","ChAdOx1_Seroconversion_Alcalde_F","ChAdOx1_Seroconversion_Alcalde_S","ChAdOx1_severity_Alcalde","ChAdOx1_susceptibility_Alcalde")
celltype_or <- c("Mono_CD14_c01.enrichment.est","Mono_CD14_c02_LUCAT1.enrichment.est",
                 "Mono_CD14_c03_ATF.enrichment.est","Mono_CD14_c04_JUN.enrichment.est",
                 "Mono_CD14_c05_FCN1.enrichment.est","Mono_CD16_c01_PILRA.enrichment.est",
                  "Mono_CD16_c02_CCL4.enrichment.est")

for(i in 1:length(phenotype)){

path <- paste0("/public/users/wangrong/home/wangrong/results/DAP/",phenotype[i])
fileNames <- dir(path) 
filePath <- sapply(fileNames, function(x){paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){read.table(x, header=T)}) 
data_new <- lapply(data, setNames, nm = c("celltype", "max","mean","min"))
data_df[[i]] <- do.call(rbind,data_new)
data_df[[i]]$pheno <- phenotype[i]
data_df[[i]] <- data_df[[i]][-c(1,12,23,34,45,56,67,78),]
rownames(data_df[[i]]) <- NULL
}
data_df_log <- do.call(rbind,data_df)
data_df_log$time <- rep(c(" day0"," day1"," day3"," day6"," day14"," day30"," day31"," day33"," day36"," day44"),times=63)
#data_df_log$enrich <- ifelse(data_df_log$mean < 0, "NO", "YES")
data_df_file <- data_df_log[,-c(2:4)]
 paste0(data_df_file$pheno,".",data_df_file$enrich)
data_df_file <- data_df_file[,-c(2,4)]


colnames(data_df_log) <- c("celltype", "mean","min","max","pheno","time")

data_df_log$time <- factor(data_df_log$time,levels=c(" day0"," day1"," day3"," day6"," day14"," day30"," day31"," day33"," day36"," day44"))

data_df_log$pheno <- factor(data_df_log$pheno,levels=c("BNT162b2_Nab_Magri_F2","BNT162b2_Nab_Magri_F2_4","BNT162b2_Nab_Magri_F4","ChAdOx1_Seroconversion_Alcalde_F","ChAdOx1_Seroconversion_Alcalde_S","BNT162b2_fatigue_Magri","BNT162b2_Lym_Magri","ChAdOx1_severity_Alcalde","ChAdOx1_susceptibility_Alcalde"))

brewer.pal(9,"Pastel1")

pheno_colors <- c("#FBB4AE",
                  "#CCEBC5",
                  "#DECBE4",
                  "#FED9A6",
                  "#FFFFCC",
                  "#E5D8BD",
                  "#FDDAEC",
                  "#66C2A5","#B3CDE3")

data_df_log1$group <- paste0(data_df_log1$celltype,"__",data_df_log1$time)
data_df_log1$celltypes <- str_replace_all(data_df_log1$celltype,"M-[:digit:]_Mono*", "Mono")
data_df_log1$type <- str_replace_all(data_df_log1$celltypes,"M-[:digit:][:digit:]_Mono*", "Mono")

ggplot(data_df_log1,aes(x =time,y=type,group = group)) +
    geom_jjPointPie(aes(pievar = p,fill=pheno,width = 1)) +
    coord_fixed(clip = 'off') +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())+ scale_fill_manual(values = pheno_colors)
save(data_df_log,data_df_log1,file="/public/users/wangrong/home/wangrong/results/new_revision/GWAS_trait.RData")
library(ggplot2)
library(jjPlot)
library(ggnewscale)

  "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF"
pheno_color <- c("#DECBE4","#CCEBC5",pal_nejm("default")(7))
filter_torus <- read_xlsx("/public/users/wangrong/home/wangrong/results/new_revision/data_df_log1.c01.c02xlsx.xlsx")
fig5a <- filter_torus %>% ggplot(.,aes(x =`time`,y=type))+
    geom_jjPointPie(aes(pievar = p,
                        fill = pheno,
                        group=group,
                        width = 1,na.rm=TRUE),color=NA)+
    scale_fill_manual(values = pheno_colors)+
    guides(fill=guide_legend(theme = theme(legend.title = element_blank())))+
    new_scale_fill()+
    theme_test()+
    theme(axis.text.y=element_text(color="black",face="bold"),
          axis.text.x=element_text(color="black",face="bold"),
          axis.ticks = element_blank())
saveRDS(fig5a,file = "/public/users/wangrong/home/wangrong/results/new_revision/Fig5a.rds")
bar_torus <- melt(data_df_log1,id=c("pheno","enrich"),measure=c("time","type"))
###bar plot
result1 <- bar_torus %>%
  group_by(pheno, variable, value) %>%
  summarize(enrich = sum(enrich, na.rm = TRUE))%>% ungroup()
result1$value <- factor(result1$value ,levels=c(" day0"," day1"," day3"," day6"," day14"," day30"," day31"," day33"," day36"," day44"))
result1$pheno <- factor(result1$pheno,levels=c("BNT162b2_fatigue_Magri","BNT162b2_Lym_Magri","BNT162b2_Nab_Magri_F2","BNT162b2_Nab_Magri_F2_4","BNT162b2_Nab_Magri_F4","ChAdOx1_Seroconversion_Alcalde_F","ChAdOx1_Seroconversion_Alcalde_S","ChAdOx1_severity_Alcalde","ChAdOx1_susceptibility_Alcalde"))

bar_plot <- ggplot(result1, aes(x = pheno, y = enrich, fill = pheno)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Trait", y = "Count", title = "Bar plot of Mount by Trait and Category") +
  theme_minimal() +
  facet_wrap(~ value, scales = "free_y") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_fill_manual(values = pheno_colors)+
  theme_linedraw() +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),  
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    axis.text.y = element_text(color = "black"), 
    axis.text.x = element_text(color = "black"), 
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),  
    axis.ticks.x = element_blank() 
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1))
saveRDS(bar_plot,"/public/users/wangrong/home/wangrong/results/new_revision/Fig5a_2.rds")
save(data_df_log,result1,file="GWAS_figure.RData")


colnames(data_df_log) <- c("Type","Enrichment","CI_low","CI_high","Trait","Time")
data_df_log$Trait <- factor(data_df_log$Trait,levels=c("BNT162b2_Nab_Magri_F2","BNT162b2_Nab_Magri_F2_4","BNT162b2_Nab_Magri_F4","ChAdOx1_Seroconversion_Alcalde_F","ChAdOx1_Seroconversion_Alcalde_S","BNT162b2_fatigue_Magri","BNT162b2_Lym_Magri","ChAdOx1_severity_Alcalde","ChAdOx1_susceptibility_Alcalde"))

#####specific timepoints
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggsci)
library(ggforestplot)
data_df_log$Time <- as.character(data_df_log$Time)
data_df_log$Type<- gsub("M-10_","",data_df_log$Type)
data_df_log$Type<- gsub("M-[:digit:]_","",data_df_log$Type)
data_df_log_Time <- data_df_log[data_df_log$Time %in% c(' day1'," day31"),]

ggplot(data_df_log_Time, aes(x = Enrichment, 
                             y = Type, 
                             color = Type)) +
    facet_wrap(~ Time+Trait)+
    geom_vline(xintercept= 0,linetype=2 ) +geom_line(size =0.3)+
    geom_point(size = 1, 
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), 
                  width = 1, 
                  position = position_dodge(width = 0.5), size = 0.75) +
    scale_color_manual(values = colors) +
    theme_linedraw() +
    theme(
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),  
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_text(color = "black"), 
        panel.background = element_rect(fill = "white", color = NA),  
        plot.background = element_rect(fill = "white", color = NA),  
        axis.ticks.x = element_blank() 
    ) +
    labs(x = NULL)



library("scales")

#######specific celltype
library(ggplot2)
library(dplyr)
table(data_df_log$Type)
data_df_log_c02_c05 <- data_df_log[data_df_log$Type=="Mono_CD14_c02_LUCA.1" |data_df_log$Type=="Mono_CD14_c05_FC.1",]
data_df_log_c02_c05$Time <- factor(data_df_log_c02$Time,levels = c(" day0"," day1"," day3"," day6"," day14"," day30"," day31"," day33"," day36"," day44"))
c02.co5.line <- ggplot(data_df_log_c02_c05, aes(x = Time, y = Enrichment, color = Type, group = Type)) +
  facet_wrap(Trait~.,scale="free_y") +
  scale_shape_manual(values = c(20,17))+
  geom_point(size = 1) +  
  geom_line(size = 0.3) + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"))+
  labs(x = "Time", y = "Log10(OR)")+
  scale_color_nejm()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title = element_text(size=20),panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 45,size = rel(1.2)))+
  theme(axis.text.y= element_text(size = rel(1.2)))+geom_hline(aes(yintercept=0), colour="black", linetype="dashed",lwd=0.3)
c02.co5.line