

library(ggplot2)
library(plyr)
library(dplyr)
library(stringr)

folder='/home/deborah/Documents/Crick/Projects/UKB_phenome_wide_ZNF516/'
  #'/home/schneid/Documents/CAMP/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results/'
  #'/home/deborah/Documents/CAMP_true/CAMP/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results/'
setwd(folder)
File=list.files(folder,pattern='_AP_')
#
df<-read.csv(File)
df$log10<-(-1)*log10(df$PVAL_adj)
df<-df[order(df$Category,df$PVAL_adj),]
df$idx<-1:nrow(df)

df$Phenotype<-unlist(lapply(df$Phenotype, function(x) str_replace_all(x,"_"," ")))
df$Phenotype<-unlist(lapply(df$Phenotype, function(x) str_replace_all(x,"rare","")))
df$Phenotype<-unlist(lapply(df$Phenotype, function(x) str_replace_all(x,"estimate","")))
df$Phenotype<-unlist(lapply(df$Phenotype, function(x) str_replace_all(x,"percentage","%")))

Thresh<-(-1)*log10(0.05)
border<-(-1)*log10(0.1)

nCat<-length(levels(df$Category))

axis.set <- df %>% 
  group_by(Category) %>% 
  summarize(center = mean(idx))#

mns.cat<-df %>% 
  group_by(Category) %>% 
  summarize(center = min(idx)-0.5)#

# library(extrafont)
# font_import()
# loadfonts()

sub=df[df$log10>=Thresh,]

# p<-ggplot(df,aes(x=idx,y=log10,label=Phenotype))+
#   geom_point(size = 2)+
#   geom_point(data=sub,aes(idx, log10), fill="white",colour="black",alpha=0.8,pch=21,size=2)+  
#   geom_text(hjust=-0.08, vjust=0, angle = 90, family="Helvetica", size=3)+
#   geom_hline(yintercept = Thresh,color='red')+
#   annotate("text", label = "p=0.05", x = -1, y =Thresh+0.05, color = "black")+
#   geom_hline(yintercept = border,color='grey')+
#   annotate("text", label = "p=0.1", x = -1, y =border+0.05, color = "black")+
#   coord_trans(y="log2")+
#   xlab(NULL)+ 
#   ylab('-log10 p-values (FDR adjusted)')+
#   scale_x_continuous(label = axis.set$Category, breaks = axis.set$center)+
#   geom_vline(xintercept =mns.cat$center,color="grey60", linetype = "dashed",alpha=0.5)+
#   theme_minimal()+
#   theme(axis.text.x=element_text(angle=60, family="Helvetica", hjust=1), text = element_text(size=9, family="Helvetica"))+  
#   ylim(NA, 35)
# p
# 
# #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCat))
# 
# #############
# pdf(file="GWA_results_all.pdf", height=10, width=15)
# p
# dev.off()

#############
# 
# df$OR<-exp(df$BETA)
# Keep<-c("GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","OR","SE","PVAL_adj","Phenotype","Category")
# df<-df[order(df$Category,df$PVAL_adj),]
# DT::datatable(df[,Keep])
# write.csv(df[,Keep],"GWA_results_all.csv")
