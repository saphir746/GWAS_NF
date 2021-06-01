
library(dplyr)
library(ggplot2)
library(readr)
library(rtracklayer)
library(tibble)
library(tidyr)
library(DT)
library(data.table)
library(ashr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)


#dir='/home/schneid/Documents/Crick/Projects/UKB_phenome_wide_ZNF516/'
dir='/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/GWAS_NF/'
Df<-read_csv(paste0(dir,"/Interim_results_allcombined.csv"))
SNPs<-Df$ID
# annotate SNPs in biomaRt
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl.37<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                    host="grch37.ensembl.org", path="/biomart/martservice", 
                    dataset = "hsapiens_gene_ensembl") 
snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
res<-getBM(attributes = c('refsnp_id','ensembl_gene_stable_id'),
           filter=c('snp_filter'),
           values = SNPs, 
         mart = snpmart)
#
ens<-unique(res$ensembl_gene_stable_id)
ens<-ens[ens!='']
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
#
Annot=as.data.frame(cbind(ens,symbols))
Annot=merge(res,Annot,by.x='ensembl_gene_stable_id',by.y='ens',all.x=TRUE)#
Annot<-Annot[!is.na(Annot$symbols),]
#
Df.2<-merge(Df,Annot,by.x='ID',by.y='refsnp_id')
Df.2.signif<-Df.2[Df.2$BH_ajd<0.05,]
#
write.csv(Df.2.signif,paste0(dir,"Verified_annotations_results_allcombined.csv"))

#
########################
#
source(paste0(dir,"LocusZooms/functions/locus_zoom.R"))
df.plot<-Df[,c('GENPOS','ID',"PVAL")]
df.plot$CHR<-18
colnames(df.plot)<-c('BP','SNP','P','CHR')
df.plot<-df.plot[!duplicated(df.plot$SNP),]
#
res.chr18<-getBM(attributes = c('chromosome_name','ensembl_gene_id', 'external_gene_name','start_position','end_position'),
           filter=c('chromosome_name','start','end'), 
           values = list(18,73569637, 74707146), 
           mart = ensembl.37 )
###############################################################################################################
################################################################################################################
############################# run this on CAMP after loading bcftools + plink2 #################################
################################################################################################################
###############################################################################################################

locus.zoom(data = df.plot,                                    # a data.frame with columns CHR, BP, SNP, and P
           region = c(18, 73569637, 74707146),                             # the chromosome region to be included in the plot
          # offset_bp = 0,                                                  # how many basepairs around the SNP / gene / region of interest to plot
           genes.data = res.chr18,                  # a file of all the genes in the region. Cols: "Gene", "Chrom", "Start", "End".
           plot.title = "results",
           file.name = paste0(dir,"Example.jpg"))                              # the name of the file to save the plot to
