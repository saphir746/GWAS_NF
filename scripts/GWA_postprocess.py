#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 17:11:31 2020

@author: schneid
"""

import pandas as pd
import re
import os
import argparse

parser = argparse.ArgumentParser()
in_opts = parser.add_argument_group()
in_opts.add_argument("--Folder", metavar="{File1}", type=str, required=True, default=None)

args = parser.parse_args()

folder=args.Folder
#"/mnt/CAMP/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results/"

files=[f for f in os.listdir(folder) if ('.csv' in f) & ('ZNF516_GWA' in f)]

Df=pd.DataFrame(columns=['ID','GENPOS', 'ALLELE0', 'ALLELE1', 'A1FREQ'])
Dft=pd.DataFrame(columns=['GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE','PVAL', 'PVAL_adj','Phenotype'])
for f in files:
    name=re.findall(r'ZNF516_GWA_(.*)_results.csv',f)[0]
    df=pd.read_csv(folder+f)
    #df=df[df.PVAL_adj<0.1]#df[df['Signif']==1]
    if len(df)>0:
         Df=pd.concat([Df,df[['ID','GENPOS', 'ALLELE0', 'ALLELE1', 'A1FREQ']]])
         ####
         dfsub=df[df.PVAL_adj==min(df.PVAL_adj)]
         if len(dfsub)>1:
              dfsub=dfsub[[re.match('rs',x) is not None for x in dfsub.ID]]
              dfsub=dfsub[dfsub.PVAL==min(dfsub.PVAL)]
         dfsub=dfsub[[c for c in dfsub.columns.tolist() if c!="Signif"]]
         dfsub['Phenotype']=name
         Dft=pd.concat([Dft,dfsub],axis=0,ignore_index=True)

###
#####
#
Dft=Dft[['GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'BETA', 'SE','PVAL', 'PVAL_adj','Phenotype']]
Dft.sort_values(by='GENPOS',inplace=True)
Dft.to_csv('ZNF516_All_GWA_signif_results.csv',index=None)
#
Signif_pheno=Dft.Phenotype.tolist()
#
pics=[f for f in os.listdir(folder) if ('.png' in f) & ('ZNF516_GWA' in f)]
pics=[p for p in pics if re.findall('|'.join(Signif_pheno),p)]

with open('Signif_phenotypes_manPlot.txt', 'w') as f:
    f.writelines('\n'.join(pics))
