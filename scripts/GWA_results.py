#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 09:29:00 2020

@author: deborah
"""

import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import argparse
import statsmodels.stats.multitest as stats

parser = argparse.ArgumentParser()
in_opts = parser.add_argument_group()
in_opts.add_argument("--RegenieRes", metavar="{File1}", type=str, required=True, default=None)


args = parser.parse_args()

phenos=args.RegenieRes
name=re.findall(r'chr18_(.*).regenie',phenos)[0]##
#
rare=re.findall(r'_rare_',phenos)
#
#"/mnt/CAMP/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results/"
#"Documents/CAMP/ZNF516/ukb_step2_BT_chr18_CongHD.regenie"#
Df=pd.read_csv(phenos,sep=' ',index_col=False )
Df.drop_duplicates(subset=['GENPOS','ID'],inplace=True)
Df.dropna(subset=['LOG10P'],inplace=True)
if len(rare)==0:
    Df=Df[(Df.A1FREQ>=0.01)&(Df.A1FREQ<=0.99)]
else:
    Df=Df[(Df.A1FREQ<=0.01)|(Df.A1FREQ>=0.99)]
    name += '_rare'

N_snps=len(Df)
Df=Df[['GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ',
       'BETA', 'SE', 'LOG10P']]
thresh=-np.log10(0.05)
thresh_signif=-np.log10(0.05/N_snps)
Df['PVAL']=[np.power(10,(-1)*x) for x in Df.LOG10P]
Signif_BH=stats.multipletests(Df['PVAL'],alpha=0.05,method='fdr_bh')
Df['PVAL_adj']=Signif_BH[1]
Df['Signif']=Signif_BH[0]###np.where(Df.PVAL_adj<0.05,1,0)

#plt.plot(Df.GENPOS, Df.LOG10P, ls='', marker='.')
#plt.axhline(y=thresh, xmin=0.01, xmax=0.99, color='orchid')
#plt.axhline(y=thresh_signif, xmin=0.01, xmax=0.99,color='red')
#plt.xlabel("SNP position on Chrm 18")
#plt.ylabel("-log10 p-value")
#plt.title(name)

Df_res=Df.sort_values(by=['LOG10P'],ascending=False)
#Df_res=Df_res.head(20)

Df_res.to_csv('ZNF516_GWA_'+name+'_results.csv',index=None)
#plt.savefig('ZNF516_GWA_'+name+'_results.png')

