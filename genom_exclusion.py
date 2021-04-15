#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:30:37 2020

@author: deborah
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
in_opts = parser.add_argument_group()
in_opts.add_argument("--qcID", metavar="{File1}", type=str, required=False, default=None)
in_opts.add_argument("--genomExcl", metavar="{File2}", type=str, required=True, default=None)

#folder="/home/deborah/Documents/Crick/Projects/UKB_phenome_wide_ZNF516/UKB_ZNF516/"
#"/home/schneid/Documents/CAMP/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/"
#"/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/"

args = parser.parse_args()

#qcID=args.qcID # qc_pass.id
excl=args.genomExcl # 'ZNF516_genetic_excl.csv'
#
#Df=pd.read_csv(qcID,sep='\t',header=None)
Excl_df=pd.read_csv(excl)
#
# ethnic grouping : 1=Caucasian, NA otherwise
Excl_df.dropna(subset=['Genetic_ethnic_grouping-0.0'],inplace=True)
# recommended exclusions: 1=poor heterozygosity/missingness, NA otherwise
Excl_df= Excl_df[Excl_df['Recommended_genomic_analysis_exclusions-0.0'].isna()]
#
Df=Excl_df#Df=pd.merge(Df,Excl_df['eid'],left_on=0,right_on='eid',how='inner')
Df['FID']=Df['eid']
Df=Df[['eid','FID']]
Df.to_csv('qc_pass2.id',sep='\t',header=None,index=None)
