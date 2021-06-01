#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 10:42:42 2021

@author: schneid
"""

import pandas as pd
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser()
in_opts = parser.add_argument_group()
in_opts.add_argument("--Phenotypes", metavar="{File1}", type=str, required=True, default=None)

#dir="/home/schneid/documents/camp/projects/swantonc/aska.przewrocka/znf516_ukb/data_dsl/ukb_znf516/regeniexnextflow/regenie_phenopre/results/"
#dir_camp="/camp/stp/babs/working/schneid/projects/projects/swantonc/aska.przewrocka/znf516_ukb/data_dsl/ukb_znf516regeniexnextflow/regenie_phenopre/results/"
#file=dir+'znf516_chd.txt'#'znf516_quant_traits.txt'
args = parser.parse_args()
File=args.Phenotypes

Df=pd.read_csv(File,sep=' ')
Vars=Df.columns.tolist()

Table=Df[Vars[2::]].describe()
Table=Table.loc[['count','mean']].T

File_name=re.match(r'.+\/(.*)$',File)[1]

def CC_ratio(X,Y):
    res=np.array([y/x for x,y in zip(X,Y)])
    res=[np.where(r>1.0,'1/'+str(int(r)),
               str(int(1/r))+'/1') for r in res]
    return res
    
if re.search('quant',File_name) is None:
    Table['Perc']=Table[['mean']]*100
    Table['mean-1']=1-Table['mean']
    Table['Case/Control']=CC_ratio(Table['mean'],Table['mean-1'])
    Table=Table[['count','Perc', 'Case/Control']]
    
Table['Phenotype']=Table.index
Table.to_csv('Summary_Table'+File_name+'.csv',index=None)