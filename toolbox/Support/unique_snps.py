# -*- coding: utf-8 -*-
"""
Created on Thu May  9 14:50:38 2024

@author: DNACONSULT


Este código é para pegar os snps únicos presentes na pasta dos arquivos preenchidos do 23andme.
Ele basicamente lê o arquivo na pasta 23andme_complete, que foi preenchido com a API do dbSNP com o código complete_snps.py


"""


import pandas as pd
import numpy as np
import os

files_23andme = os.listdir('../Out/dbsnp_api/23andme_complete')

total_snps = pd.DataFrame()

for k in range(0,len(files_23andme)):
    file_name = '../Out/dbsnp_api/23andme_complete/'+files_23andme[k]
    print(file_name)
    
    source_sheet = pd.read_csv(file_name, header=None, sep = '\t', comment='#')
    snps = source_sheet[0]
    print(len(snps))
    total_snps = pd.concat([total_snps,snps], ignore_index=True)
    
total_snps = total_snps.rename(columns={0:"rsNum"})

unique_snps = total_snps['rsNum'].unique()
unique_snps = unique_snps[unique_snps!='.']

np.savetxt('../Out/dbsnp_api/Unique_SNPs_so_far.txt', unique_snps, fmt="%s")
