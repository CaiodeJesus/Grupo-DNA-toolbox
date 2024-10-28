# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:53:27 2024

@author: DNACONSULT
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fleury = pd.read_csv('fleury_csv.csv')

snps = np.array(fleury['SNP Name'])

bla = len(snps)*[' ']

in_array=['r','G']
for i in range(len(snps)):
    snp = snps[i]
    if snp[0] == 'r':    
        bla[i]  = snp
    elif snp[0] == 'G' and '-' in snp:
        if snp.split('-')[1][0] == 'r':
            bla[i] = snp.split('-')[1]
        
bla = np.array(bla)
bla = np.unique(bla)

with open('SNPsFleury.txt','w') as f:
    f.write('\n'.join(bla))
    
        


with open('../SNPsArray.txt') as f:
    snps_array = f.readlines()
    snps_array = [snp.strip('\n') for snp in snps_array]

#Não tem nenhum que não começa com rs
for snp in bla:
    if snp[0] != 'r':
        print(snp+'bla')

