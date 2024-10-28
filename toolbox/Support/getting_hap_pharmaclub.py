# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:19:15 2024

@author: DNACONSULT

Este é um código para separar os haplótipos que estamos usando no pharmaclub,
para uma análise de frequências.

"""

import numpy as np
import pandas as pd
import os

folder_drugs = '../../geradorFarmaco/PGX/queries/drugs'
out_folder = '../Out/pharmaclub/'

drugs = os.listdir(folder_drugs)

list_of_haplotypes = []
drug_name = []

empty_folders = []
for drug in drugs[:]:
    
    haplotypes_name = 'cpic_clean.csv'
    
    if 'cpic_clean.csv' in os.listdir(folder_drugs+'/'+drug):
        haplotypes = pd.read_csv(folder_drugs+'/'+drug+'/cpic_clean.csv')
        list_of_haplotypes.append(haplotypes[["Haplotypes","rsid","nucleotide","Gene"]])
        
        drug_name += len(haplotypes)*[drug]
        
        
    else:
        print("mv drugs/'%s'/ empty/'%s'/"%(drug, drug))
        empty_folders.append(drug)
        # (folder_drugs+'/'+drug,'../../geradorFarmaco/PGX/empty/drugs')
        continue
    
mega_haplotypes = pd.concat(list_of_haplotypes)


hap_drug_pairs = pd.DataFrame()
hap_drug_pairs['Haplotypes'] = mega_haplotypes['Haplotypes']
hap_drug_pairs['Gene'] = mega_haplotypes['Gene']
hap_drug_pairs['drug'] = drug_name


mega_haplotypes.drop_duplicates(inplace=True)
mega_haplotypes.to_csv(out_folder + 'All_the_haplotypes.csv', index=False)


unique_hap = np.array(mega_haplotypes['Haplotypes'])

drugs_that_use_this_hap = []

for hap in unique_hap:
    query = hap_drug_pairs.query("Haplotypes == '%s'"%hap)
    drugs_that_use_this_hap.append(','.join(query['drug']))
    
    
hap_drug_clean = pd.DataFrame()
hap_drug_clean['Haplotype'] = unique_hap
hap_drug_clean['Drugs that use this haplotype'] = drugs_that_use_this_hap

hap_drug_clean.to_csv(out_folder + 'Haplotype_drug_pairs.csv')

print('Feito!')
print('Out files in ',out_folder)



