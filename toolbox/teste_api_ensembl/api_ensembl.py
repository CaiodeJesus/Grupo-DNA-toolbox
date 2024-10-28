# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 11:45:55 2024

@author: DNACONSULT
"""

import pandas as pd
import numpy as np
import os
import wget
import json

rsid = input("Qual o SNP que você quer pesquisar?\n") #'rs2835370'

if __name__ == '__main__':

    if rsid not in os.listdir():
        os.mkdir(rsid)
        
    
    query_out_name = 'query_%s.txt'%rsid

    
    print('Vamos baixar!')
    
    if query_out_name in os.listdir(rsid):
        os.remove(rsid+'/'+query_out_name)
        
        
    url  = 'https://rest.ensembl.org/variation/human/%s?'%rsid
    url+= 'pops=1;content-type=application/jsonp'           
    
    try:
        wget.download(url, out = rsid+'/'+query_out_name)
    except:
        print(rsid, 'não foi encontrado pela API :(')
        os.removedirs(rsid)
        quit()

    print()
    print('Baixado!')
    
    
    print('Organizando as coisas...')
    with open(rsid+'/' + query_out_name, 'r') as file:
        lines = file.readlines()
        
        
    columns = ['allele',
               'allele_count',
               'frequency',
               'population']
    
    alleles = []
    allele_count = []
    frequency = []
    population = []
        
    arrays = [alleles,
              allele_count,
              frequency,
              population]
    
    for line in lines:
        if line[0] == ' ':
            line = line.strip(' ').split(':')
            
            for i in range(len(columns)-1):
                if line[0] == columns[i]:
                    arrays[i].append(line[1].strip('\n').strip(' '))
                    info = line[1].strip('\n').strip(' ')
                    
            if line[0] == 'population':
                info = ':'.join(line[1:]).strip('\n').strip(' ')
                population.append(info)
                if len(population) != len(allele_count):
                    allele_count.append('0')
            
    #allele_count pra inteiro
    #tirar coisas estranhas do frequency - frequency pra float
    #dividir população e subpopulação
    
        
    
    arrays = np.array(arrays).T
    df = pd.DataFrame(arrays,columns = columns)
    
    #transformando allele_count em inteiro
    df['allele_count'] = df['allele_count'].astype(int)
    
    
    #ajustando a frequencia pra float
    frequency_float = []
    for f in frequency:
        f = f.strip('&#39;')
        frequency_float.append(float(f))
    
    df['frequency'] = frequency_float
    
    
    #dividindo em populações e subpopulações
    dataset = []
    subpop = []
    
    for pop in population:
        info = pop.split(':')
        if len(info) < 2:
            dataset.append(info[0])
            subpop.append('-')
        else:
            dataset.append(info[0])
            subpop.append(':'.join(info[1:]))
            
    df.pop('population')
    df['Dataset'] = dataset
    df['Subpopulation'] = subpop
    
    print('Tudo pronto! Deve ter um .csv e o arquivo de query em uma pasta com o nome do snp.')
    df.to_csv(rsid+'/Data_from_rsid_%s.csv'%(rsid))
    
