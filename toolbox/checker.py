# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:13:10 2024

@author: DNACONSULT


Código para conferir se os SNPs estão no chip, dá pra configurar pro chip da Illumina ou pro da Fleury.

"""

import numpy as np
import pandas as pd
import sys
import os





folder_gwas_files = './Data/Checker/gwas_files/'

if len(sys.argv) > 1:
    if sys.argv[1] == '-p':
        print('Execução padrão, lendo SNPs do arquivo de texto snps_to_check.txt')
        flag_p = True
        origin = 'txt'
    else:
        print('Opção inválida, ou escreve só\npython3 checker.py ou\npython3 checker.py -p')
        
else:    

    origin = input('Os SNPs são de GWAS ou de um arquivo txt? <gwas/txt>\n')
    flag_p = False

print()
if origin == 'gwas':
    name_gwas_file = input('Qual o nome do arquivo? Sem o .tsv\n')

    print('Procurando em ./Data/Checker/gwas_files/')
    if name_gwas_file+'.tsv' in os.listdir('./Data/Checker/gwas_files/'):
        print('Arquivo encontrado!')
        
        
        df_gwas = pd.read_csv(folder_gwas_files+name_gwas_file+'.tsv', sep = '\t')
        search_snps = df_gwas['SNPS']
        search_snps = np.array(search_snps)
        
        
    else:
        print('Não achei :(')
        quit()     

elif origin == 'txt':
    print('Procurando o arquivo snps_to_check.txt nesta pasta')
    if 'snps_to_check.txt' in os.listdir():
        print('Arquivo encontrado!')
        search_snps = np.loadtxt('snps_to_check.txt', dtype = str)
    else:
        print('Não achei :(')
        quit()        

else:
    print('Tente novamente! Escreva "gwas" se quiser um arquivo saído do GWAS Catalog; ou "txt" se os seus SNPs de busca estão em snps_to_check.txt')
    quit()
    
    
print()   

if not flag_p:
    fleury = input('Chip da Fleury? <S/N>\n')
    
    if fleury.lower() == 's':
        fleury = True
    elif fleury.lower() =='n':
        fleury = False
    else:
        print('Erro de digitação, tente novamente!')
        quit()
else:
    fleury = False




if fleury:
    snps_array_source = './Data/Checker/fleury/SNPsFleury.txt'    
else:
    snps_array_source = './Data/Checker/SNPsArray.txt'

  
    
snps_in_chip = np.loadtxt(snps_array_source, dtype = str)

search_snps = np.unique(search_snps)

print('\nTotal SNPs')
print('\n'.join(search_snps))
print()


count = 0
search_in_chip = ''
for snp in search_snps:
    if snp in snps_in_chip:
        print(snp, 'está no chip!')
        count+=1
        search_in_chip+=snp+'\n'
    else:
        print(snp,'não está no chip.')
    
if fleury:
    print('\nSNPs no chip da Fleury')
else:
    print('\nSNPs no chip da Eurofins')
    
print(search_in_chip)

file_result = open('./result.txt', 'w')
file_result.write(search_in_chip)
file_result.close()

print('%.d dos %.d SNPs de procura tão no chip, o que significa %.2f%%.'%(count, len(search_snps),100*count/len(search_snps) ))

print('SNPs presentes no chip foram escritos em result.txt')
