# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:28:47 2024

@author: DNACONSULT


Código para buscar os SNPs na API do dbSNP

"""

import sys
sys.path.append('Support/')

import query
import os

if __name__ == '__main__':
    

    by_rsid = input('Pelo rsid? <S/N>\n')
    
    if by_rsid.lower() == 's':
        by_rsid = True
    elif by_rsid.lower() =='n':
        by_rsid = False
    else:
        print('Erro de digitação, tente novamente!')
        quit()






    

    is_many = input('Vários snps num arquivo? <S/N>\n')
    
    if is_many.lower() == 's':
        is_many = True
    elif is_many.lower() =='n':
        is_many = False
    else:
        print('Erro de digitação, tente novamente!')
        quit()


    
    if by_rsid:

    
        if is_many:
            print('Buscando SNPs listados em search_snps.txt')
            with open('search_snps.txt') as f:
                rsid_list = f.readlines()      
                
                
                
            log = 'rsid de busca\trsid da query\tchromossomo\tposição\tnucleotideo\tgene'
            print(log)
            log+='\n'
            for rsid in rsid_list[:]:
                rsid = rsid.strip('\n')
                # print(rsid)
                result = query.dbsnp37_from_rsid(rsid,
                                                 query_folder = 'Out/dbsnp_api/queries/', 
                                                 do_print=True,
                                                 return_result = False,
                                                 return_log = True)
                log+=result+'\n'
                
                    
                    
            with open('./Out/dbsnp_api/log/out_dbsnp_api.txt','w') as f:
                f.write(log)
                
                
        else:
            rsid = input('Qual o SNP que você quer pesquisar?\n')
            
            if rsid[:2] != 'rs':
                print(rsid[:2])
                quit()
                
                
                
            result = query.dbsnp37_from_rsid(rsid,
                                             query_folder = 'Out/dbsnp_api/queries/', 
                                             do_print=True,
                                             return_result = False,
                                             return_log = True)
            
            
            
            print('Query completa em ./Out/dbsnp_api/queries/%s'%(rsid))



            

        
            
    else:
        
        if is_many:

            
            file_name = input('Digite aqui o nome do arquivo:\n') #aline_santana.txt'
            
            print('Procurando em ./Data/dbsnp_api/Multifiles/')
            if file_name in os.listdir('./Data/dbsnp_api/Multifiles/'):
                print('Arquivo encontrado!')
                
                
            with open('./Data/dbsnp_api/Multifiles/'+file_name) as f:
                 rsid_list = f.readlines()    
                 
            out_name = 'out_'+file_name+'.txt'


            log = 'rsid de busca\trsid da query\tchromossomo\tposição\tnucleotideo\tgene\n'
            for line in rsid_list[:]:
                line = line.strip('\n')
                line = line.split()


                if len(line) != 4:
                    continue

                print('\t'.join(line))

                if line[1] not in ['X','Y']:
                    chromossome_num = int(line[1])
                else:
                    chromossome_num = line[1]
                    
                position = int(line[2])
                search_alleles = line[3]
                if search_alleles[0] == search_alleles[1]:
                    search_alleles = [search_alleles[0]]
                else:
                    search_alleles = [search_alleles[0],search_alleles[1]]
                
                
                result = query.dbsnp37(search_alleles,chromossome_num, position,
                                       query_folder = './Out/dbsnp_api/queries/', do_print=True, return_result = False, return_log = True)
        
                log+=result
                print('###################')
                
                
                
            with open('./Out/dbsnp_api/log/'+out_name,'w') as f:
                f.write(log)


                 
        else:           
            print('Por favor, coloque uma linha com as informações a serem buscadas da seguinte forma:')
            print('. chr_num    pos genotipo')
            print('ex:\n.   1   11875   AA')

            line = input("Escreva aqui:\n")
            
            #.   1   11875   AA
            
            
            line = line.split()
            if len(line) != 4:
                print('Você colocou do jeito certo?\nTente novamente.')
                quit()

            if line[1] not in ['X','Y']:
                chromossome_num = int(line[1])
            else:
                chromossome_num = line[1]
            
            position = int(line[2])
            search_alleles = line[3]
            if search_alleles[0] == search_alleles[1]:
                search_alleles = [search_alleles[0]]
            else:
                search_alleles = [search_alleles[0],search_alleles[1]]
            
            
            # print('rsid da query\tposição\tnucleotideo\n')
            
            return_test = query.dbsnp37(search_alleles,chromossome_num, position, 
                                        do_print=True,
                                        query_folder = './Out/dbsnp_api/queries/',
                                        return_result = True,
                                        return_log = False)

            print('Query completa em ./Out/dbsnp_api/queries/'+'chr%s_%d'%(chromossome_num,position))

    
        
