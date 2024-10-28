# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:01:50 2024

@author: DNACONSULT

Esse código funciona basicamente como uma biblioteca. Você importa ele e ele busca coisas na API.
"""


# ['chr1',1265460,'T','C']
# ['chr1',1265505,G,A
# chr1,11854457,G,A
# chr1,21785981,G,T

import wget
import os
import numpy as np



def dbsnp37_from_rsid(rsid, query_folder = 'queries/', do_print=True, return_result=True, return_log=False, exact = True):

    """Pesquisa as posições do SNP a partir do rsid no genoma de referência 37
    Parâmetros
    ----------
    rsid: string
        Rsid que se deseja pesquisar
        
    query_folder: string
        Caminho relativo da pasta onde se deseja salvar os arquivos das buscas.
    
    do_print: bool
        Variável que indica se o usuário deseja printar os resultados ou não.
    
    return_result: bool
        Variável que indica se o usuário deseja que se retorne o resultado da busca.
    
    return_log: bool
        Variável que indica se o usuário deseja que se retorne um log de saída sobre a busca.
    
    exact: bool
        Variável que indica se o usuário deseja uma busca de um rsid somente de valores de rsid iguais ao buscado.

        
    O que retorna
    ------
    snps_result: array 
        Array com os rsids encontrados.
        Retornado se return_result=True.
        
    pos_result: array
        Array com as posições dos rsids encontrados.
        Retornado se return_result=True.

    log: string
        String com o log da busca.
        Retornado se return_log=True.
        

    Referência
    ------
    API for SNPs, NIH: clinicaltables.nlm.nih.gov/apidoc/snps/v3/doc.html

    """   
    
    query_out_name = '%s.txt'%(rsid)
    if query_out_name in os.listdir(query_folder):
        os.remove(query_folder + query_out_name)
    
    url  = 'https://clinicaltables.nlm.nih.gov/api/snps/v3/search?'
    
    url+= 'terms=%s'%rsid           
    
    wget.download(url, out = query_folder+query_out_name)
    
    file = open(query_folder+query_out_name)
    query = file.readlines()
    # print(query)
    query = query[0].replace("\"","").replace("]","").replace(",","\t")
    query = query.split('[')
    file.close()
    
    
    # print(query)
    
    num_results = int(query[1].strip('\t'))
    snps_result = []
    chr_result = []
    pos_result = []
    alleles_result = []
    gene_result = []
    
    log = ''
    
    if num_results==1:

        # log+='A busca encontrou %d resultado'%num_results+'\n'
        
        print('A busca encontrou %d resultado'%num_results)
        
    else:
        # log+='A busca encontrou %d resultados'%num_results+'\n'
        print('A busca encontrou %d resultados'%num_results)
     
        
    if num_results>10:
        num_results=10
    
    # print(query)
    for i in range(num_results):
        result_i = query[-1-i].split('\t')

        # print(len(result_i), result_i)
        if len(result_i) == 1:
            continue
        elif 'rs' not in result_i[0]:
            continue
        elif 'rs' in result_i[1]:
            continue
    
        snps_result.append(result_i[0])
        chr_result.append(result_i[1])
        pos_result.append(result_i[2])
        alleles_result.append(result_i[3])
        gene_result.append(result_i[4])

        # print('passed\n')
    
    
    snps_result = np.array(snps_result)
    chr_result = np.array(chr_result)
    pos_result = np.array(pos_result)
    alleles_result = np.array(alleles_result)
    gene_result = np.array(gene_result)
    
    # print('snps_res', snps_result)
    # print('chr_res',chr_result)
    # print('pos_res', pos_result)
    # print('alleles_res', alleles_result)
    # print('gene_res',gene_result)

    if exact:
        mask = snps_result == rsid

        #if mask.sum() < 1:
        #    print('not found')
        #    return
        snps_result = snps_result[mask]
        chr_result = chr_result[mask]
        pos_result = pos_result[mask]
        alleles_result = alleles_result[mask]
        gene_result = gene_result[mask]


        # print(query)

    # print()
    # print('snps_res', snps_result)
    # print('chr_res',chr_result)
    # print('pos_res', pos_result)
    # print('alleles_res', alleles_result)
    # print('gene_res',gene_result)



    if exact:
        log += '\t'.join([rsid,snps_result[0],chr_result[0],pos_result[0],alleles_result[0],gene_result[0]])
    else:
        log+='\t'.join(snps_result)+'\n'
        log+='\t'.join(pos_result)+'\n'
        log+='\t'.join(alleles_result)+'\n'

    if do_print:
        # print()
        # print(do_print)
        print(log)

    # result = [snps_result,pos_result]
        
    if return_result and return_log:
        return snps_result,pos_result, log
    elif return_result:
        return snps_result,pos_result
    elif return_log:
        return log



def dbsnp37(search_alleles,chromossome_num, position, query_folder = 'queries/', do_print=True, return_result=True, return_log=False):

    """Pesquisa os rsids a partir da número do cromossomo, posição, e alelos a partir do rsid no genoma de referência 37
    Parâmetros
    ----------
    search_alleles: list
        Lista com os alelos possíveis para o SNP.
        
    chromossome_num: str
        Número cromossômico da variante de acordo com o GRCh37.
        
    position: int
        Posição crômossômica da busca.
        
    query_folder: string
        Caminho relativo da pasta onde se deseja salvar os arquivos das buscas.
    
    do_print: bool
        Variável que indica se o usuário deseja printar os resultados ou não.
    
    return_result: bool
        Variável que indica se o usuário deseja que se retorne o resultado da busca.
    
    return_log: bool
        Variável que indica se o usuário deseja que se retorne um log de saída sobre a busca.
    

    O que retorna
    ------
    snps_result: array 
        Array com os rsids encontrados.
        Retornado se return_result=True.
        
    pos_result: array
        Array com as posições dos rsids encontrados.
        Retornado se return_result=True.

    log: string
        String com o log da busca.
        Retornado se return_log=True.

    Referência
    ------
    API for SNPs, NIH: clinicaltables.nlm.nih.gov/apidoc/snps/v3/doc.html


    """   
    
    
    query_out_name = 'chr%s_%d.txt'%(chromossome_num,position)
    if query_out_name in os.listdir(query_folder):
        os.remove(query_folder + query_out_name)
    
    url  = 'https://clinicaltables.nlm.nih.gov/api/snps/v3/search?'
    
    if len(search_alleles)>1:
        url+= 'terms=%s&sf=37.alleles&df=rsNum,37.pos,37.alleles'%','.join(search_alleles)
    else:
        url+= 'terms=%s&sf=37.alleles&df=rsNum,37.pos,37.alleles'%search_alleles[0]
            
    url+='&q=(37.pos:['+'%d'%(position-1)+'%20TO%20'+'%d'%(position+1)+'])AND(37.chr:%s)'%chromossome_num
    
    wget.download(url, out = query_folder+query_out_name)
    
    file = open(query_folder+query_out_name)
    query = file.readlines()
    # print(query)
    query = query[0].replace("\"","").replace("]","").replace(",","\t")
    query = query.split('[')
    file.close()
    
    
    # print(query)
    
    num_results = int(query[1].strip('\t'))
    snps_result = []
    pos_result = []
    alleles_result = []
    
    log = ''
    if num_results==1:

        # log+='A busca encontrou %d resultado'%num_results+'\n'
        
        print('A busca encontrou %d resultado'%num_results)
        
    else:
        # log+='A busca encontrou %d resultados'%num_results+'\n'
        print('A busca encontrou %d resultados'%num_results)

    # print()
    print('rsid da query\tposição\tnucleotideo')


    if num_results>7:
        num_results=7
        
    for i in range(num_results):
        result_i = query[-1-i].split('\t')
        # print(result_i)
        snps_result.append(result_i[0])
        pos_result.append(result_i[1])
        alleles_result.append(','.join(result_i[2:]))
        # print(result_i)
    

            
        # print(query)
    log+=''
    for i in range(len(snps_result)):
        log+= '\t'.join([snps_result[i],pos_result[i],alleles_result[i]]) + '\n'

        
    # log+='\t'.join(snps_result)+'\n'
    # log+='\t'.join(pos_result)+'\n'
    # log+='\t'.join(alleles_result)+'\n'
    
    
    if do_print:
        # print(do_print)
        print(log)

    # result = [snps_result,pos_result]
        
    if return_result and return_log:
        return snps_result,pos_result, log
    elif return_result:
        return snps_result,pos_result
    elif return_log:
        return log





if __name__ == '__main__':

    #search_alleles = ['A']
    #chromossome_num = 1
    #position = 98352053
    print('Por favor, coloque uma linha com as informações a serem buscadas da seguinte forma:')
    print('. chr_num    pos genotipo')
    print('ex:\n.   1   11875   AA')

    # line = input("Escreva aqui:")
    line = ".	1	11854457	AA"
    line = line.split()
    # if len(line) != 4:
    #     print('Você colocou do jeito certo?\nTente novamente.')
    #     quit()

    chromossome_num = int(line[1])
    position = int(line[2])
    search_alleles = line[3]
    if search_alleles[0] == search_alleles[1]:
        search_alleles = [search_alleles[0]]
    else:
        search_alleles = [search_alleles[0],search_alleles[1]]


#    print(chromossome_num)
#    print(position)
#    print(search_alleles)

    #return_test = dbsnp37(search_alleles,chromossome_num, position, do_print=True, return_result = True, return_log = False)
    #print(return_test[0])

    #for rsid in return_test[0]:
    #    return_test = dbsnp37_from_rsid(rsid, do_print=True, return_result = True, return_log = False)

    rsid = "rs7554936"
    return_test = dbsnp37_from_rsid(rsid, query_folder = '../Out/dbsnp_api/queries/',do_print=True, return_result = True, return_log = True)


   # print(return_test)
