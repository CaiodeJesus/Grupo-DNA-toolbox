# Curadorias ToolBox: Busca por rsid

## 📋 Descrição

É um código que usa a API do dbSNP dataset para buscar rsids a partir do número do cromossomo e localização e vice-versa.

**query.py**: pesquisar um SNP só pelo rsid
**probable_snps.py**: completar um arquivo .txt de genotipagem com .


## 🚀 Tecnologias

Este projeto foi desenvolvido com as seguintes tecnologias:

- **Python 3.10.12**
- **API do dbSNP** para pesquisar sobre os SNPs. A documentação está disponível em clinicaltables.nlm.nih.gov/apidoc/snps/v3/doc.html
- **wget** para baixar coisas da internet
- **Numpy**  para manipulação de dados

## 📦 Como instalar e rodar o projeto

### Pré-requisitos

Antes de começar, você vai precisar ter instalado em sua máquina as seguintes ferramentas. Em parênteses temos as instruções de como baixar::

- **[Python 3.10](https://www.python.org/downloads/)**
- **[wget](pip install wget)**


### Uso query.py

exemplo de uso:
> python3 query.py
Por favor, coloque uma linha com as informações a serem buscadas da seguinte forma:
. chr_num    pos genotipo
ex:
.   1   11875   AA
Escreva aqui:.   1   11854457    AA
100% [..................................................................................] 107 / 107
A busca encontrou 2 resultados
rs191183052     rs754015864
118544572       11854457
A/C     A/G,


### Uso Probable SNPs

- o uso é simples, é só vc colocar o arquivo sem preenchimento na pasta 23andme
e digitar
	python3 probable_snps.py <nome_do_arquivo>.txt 1
se não quiser que apareçam as mensagens, é só trocar o 1 por 0

se tudo der certo, (se não der pode dar um toque e me mandar um print do erro),
é pra aparecer um arquivo <nome_do_arquivo>.log na pasta logs/ e o com os pontos substituídos na pasta 23andme_complete/


