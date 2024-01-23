# Treinamento-Analise-dos-Dados-COVID-Longa-UFES
Esse repositório contém instruções para análise dos dados do projeto COVID Longa da Universidade Federal do Espírito Santo (UFES).
Todo o treinamento foi realizado no sistema operacional Linux na distribuição Ubuntu.

## 1. Instalação das ferramentas
### 1.1. Instalação do Anaconda
1. Baixe o Anaconda pelo link (Certifique-se de baixar a versão para Linux): https://www.anaconda.com/products/individual
2. Abra o terminal e digite os seguintes comandos:
```sh
cd Downloads/
chmod +x Anaconda3-2020.11-Linux-x86_64.sh
./Anaconda3-2020.11-Linux-x86_64.sh
```
3. Siga as instruções de instalação do Anaconda.
4. Após a instalação, feche o terminal e abra novamente.
5. Digite o seguinte comando para verificar se o Anaconda foi instalado corretamente:
```sh
conda --version
```
### 1.2 Instalação dos pacotes necessários
1. Instale os pacotes básicos:
```sh
sudo apt install g++
sudo apt install git
sudo apt install openjdk-17-jre
sudo apt install wget
sudo apt install xdg-utils
```
2. Crie e se conecte a um ambiente virtual do Anaconda:
```sh
conda create -n covid-longa
conda activate covid-longa
```
3. Instale os pacotes de bioinformática no ambiente virtual:
```sh
conda install -c bioconda bcftools
conda install -c bioconda vcftools
conda install -c bioconda plink
conda install -c bioconda admixture
```
### 1.3. Instalação do DISCVRSeq
1. Abra o terminal no diretório do projeto e digite os seguintes comandos (ou baixe o arquivo pelo mesmo link e coloque-o no diretório do projeto):
```sh
wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.62/DISCVRSeq-1.3.62.jar
```
2. Abra o terminal e digite o seguinte comando:
```sh
java -jar DISCVRSeq-1.3.62.jar
```
### 1.4. Instalação do NgsRelate
1. Abra o terminal no diretório do projeto e digite os seguintes comandos:
```sh
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/
cd ..
```

## 2. Análise dos dados
### 2.1. Configurações iniciais
1. Abra o terminal no diretório do projeto.
2. Ative o ambiente virtual Conda (se ainda não estiver ativado):
```sh
conda activate covid-longa
```
### 2.2. Merge dos arquivos VCF
1. Crie um arquivo de texto 'merge.txt' com os nomes de todos VCF:
```sh
ls *.vcf > merge.txt
```
2. Faça o merge dos arquivos VCF:
```sh
bcftools merge -l merge.txt -Oz -o merged.vcf.gz
bcftools index -t merged.vcf.gz
```
### 2.3. Análise de qualidade de variantes
1. Gere o arquivo HTML com a análise de qualidade de variantes:
```sh
java -jar DISCVSeq-1.3.62.jar VariantQC -R hg38.fa -V merged.vcf.gz -O VCF_quality.html
```
> [!IMPORTANT]
>  O arquivo FASTA (.fa) precisa ter a mesma versão do genoma da chamada de variantes

2. Abra o arquivo VCF_quality.html no navegador para visualizar a análise de qualidade de variantes.
```sh
open VCF_quality.html
```
### 2.4. Análise da Relação de Parentesco
1. Filtre as variantes:
```sh
vcftools --gzvcf merged.vcf.gz --remove-indels --maf 0.05 --minQ 20 --minDP 5 --min-alleles 2 --max-alleles 2 --hwe 1e-5 --recode --stdout | gzip -c > merged_filtered_no_indels.vcf.gz
```
2. Gere o arquivo de Relação de Parentesco:
```sh
./ngsRelate/ngsRelate -h merged_filtered_no_indels.vcf.gz -O results_relatedness.txt 
```
<!-- TODO: Adicionar o script para plotar o gráfico da Relação de Parentesco -->
3. Plote o gráfico da Relação de Parentesco:
```sh
Rscript plot_relatedness.R
```
### 2.5. Análise de Mistura Genética
1. Filtre as variantes para o ADMIXTURE:
```sh
zcat merged_filtered_no_indels.vcf.gz | grep -E “^(chr[1-9]*($'\t')*)|(^#*)” | grep -v “_alt” | grep -v “Un_” | grep -v “HLA” | grep -v “random” | grep -E -v “ID\=X” | grep -E -v “ID\=Y” | grep -E -v “ID\=M” | grep -E -v “EBV” | sed s'/chr//'g > merged_filtered_plink.vcf
```
2. Converta o arquivo VCF para o formato PLINK:
```sh
plink --vcf merged_filtered_plink.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --recode --out my_plink
plink --file my_plink --make-bed --out my_plink
plink --bfile my_plink --geno 0.90 --make-bed --out my_plink_missing
```
3. Gere o arquivo de Mistura Genética:
```sh
for K in 2 3 4 5; do admixture --cv my_plink_missing.bed $K | tee log${K}.out; done
```
<!-- TODO: Adicionar o script para plotar o gráfico da Mistura Genética -->
4. Plote o gráfico da Mistura Genética:
```sh
Rscript plot_admixture.R
```
