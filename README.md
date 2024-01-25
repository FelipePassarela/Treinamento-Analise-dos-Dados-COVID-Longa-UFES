# Treinamento de Análise dos Dados do Projeto COVID Longa UFES
Esse repositório contém instruções para análise dos dados do projeto COVID Longa da Universidade Federal do Espírito Santo (UFES).

O treinamento abrange a instalação e uso de várias ferramentas de bioinformática e ciência de dados, incluindo Conda, bcftools, vcftools, PLINK, Admixture, e bibliotecas Python como numpy, pandas, scikit-learn, matplotlib e rpy2. 

Todo o treinamento foi realizado no sistema operacional Linux, em específico na distribuição Ubuntu.
<br>

<br>

# 1. Instalação das ferramentas
## 1.1. Instalação do Conda

Conda é um gerenciador de pacotes e sistema de gerenciamento de ambiente de código aberto, que inclui uma grande variedade de ferramentas úteis para análise de dados e bioinformática. O Anaconda é uma distribuição do Conda que inclui o Python e mais de 1500 outros pacotes. Para instalá-lo, siga os seguintes passos:

1. Baixe o Conda pelo link (Certifique-se de baixar a versão para Linux): https://www.anaconda.com/products/individual
2. Abra o terminal e digite os seguintes comandos:
```sh
cd Downloads/
chmod +x Anaconda3-2020.11-Linux-x86_64.sh
./Anaconda3-2020.11-Linux-x86_64.sh
```
3. Siga as instruções de instalação do Conda.
4. Após a instalação, feche o terminal e abra novamente.
5. Digite o seguinte comando para verificar se o Conda foi instalado corretamente:
```sh
conda --version
```

## 1.2 Instalação dos pacotes necessários

Depois de instalar o Conda, você precisará instalar todas as dependências necessárias para executar as ferramentas e scripts de análise de dados. Para isso, siga os passos abaixo:

### 1.2.1 Instale os pacotes básicos
Digite os seguintes comandos no terminal:
```sh
sudo apt update
sudo apt install g++
sudo apt install git
sudo apt install openjdk-17-jre
sudo apt install wget
sudo apt install xdg-utils
sudo apt install python3
sudo apt install python3-pip
python3 -m pip install --upgrade pip
```
### 1.2.2 Ambiente virtual do Conda
1. Crie e se conecte a um ambiente virtual Conda:
```sh
conda create -n covid-longa python=3.11
conda activate covid-longa
```
2. Instale os pacotes de bioinformática no ambiente virtual:
```sh
conda install -c bioconda bcftools
conda install -c bioconda vcftools
conda install -c bioconda plink
conda install -c bioconda plink2
conda install -c bioconda admixture
pip install pyplink
```
3. Instale os pacotes de ciência de dados no ambiente virtual:
```sh
conda install numpy
conda install pandas
conda install scikit-learn
conda install matplotlib
conda install rpy2
conda install -c r r-devtools
```
4. Inicialize o R:
```sh
R
```
5. Instale o pacote `hdpca`:
```r
devtools::install_github("cran/hdpca", dependencies = TRUE)
```
6. Saia do R:
```r
q()
```

## 1.3 Instalação das ferramentas de análise de dados e bioinformática

Algumas ferramentas científicas necessárias não estão disponíveis no Conda, portanto, você precisará instalá-las manualmente. Para isso, siga os passos abaixo:

### 1.3.1 Instalação do DISCVRSeq
Abra o terminal no diretório do projeto e digite os seguintes comandos (ou baixe o arquivo pelo mesmo link e coloque-o no diretório do projeto):
```sh
wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.62/DISCVRSeq-1.3.62.jar
```
### 1.3.2 Instalação do NgsRelate
Abra o terminal no diretório do projeto e digite os seguintes comandos:
```sh
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/
cd ..
```
### 1.3.3 Instalação do FRAPOSA
Abra o terminal no diretório do projeto e digite os seguintes comandos:
```sh
git clone https://github.com/daviddaiweizhang/fraposa.git
```
<br>

# 2. Análise dos dados

Depois de instalar todas as ferramentas e dependências necessárias, você pode iniciar a análise dos dados.

## 2.1 Configurações iniciais

Primeiro, ative o ambiente virtual Conda (se não estiver ativo) e faça o merge dos arquivos VCF em um único só.

### 2.1.1 Ativação do ambiente virtual Conda
Abra o terminal no diretório do projeto e digite o seguinte comando:
```sh
conda activate covid-longa
```

### 2.1.2 Merge dos arquivos VCF
1. Crie um arquivo de texto `merge.txt` com os nomes de todos VCF:
```sh
ls *.vcf > merge.txt
```
2. Faça o merge dos arquivos VCF:
```sh
bcftools merge -l merge.txt -Oz -o merged.vcf.gz
bcftools index -t merged.vcf.gz
```
<br>

## 2.2. Execução dos programas de análise de amostras

Agora você pode executar os programas de análise de amostras. Para isso, siga os passos abaixo:

### 2.2.1 Análise de qualidade das amostras
1. Gere o arquivo HTML `VCF_quality.html` com a análise de qualidade de variantes e abra no navegador para visualiza-lo:
```sh
java -jar DISCVSeq-1.3.62.jar VariantQC -R hg38.fa -V merged.vcf.gz -O VCF_quality.html
open VCF_quality.html
```
> [!IMPORTANT]
> O arquivo FASTA `.fa` precisa ter a mesma versão do genoma da chamada de variantes

### 2.2.2 Análise da relação de parentesco
1. Filtre as variantes:
```sh
vcftools --gzvcf merged.vcf.gz --remove-indels --maf 0.05 --minQ 20 --minDP 5 --min-alleles 2 --max-alleles 2 --hwe 1e-5 --recode --stdout | gzip -c > merged_filtered_no_indels.vcf.gz
```
2. Gere o arquivo de relação de parentesco e plote o gráfico:
<!-- TODO: Adicionar o script para plotar o gráfico da Relação de Parentesco -->
```sh
./ngsRelate/ngsRelate -h merged_filtered_no_indels.vcf.gz -O results_relatedness.txt 
Rscript plot_relatedness.R
```

### 2.2.3 Análise de mistura genética
1. Filtre as variantes para o ADMIXTURE:
```sh
zcat merged_filtered_no_indels.vcf.gz | grep -E "^(chr[1-9]*($'\t')*)|(^#*)" | grep -v "_alt" | grep -v "Un_" | grep -v "HLA" | grep -v "random" | grep -E -v "ID\=X" | grep -E -v "ID\=Y" | grep -E -v "ID\=M" | grep -E -v "EBV" | sed s'/chr//'g > merged_filtered_plink.vcf
```
2. Converta o arquivo VCF para o formato BED:
```sh
plink --vcf merged_filtered_plink.vcf --double-id --allow-extra-chr --set-missing-var-ids @: --indep-pairwise 50 10 0.1 --recode --out my_plink
plink --file my_plink --make-bed --out my_plink
plink --bfile my_plink --geno 0.90 --make-bed --out my_plink_missing
```
3. Gere o arquivo de mistura genética e plote o gráfico:
<!-- TODO: Adicionar o script para plotar o gráfico da Mistura Genética -->
```sh
for K in 2 3 4 5; do admixture --cv my_plink_missing.bed $K | tee log${K}.out; done
Rscript plot_admixture.R
```

### 2.2.4 Análise de ancestralidade
1. Baixe o genoma de referência com o `wget` (ou acesse https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3, baixe os três arquivos manualmente e renomeie-os para `all_hg38.psam`, `all_hg38.pgen.zst` e `all_hg38.pvar.zst`):
```sh
wget -O all_hg38.psam "https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1"
wget -O all_hg38.pgen.zst "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
wget -O all_hg38.pvar.zst "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"
```
2. Aplique os filtros de qualidade e gere os arquivos necessários para o FRAPOSA:
```sh
plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
plink2 --pfile all_hg38 vzs --chr 1-22 --output-chr 26 --max-alleles 2 --rm-dup exclude-mismatch --set-all-var-ids @:# --make-pgen --out all_hg38_autosomes
plink2 --pfile all_hg38 vzs --min-alleles 2 --max-alleles 2 --allow-extra-chr --maf 0.005 --make-bed --out all_hg38_autosomes
plink2 --bfile all_hg38_autosomes allow-extra-chr --set-all-var-ids @:# --make-bed --out thousandGenomes_renamed
cat all_hg38.psam | awk '{print $1"\t"$1"\t"$5}' > outputThousand.popu
```
3. Gere os arquivos de análise de ancestralidade:
```sh
./fraposa/commvar.sh thousandGenomes_renamed my_plink_missing outputThousand outputMySamples
./fraposa/fraposa_runner.py --stu_filepref outputMySamples outputThousand
./fraposa/predstupopu.py --nneighbors 20 --weights uniform outputThousand outputMySamples
```
> [!IMPORTANT]
> Se você executou o FRAPOSA anteriormente usando o mesmo conjunto de referência ou diferentes configurações de parâmetros,
> você precisa excluir os arquivos intermediários `.dat`, caso contrário o FRAPOSA irá gerar um erro.
4. Plote o gráfico da análise de ancestralidade:
```sh	
./fraposa/plotpcs.py outputThousand outputMySamples
```
