# Treinamento de Análise dos Dados do Projeto COVID Longa UFES
Esse repositório contém instruções para análise dos dados do projeto COVID Longa da Universidade Federal do Espírito Santo (UFES).

Todo o treinamento foi realizado no sistema operacional Linux, na distribuição Ubuntu.


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
2. Crie e se conecte a um ambiente virtual do Anaconda:
```sh
conda create -n covid-longa python=3.11
conda activate covid-longa
```
3. Instale os pacotes de bioinformática no ambiente virtual:
```sh
conda install -c bioconda bcftools
conda install -c bioconda vcftools
conda install -c bioconda plink
conda install -c bioconda plink2
conda install -c bioconda admixture
pip install pyplink
```
4. Instale os pacotes de ciência de dados no ambiente virtual:
```sh
conda install numpy
conda install pandas
conda install scikit-learn
conda install matplotlib
conda install rpy2
conda install -c r r-devtools
```
5. Inicialize o R:
```sh
R
```
6. Instale o pacote hdpca:
```r
devtools::install_github("cran/hdpca", dependencies = TRUE)
```
7. Saia do R:
```r
q()
```

### 1.3. Instalação do DISCVRSeq
1. Abra o terminal no diretório do projeto e digite os seguintes comandos (ou baixe o arquivo pelo mesmo link e coloque-o no diretório do projeto):
```sh
wget https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.62/DISCVRSeq-1.3.62.jar
```

### 1.4. Instalação do NgsRelate
1. Abra o terminal no diretório do projeto e digite os seguintes comandos:
```sh
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/
cd ..
```

### 1.5 Instalação do FRAPOSA
1. Abra o terminal no diretório do projeto e digite os seguintes comandos:
```sh
git clone https://github.com/daviddaiweizhang/fraposa.git
```

## 2. Análise dos dados
### 2.1. Configurações iniciais
1. Abra o terminal no diretório do projeto.
2. Ative o ambiente virtual Conda (se ainda não estiver ativado):
```sh
conda activate covid-longa
```

### 2.2. Merge dos arquivos VCF
1. Crie um arquivo de texto `merge.txt` com os nomes de todos VCF:
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
> O arquivo FASTA `.fa` precisa ter a mesma versão do genoma da chamada de variantes

2. Abra o arquivo `VCF_quality.html` no navegador para visualizar a análise de qualidade de variantes:
```sh
open VCF_quality.html
```

### 2.4. Análise da Relação de Parentesco
1. Filtre as variantes:
```sh
vcftools --gzvcf merged.vcf.gz --remove-indels --maf 0.05 --minQ 20 --minDP 5 --min-alleles 2 --max-alleles 2 --hwe 1e-5 --recode --stdout | gzip -c > merged_filtered_no_indels.vcf.gz
```
2. Gere o arquivo de Relação de Parentesco e plote o gráfico:
<!-- TODO: Adicionar o script para plotar o gráfico da Relação de Parentesco -->
```sh
./ngsRelate/ngsRelate -h merged_filtered_no_indels.vcf.gz -O results_relatedness.txt 
Rscript plot_relatedness.R
```

### 2.5. Análise de Mistura Genética
1. Filtre as variantes para o ADMIXTURE:
```sh
zcat merged_filtered_no_indels.vcf.gz | grep -E "^(chr[1-9]*($'\t')*)|(^#*)" | grep -v "_alt" | grep -v "Un_" | grep -v "HLA" | grep -v "random" | grep -E -v "ID\=X" | grep -E -v "ID\=Y" | grep -E -v "ID\=M" | grep -E -v "EBV" | sed s'/chr//'g > merged_filtered_plink.vcf
```
2. Converta o arquivo VCF para o formato BED:
```sh
plink --vcf merged_filtered_plink.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --recode --out my_plink
plink --file my_plink --make-bed --out my_plink
plink --bfile my_plink --geno 0.90 --make-bed --out my_plink_missing
```
3. Gere o arquivo de Mistura Genética e plote o gráfico:
<!-- TODO: Adicionar o script para plotar o gráfico da Mistura Genética -->
```sh
for K in 2 3 4 5; do admixture --cv my_plink_missing.bed $K | tee log${K}.out; done
Rscript plot_admixture.R
```

### 2.5 Previsão de ancestralidade
:construction: **AVISO: Essa seção atualmente está em construção e pode conter erros.** :construction:
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
3. Gere os arquivos de Previsão de ancestralidade:
```sh
./fraposa/commvar.sh thousandGenomes_renamed my_plink_missing outputThousand outputMySamples
./fraposa/fraposa_runner.py --stu_filepref outputMySamples outputThousand
./fraposa/predstupopu.py --nneighbors 20 --weights uniform outputThousand outputMySamples
```
> [!IMPORTANT]
> Se você executou o FRAPOSA anteriormente usando o mesmo conjunto de referência ou diferentes configurações de parâmetros
> você precisa excluir os arquivos intermediários `.dat`, caso contrário o FRAPOSA irá gerar um erro.
4. Plote o gráfico da Previsão de ancestralidade:
```sh	
./fraposa/plotpcs.py outputThousand outputMySamples
```