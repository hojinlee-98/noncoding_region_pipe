# noncoding region pipe
this repository includes scripts for analyzing noncoding regions using WGS data.

## decompose_normalization_Nov292022_hj.sh
1. decompose
split multi-allelic variants to biallelic variants. 
![image](https://user-images.githubusercontent.com/121307215/209567029-e8f992f9-af21-442e-8ff9-be6327434263.png)

2. left-normalize
left-norm expresses variants concisely.
![image](https://user-images.githubusercontent.com/121307215/209567058-0ccba6d7-ca1f-4fd9-a655-13dfe58d6e45.png)

## annovar_annotation_Nov292022_hj.sh
annotate vcf file with ANNOVAR that uses diverse databases indexed with Perl script.(the script is gaved in ANNOVAR github)

## noncoding_preprocess_Jul232022_20221119_hj.sh
the method that is established for region-based annotation use *.bed as the input file.
but there is a trivial thing we need to deal with, that is they do not give names of custom bed files to the header and info columns of vcf file.

## lcr file preprocess
download the following LCR-hs37d5.bed.gz file including low complexicity regions.
https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs37d5.bed.gz
#### LCR bed file preprocess
```shell
# 파일의 끝 부분을 확인하였더니, 기본적인 chromosome 외에 다른 것들이 존재하였다.
# 따라서 col1에 대하여 어떤 것들이 있는지 먼저 확인해 준다. 
(base) hojin@vpn1722513983 test_hail3 % tail LCR-hs37d5.bed 
hs37d5	35465535	35465560
hs37d5	35466401	35466455
hs37d5	35466559	35466595
hs37d5	35466595	35466616
hs37d5	35467526	35467558
hs37d5	35468025	35468056
hs37d5	35469953	35470048
hs37d5	35473919	35473954
hs37d5	35474399	35474448
hs37d5	35475871	35475898

# 다음을 통하여 어떠한 chromosome 이 있는지 확인하였다. 
(base) hojin@vpn1722513983 test_hail3 % awk -F '\t' '{print $1}' LCR-hs37d5.bed | sort -V -k1,1 | uniq       
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
GL000191.1
GL000192.1
GL000193.1
GL000194.1
GL000195.1
GL000196.1
GL000197.1
GL000198.1
GL000199.1
GL000200.1
GL000201.1
GL000202.1
GL000203.1
GL000204.1
GL000205.1
GL000206.1
GL000207.1
GL000208.1
GL000209.1
GL000210.1
GL000211.1
GL000212.1
GL000213.1
GL000214.1
GL000215.1
GL000216.1
GL000217.1
GL000218.1
GL000219.1
GL000220.1
GL000221.1
GL000222.1
GL000223.1
GL000224.1
GL000225.1
GL000227.1
GL000228.1
GL000229.1
GL000230.1
GL000231.1
GL000232.1
GL000233.1
GL000234.1
GL000235.1
GL000236.1
GL000237.1
GL000238.1
GL000239.1
GL000240.1
GL000241.1
GL000242.1
GL000243.1
GL000244.1
GL000245.1
GL000246.1
GL000247.1
GL000248.1
GL000249.1
NC_007605
X
Y
hs37d5

# 기본적인 chromosome (autosomal, sex) 들만 남겨주고 파일을 LCR_GRCh37_Nov182022_hj.bed에 저장하였다. 
(base) hojin@vpn1722513983 test_hail3 % cat LCR-hs37d5.bed | grep -v '^GL\|^NC\|^hs' | sort -V -k1,1 -k2,2 > LCR_GRCh37_Nov182022_hj.bed
```

### split lcr bed by chrom 
```shell
[jc2545@ruddle2 v8]$ seq 1 22 > chrom.txt ; echo -e 'X\nY' >> chrom.txt
[jc2545@c14n07 v8]$ while read line; do cat LCR_GRCh37_Nov182022_hj.bed | awk -F '\t' '{if ($1 == "'${line}'") print $0}' > ./split_bed/LCR_GRCh37_chr${line}_Nov182022_hj.bed; done < chrom.txt
[jc2545@c13n07 split_test]$ ls
LCR_GRCh37_chr10_Nov182022_hj.bed  LCR_GRCh37_chr15_Nov182022_hj.bed  LCR_GRCh37_chr1_Nov182022_hj.bed   LCR_GRCh37_chr3_Nov182022_hj.bed  LCR_GRCh37_chr8_Nov182022_hj.bed  split.sh
LCR_GRCh37_chr11_Nov182022_hj.bed  LCR_GRCh37_chr16_Nov182022_hj.bed  LCR_GRCh37_chr20_Nov182022_hj.bed  LCR_GRCh37_chr4_Nov182022_hj.bed  LCR_GRCh37_chr9_Nov182022_hj.bed  split.test.sh
LCR_GRCh37_chr12_Nov182022_hj.bed  LCR_GRCh37_chr17_Nov182022_hj.bed  LCR_GRCh37_chr21_Nov182022_hj.bed  LCR_GRCh37_chr5_Nov182022_hj.bed  LCR_GRCh37_chrX_Nov182022_hj.bed  split_bed
LCR_GRCh37_chr13_Nov182022_hj.bed  LCR_GRCh37_chr18_Nov182022_hj.bed  LCR_GRCh37_chr22_Nov182022_hj.bed  LCR_GRCh37_chr6_Nov182022_hj.bed  LCR_GRCh37_chrY_Nov182022_hj.bed
LCR_GRCh37_chr14_Nov182022_hj.bed  LCR_GRCh37_chr19_Nov182022_hj.bed  LCR_GRCh37_chr2_Nov182022_hj.bed   LCR_GRCh37_chr7_Nov182022_hj.bed  chrom.txt
[jc2545@c13n07 split_test]$ mv *.bed ./split_bed/.
```

### split vcf file by chrom
```shell
# indexing
[jc2545@c13n10 v7]$ cat indexing.sh 
#! /bin/bash
#SBATCH -J hj_indexing
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 120:00:00

module purge
module load SAMtools/1.16-GCCcore-10.2.0

bgzip -c thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf > thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf.gz
tabix -p vcf thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf.gz

# split
[jc2545@ruddle2 split_test]$ awk '{print "source ~/.bashrc; module purge; module load GATK/4.2.6.1-GCCcore-10.2.0-Java-11; gatk SelectVariants -R /home/jc2545/ref_data/h_sapiens/1000genomes/2.5/b37/human_g1k_v37_decoy.fasta -V ../../v7/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf.gz -L "$1" -O thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit_chr"$1".vcf.gz"}' < chrom.txt > split.sh
[jc2545@ruddle2 split_test]$ dsq --job-file split.sh --mem-per-cpu 100g -c 1 -t 120:00:00 -J hj_split
Batch script generated. To submit your jobs, run:
 sbatch dsq-split-2022-11-20.sh
[jc2545@ruddle2 split_test]$ sbatch dsq-split-2022-11-20.sh
```


## enhancer_promoter_burden_preprocess.splitlcr.dom.Nov212022_hj.py
domiant 
```shell
python enhancer_promoter_burden_preprocess.splitlcr.dom.Nov212022_hj.py [input_file] [lcr_bed] [MAF]
```
input_file : vcf file split by chromosome
lcr_ved : low complexity region file split by chromosome (*.bed)
MAF : minor allele frequency, this parameter is used to filter rare variants using diverse frequency. 

## enhancer_promoter_burden_preprocess.splitlcr.rec.Nov212022_hj.py
recessive
```shell
python enhancer_promoter_burden_preprocess.splitlcr.rec.Nov212022_hj.py [input_file] [lcr_bed] [MAF]
```
input_file : vcf file split by chromosome
lcr_ved : low complexity region file split by chromosome (*.bed)
MAF : minor allele frequency, this parameter is used to filter rare variants using diverse frequency. 
