#! /bin/bash
#SBATCH --job-name=hj_noncoding
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH -t 120:00:00

cd /home/jc2545/scratch60/hj/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v8

#cat thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno.vcf | sed 's/ID=bed,/ID=gnomad_constraint_zscore,/g' | sed 's/ID=bed2,/ID=encode_CTCF_only,/g' | sed 's/ID=bed3,/ID=encode_DNase_H3K4me3,/g' | sed 's/ID=bed4,/ID=encode_PLS,/g' | sed 's/ID=bed5,/ID=encode_cCREs,/g' | sed 's/ID=bed6,/ID=encode_dELS,/g' | sed 's/ID=bed7,/ID=encode_pELS,/g' | sed 's/ID=bed8,/ID=ChromHMM_GM12878,/g' | sed 's/ID=bed9,/ID=ChromHMM_H1_hESC,/g' | sed 's/ID=bed10,/ID=ChromHMM_NHLF,/g' | sed 's/ID=bed11,/ID=vista_enhancers,/g' | sed 's/ID=bed12,/ID=fantom5_enhancers,/g' | sed 's/bed=/gnomad_constraint_zscore=/g' | sed 's/bed2=/encode_CTCF_only=/g' | sed 's/bed3=/encode_DNase_H3K4me3=/g' | sed 's/bed4=/encode_PLS=/g' | sed 's/bed5=/encode_cCREs=/g' | sed 's/bed6=/encode_dELS=/g' | sed 's/bed7=/encode_pELS=/g' | sed 's/bed8=/ChromHMM_GM12878=/g' | sed 's/bed9=/ChromHMM_H1_hESC=/g' | sed 's/bed10=/ChromHMM_NHLF=/g' | sed 's/bed11=/vista_enhancers=/g' | sed 's/bed12=/fantom5_enhancers=/g' | sed 's/Name\\x3d//g' > thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf

module purge
module load miniconda/4.12.0

unset PYTHONPATH
conda activate hail

gunzip /home/jc2545/scratch60/hj/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v7/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf.gz

python enhancer_promoter_burden_preprocess_Nov182022_hj.py /home/jc2545/scratch60/hj/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v7/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf 0.001

gzip /home/jc2545/scratch60/hj/thyroiditis/thyroiditis_noncoding/enhancer_promoter/v7/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_edit.vcf.gz

conda deactivate
