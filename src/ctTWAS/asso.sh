#!/bin/bash
#SBATCH --job-name=UVA_asso
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0

sample=$1
biosoft/MetaXcan/software/SPrediXcan.py --model_db_path ${sample}.db --covariance ${sample}_cov.txt.gz --gwas_folder Cancer_GWAS_SS/CRC_GWAS_SS/MetaXcan-UK-US-180K/ --gwas_file_pattern ".*gz" --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA  --pvalue_column P --output_file  ${sample}.TWAS.csv --verbosity 1
