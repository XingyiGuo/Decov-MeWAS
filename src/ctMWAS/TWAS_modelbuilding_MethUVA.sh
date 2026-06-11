#!/bin/bash
#SBATCH --job-name=UVA
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7-0:0:0

Rscript TWAS_modelbuilding_MethUVA.R $1 $2

