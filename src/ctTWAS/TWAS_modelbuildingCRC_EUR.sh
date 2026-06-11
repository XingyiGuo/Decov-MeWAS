#!/bin/bash
#SBATCH --job-name=TWAS_EUR
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-0:0:0

Rscript TWAS_modelbuildingCRC_EUR.R $1 $2
