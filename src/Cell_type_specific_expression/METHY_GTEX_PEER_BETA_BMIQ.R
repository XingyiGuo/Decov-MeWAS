#!/usr/bin/env Rscript

# conda activate peer
# 2023-10-02
# conda activate peer
# expression matrix must be n x e
# covs_pc_gender_age_num must be n x (5+1+1)
# -------------------------------------------
# SETUP
library(peer)
library(data.table)

# print usage
usage <- function() {
  cat(
    'usage: Rscript peer.R tissue_name
    author: Qing Li (liqingbioinfo@gmail.com)
    R version: 3.4+, 3.6.2 recommanded')
}

# Check input args
args <- commandArgs(trailingOnly=T)
cts <- c("Enteriendocrine", "Enterocyte","Goblet","Progenitor")
ct <- cts[as.numeric(args[1])]
print(ct)

wk="/nobackup/sbcs/keep/data/GuoLab/backup/liq17/MethExp/Meth_GTEx/"
#Read in expression matrices
mat <- as.data.frame(fread(paste0(wk,'Methy_CellTypes/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.',ct,'.qn.csv'), header=TRUE))
row_num= dim(mat)[1]
col_num= dim(mat)[2]
methy_BETA <- as.matrix(mat[,4:col_num])
rownames(methy_BETA)<-mat[,1]
class(methy_BETA) <- 'numeric'

## ----------------------------------------
#Load genotype PCs, age, gender
covs_pc_df = as.data.frame(fread(paste0(wk,'838_covariates_age_sex_5pc.csv'), header=TRUE))
covs_pc_samples = colnames(covs_pc_df)
methy_BETA_samples = colnames(methy_BETA)

common_samples_ordered <- intersect(methy_BETA_samples, covs_pc_samples)

reordered_methy_BETA <- t(methy_BETA[, common_samples_ordered]) # N X C
reordered_covs_pc <- t(covs_pc_df[, common_samples_ordered]) # N X C
colnames(reordered_covs_pc) <- covs_pc_df$SUBJID

# set number of peer factors
## https://www.nature.com/articles/s41588-022-01248-z#code-availability
if (nrow(reordered_methy_BETA) < 50) {
  NumPeerFactorS <- 5
} else if (nrow(reordered_methy_BETA) >= 100) {
  NumPeerFactorS <- 20
}
print(paste0("Samples size is ", nrow(reordered_methy_BETA), ". PEER factors according to GTEx should be ", NumPeerFactorS))

fwrite(reordered_methy_BETA, paste0(wk,'Methy_CellTypes/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.',ct,'.qn.peertool.csv'),row.names=FALSE, col.names=FALSE)
fwrite(reordered_covs_pc, paste0(wk,'838_covariates_age_sex_5pc.peertool.csv'), row.names=FALSE, col.names=FALSE)

fwrite(reordered_methy_BETA, paste0(wk,'Methy_CellTypes/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.',ct,'.qn.peertool.names.csv'),row.names=TRUE, col.names=TRUE)
fwrite(reordered_covs_pc, paste0(wk,'838_covariates_age_sex_5pc.peertool.names.csv'), row.names=TRUE, col.names=TRUE)


# #Correct PEER Factors
# model <- PEER()
# print(dim(reordered_methy_BETA))
# PEER_setPhenoMean(model, as.matrix(reordered_methy_BETA)) # Sample X Genes
# PEER_setNk(model,NumPeerFactorS)
# PEER_getNk(model)

# #Nmax_iterations = 1
# #PEER_setNmax_iterations(model,Nmax_iterations)
# #Correct gentoype PCs, age, gender
# print(dim(reordered_covs_pc))
# PEER_setCovariates(model, as.matrix(reordered_covs_pc)) # Sample X covariates
# PEER_update(model)

# ## ----------------------------------------
# # Plots
# pdf(paste0(wk,'Methy_CellTypes/',ct,'_peer.diag.qn.inv.pdf'), width=6, height=8)
# PEER_plotModel(model)
# dev.off()

# # Outputs
# factors = t(PEER_getX(model))
# fwrite(factors, file=paste0(wk,'Methy_CellTypes/',ct,'.qn.inv.PEERFactors.csv'), row.names = FALSE, col.names = TRUE,sep=",",quote=FALSE)

# #creating csv files for the residuals after accounting for the factors (residuals.csv, NxG matrix), the inferred factors (X.csv, NxK), 
# #the weights of each factor for every gene (W.csv, GxK), and the inverse variance of the weights (Alpha.csv, Kx1).

# weights = PEER_getW(model)
# precision = PEER_getAlpha(model)
# residuals = t(PEER_getResiduals(model)) # tranpose back to e x n
# if(ncol(residuals) == ncol(methy_BETA[, common_samples_ordered]) ){
	# colnames(residuals) <- colnames(methy_BETA[, common_samples_ordered])
# }
# residuals_final = cbind(mat[,1:3],residuals)

# output_filename=paste0(wk, 'Methy_CellTypes/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.',ct,'.qn.inv.peerresidual.csv')
# fwrite(residuals_final, file=output_filename, row.names = FALSE, col.names = TRUE,sep=",",quote=FALSE)
# print(paste0("PEER correction of PEER factors, gentoype PCs, age and gender is finished"))
