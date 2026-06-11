#!/usr/bin/env Rscript

# Usage: Rscript METHY_PEER.R <dataset> <cell_type_index> <working_dir>
#
# <dataset>         : "UVA" or "GTEx"
# <cell_type_index> : 1=Enteriendocrine, 2=Enterocyte, 3=Goblet, 4=Progenitor
# <working_dir>     : path to working directory
#
# Prepares PEER input files for cell-type-specific methylation beta values.
# PEER factor selection follows GTEx methylation guidelines:
#   https://www.nature.com/articles/s41588-022-01248-z
#   N < 50 -> 5, N >= 100 -> 20

library(peer)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
dataset  <- args[1]
cts      <- c("Enteriendocrine", "Enterocyte", "Goblet", "Progenitor")
ct       <- cts[as.numeric(args[2])]
wk       <- args[3]

if (!dataset %in% c("UVA", "GTEx")) stop("dataset must be 'UVA' or 'GTEx'")
print(paste0("Dataset: ", dataset, " | Cell type: ", ct))

# ── Dataset-specific file name patterns ──────────────────────────────────────
if (dataset == "UVA") {
  methy_prefix <- "CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38"
  cov_file     <- "146_covariates_Normal_Mucosa.csv"
  id_col       <- "SUBJID"
} else {
  methy_prefix <- "GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38"
  cov_file     <- "838_covariates_age_sex_5pc.csv"
  id_col       <- "SUBJID"
}

# ── Load methylation matrix ──────────────────────────────────────────────────
mat <- as.data.frame(fread(paste0(wk, "Methy_CellTypes/", methy_prefix, ".", ct, ".qn.csv"), header = TRUE))
col_num <- dim(mat)[2]
methy_BETA <- as.matrix(mat[, 4:col_num])
rownames(methy_BETA) <- mat[, 1]
class(methy_BETA) <- "numeric"

# ── Load covariates ──────────────────────────────────────────────────────────
covs_pc_df <- as.data.frame(fread(paste0(wk, cov_file), header = TRUE))
common_samples <- intersect(colnames(methy_BETA), colnames(covs_pc_df))

reordered_methy_BETA <- t(methy_BETA[, common_samples])  # N x CpG
reordered_covs_pc    <- t(covs_pc_df[, common_samples])  # N x covariates
colnames(reordered_covs_pc) <- covs_pc_df[[id_col]]

# ── Set number of PEER factors ───────────────────────────────────────────────
n <- nrow(reordered_methy_BETA)
NumPeerFactorS <- if (n < 50) 5 else 20
print(paste0("Sample size: ", n, ". PEER factors: ", NumPeerFactorS))

# ── Save PEER input files ────────────────────────────────────────────────────
# Without header/rownames (for PEER tool)
fwrite(reordered_methy_BETA,
       paste0(wk, "Methy_CellTypes/", methy_prefix, ".", ct, ".qn.peertool.csv"),
       row.names = FALSE, col.names = FALSE)
fwrite(reordered_covs_pc,
       paste0(wk, cov_file, ".peertool.csv"),
       row.names = FALSE, col.names = FALSE)

# With header/rownames (for reference)
fwrite(reordered_methy_BETA,
       paste0(wk, "Methy_CellTypes/", methy_prefix, ".", ct, ".qn.peertool.names.csv"),
       row.names = TRUE, col.names = TRUE)
fwrite(reordered_covs_pc,
       paste0(wk, cov_file, ".peertool.names.csv"),
       row.names = TRUE, col.names = TRUE)

print(paste0("PEER input files saved for ", dataset, " | ", ct))
