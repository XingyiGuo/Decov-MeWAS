#!/usr/bin/env Rscript

# Usage: Rscript EpiSCORE_TCA.R <project_index> <root_dir>
# project_index: 1=UVA, 2=GTEx
# Steps 1-3 (gene-level methylation aggregation and cell fraction estimation via EpiSCORE)
# are run separately. This script runs Step 4: TCA deconvolution.
#
# Step 4: TCA (Tensor Composition Analysis)
# Estimated runtime: ~180 min; memory: ~20G (UVA), ~200G (GTEx)

library(TCA)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
projects <- c("UVA", "GTEx")
project <- projects[as.numeric(args[1])]
root <- args[2]
print(project)

# Load cell fraction estimates and beta matrix
if (project == "UVA") {
  CellsFraction <- as.data.frame(fread(paste0(root, "/MethExp/Meth_UVA/CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38.GenesMethy.CellFraction_noIC.csv")))
  beta.m <- as.data.frame(fread(paste0(root, "/MethExp/Meth_UVA/CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38.csv")))
} else if (project == "GTEx") {
  CellsFraction <- as.data.frame(fread(paste0(root, "/MethExp/Meth_GTEx/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.GenesMethy.CellFraction_noIC.csv")))
  beta.m <- as.data.frame(fread(paste0(root, "/MethExp/Meth_GTEx/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.csv")))
}

rownames(CellsFraction) <- CellsFraction$V1
CellsFraction$V1 <- NULL
rownames(beta.m) <- beta.m$cpg
beta.m$cpg <- NULL
beta.m$CHR <- NULL
beta.m$MAPINFO_hg38 <- NULL

if (project == "GTEx") {
  rownames(CellsFraction) <- gsub("GTEX\\.", "GTEX-", rownames(CellsFraction))
}

# Load covariates
if (project == "UVA") {
  covs_pc_df <- as.data.frame(fread(paste0(root, "/MethExp/Meth_UVA/146_covariates_Normal_Mucosa.csv"), header = TRUE))
  rownames(covs_pc_df) <- covs_pc_df$id_clx
  covs_pc_df$id_clx <- NULL
} else if (project == "GTEx") {
  covs_pc_df <- as.data.frame(fread(paste0(root, "/MethExp/Meth_GTEx/838_covariates_age_sex_5pc.csv"), header = TRUE))
  rownames(covs_pc_df) <- covs_pc_df$SUBJID
  covs_pc_df$SUBJID <- NULL
}

overlapped_samples_OrderBeta <- colnames(beta.m)[colnames(beta.m) %in% intersect(colnames(beta.m), colnames(covs_pc_df))]
covs_pc_df_refine_T <- t(covs_pc_df[, overlapped_samples_OrderBeta])

if (project == "UVA") {
  c1 <- covs_pc_df_refine_T[, c("sex_num", "age", "Normal")]
  colnames(c1) <- c("gender", "age", "tissue")
} else if (project == "GTEx") {
  c1 <- covs_pc_df_refine_T[, c("SEX", "AGE")]
  colnames(c1) <- c("gender", "age")
}

# Run TCA
tca.col.methy <- tca(X = beta.m[, overlapped_samples_OrderBeta],
                     W = CellsFraction[overlapped_samples_OrderBeta, ],
                     C1 = c1)
output <- tensor(X = as.matrix(beta.m[, overlapped_samples_OrderBeta]), tca.col.methy)

# Save cell-type-specific methylation matrices
cts <- colnames(tca.col.methy$mus_hat)
for (k in seq_along(output)) {
  print(cts[k])
  if (project == "UVA") {
    fwrite(output[[k]],
           paste0(root, "/MethExp/Meth_UVA/Methy_CellTypes/CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38.", cts[k], ".csv"),
           col.names = TRUE, row.names = TRUE)
  } else if (project == "GTEx") {
    fwrite(output[[k]],
           paste0(root, "/MethExp/Meth_GTEx/Methy_CellTypes/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.", cts[k], ".csv"),
           col.names = TRUE, row.names = TRUE)
  }
}
