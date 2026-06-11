#!/usr/bin/env Rscript

# Usage: Rscript CRC_PEER.R <dataset> <index1> <index2> <working_dir>
#
# <dataset>     : "EUR" or "GTEx"
# <index1>      : cell type index (EUR), or tissue index (GTEx: 1=Colon, 2=Breast, 3=Lung, 4=Prostate)
# <index2>      : cell type index within tissue (GTEx only; omit or set to NA for EUR)
# <working_dir> : path to working directory containing expression matrices and covariate files
#
# PEER factor selection follows GTEx guidelines:
#   N < 150 -> 15, 150-249 -> 30, 250-349 -> 45, >= 350 -> 60

library(peer)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]

# ── Cell type lists ──────────────────────────────────────────────────────────
EUR_CTlist <- c("ABS", "CT", "EE", "GOB", "TAC", "TUF", "STM")

GTEx_CTlist <- list(
  c("ABS", "CT", "EE", "GOB", "STM", "TAC", "TUF"),          # 1: Colon transverse
  c("cluster1", "cluster3", "cluster12", "cluster17", "cluster35", "cluster39"),  # 2: Breast
  c("cluster3", "cluster5", "cluster6", "cluster8", "cluster9", "cluster10", "cluster17", "cluster39"),  # 3: Lung
  c("cluster3", "cluster4", "cluster8", "cluster10", "cluster12", "cluster17", "cluster34", "cluster39") # 4: Prostate
)
GTEx_cov_prefix <- c("colontrans_368", "breast_151", "lung_515", "prostate_221")

# ── Parse args and set paths ─────────────────────────────────────────────────
if (dataset == "EUR") {
  CT  <- EUR_CTlist[as.numeric(args[2])]
  wk  <- args[3]
} else if (dataset == "GTEx") {
  ts_index <- as.numeric(args[2])
  CT  <- GTEx_CTlist[[ts_index]][as.numeric(args[3])]
  wk  <- args[4]
  if (ts_index > 4) stop("GTEx tissue index must be 1-4 (Colon/Breast/Lung/Prostate)")
} else {
  stop("dataset must be 'EUR' or 'GTEx'")
}
print(paste0("Dataset: ", dataset, " | Cell type: ", CT))

# ── Load expression matrix ───────────────────────────────────────────────────
mat <- as.data.frame(fread(paste0(wk, CT, ".10k.tpm.qn.invern.csv")))
col_num <- dim(mat)[2]
mat_gene_expression <- as.matrix(mat[, 4:col_num])
class(mat_gene_expression) <- "numeric"

# ── Load covariates ──────────────────────────────────────────────────────────
if (dataset == "EUR") {
  covs_raw_df <- as.data.frame(fread(paste0(wk, "AgeSexCovariate.tsv"), header = TRUE))
  covs_age_sex <- covs_raw_df[, 2:4] %>%
    mutate(sex = recode(gender, "Male" = 0, "Female" = 1))
  pc <- as.data.frame(fread(paste0(wk, "samples_pruned_merge.pca.eigenvec"), header = FALSE))
  pc_5 <- pc[, c("V1", "V3", "V4", "V5", "V6", "V7")]
  names(pc_5) <- c("VMID", "pc1", "pc2", "pc3", "pc4", "pc5")
  covs_merged <- inner_join(covs_age_sex, pc_5, by = "VMID")
  covs_common <- t(covs_merged[covs_merged$VMID %in% colnames(mat_gene_expression), c(2, 4:9)])
  colnames(covs_common) <- colnames(mat_gene_expression)
  class(covs_common) <- "numeric"
  if (ncol(covs_common) != ncol(mat_gene_expression)) stop("Sample mismatch between covariates and expression matrix")

} else {  # GTEx
  cov_raw <- as.data.frame(fread(paste0(wk, GTEx_cov_prefix[ts_index], "_covariats.tsv"), header = TRUE))
  covs_common <- t(cov_raw)
  if (ncol(covs_common) != ncol(mat_gene_expression)) {
    mat_gene_expression <- mat_gene_expression[, colnames(mat_gene_expression) %in% cov_raw$ID]
    print(paste0("Subsetted to ", ncol(mat_gene_expression), " overlapping samples"))
  }
}
print(paste0("Samples: ", ncol(mat_gene_expression)))

# ── Set number of PEER factors ───────────────────────────────────────────────
n <- ncol(mat_gene_expression)
NumPeerFactorS <- if (n < 150) 15 else if (n < 250) 30 else if (n < 350) 45 else 60
print(paste0("PEER factors: ", NumPeerFactorS))

# ── Run PEER ─────────────────────────────────────────────────────────────────
model <- PEER()
PEER_setPhenoMean(model, t(as.matrix(mat_gene_expression)))  # n x e
PEER_setNk(model, NumPeerFactorS)
PEER_getNk(model)
PEER_setCovariates(model, t(as.matrix(covs_common)))
PEER_update(model)

# ── Diagnostic plot ──────────────────────────────────────────────────────────
pdf(paste0(wk, CT, ".peer.diag.pdf"), width = 6, height = 8)
PEER_plotModel(model)
dev.off()

# ── Save residuals ───────────────────────────────────────────────────────────
residuals <- t(PEER_getResiduals(model))  # e x n
colnames(residuals) <- colnames(mat_gene_expression)
residuals_final <- cbind(mat[, 1:3], residuals)
fwrite(residuals_final,
       file = paste0(wk, CT, ".10k.tpm.qn.invern.afterpeer.csv"),
       row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
print(paste0("PEER correction finished for ", CT))
