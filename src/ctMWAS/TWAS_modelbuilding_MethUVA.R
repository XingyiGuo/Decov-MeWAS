library("methods")
library("glmnet")
library("stringr")
library("data.table")
library("dplyr")
args <- commandArgs(trailingOnly=TRUE)

###### I. paramter ###### 
cts<- c("Enteriendocrine", "Enterocyte","Goblet","Progenitor")
ct<-cts[as.numeric(args[1])]

chrom <- as.character(args[2])
print(chrom)

###### II. input files ######
root="/data/l2_bioinfo1/liq17/"
prefix <- paste0(root, "/MethExp/Meth_UVA/Methy_CellTypes/TFTWAS_RES_add_corr/", ct)

snp_annot_file <- paste0(root,"/TF_SNPs_USED/UVA/chr",chrom,"_CRC_TF_used.bed.hg38.in.UK_US.GWAS")

genotype_file <- paste0(root,"/MethExp/Meth_UVA/SNPsImputation2019/plink/chr1_22.dose.filter.recode.R03.rs.",chrom,".hg38.gt.gt")

meth_file <-  paste0(root, "/MethExp/Meth_UVA/Methy_CellTypes/CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38.",ct,".qn.peerresidual.inv.csv")

# vcf head
vcfhead <- read.table(paste0(root, "/MethExp/Meth_UVA/SNPsImputation2019/plink/genotype.header.txt"),header = F, stringsAsFactors = F)
vcfhead <- gsub("(VM\\d+)_\\S+","\\1",vcfhead,perl=T)

###### III. Set function  ######
source(paste0(root,"/GenoPhenoBasic.R"))

do_elastic_net <- function(cis_gt, expr_adj, n_folds, cv_fold_ids, n_times, alpha) {
    cis_gt <- as.matrix(cis_gt)
    fit <- cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
    lambda_seq <- fit$lambda
    cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
    fits <- list()
    fits[[1]] <- fit
    cvms <- matrix(nrow = 100, ncol = n_times)
    cvms[1:length(fit$cvm),1] <- fit$cvm
    for (i in 2:(n_times)) {
      fit <- cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
      fits[[i]] <- fit
      cvms[1:length(fit$cvm),i] <- fit$cvm
    }
    avg_cvm <- rowMeans(cvms)
    best_lam_ind <- which.min(avg_cvm)
    best_lambda <- lambda_seq[best_lam_ind]
    out <- list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda)
    out
}

evaluate_performance <- function(cis_gt, expr_adj, fit, best_lam_ind, best_lambda, cv_fold_ids, n_folds) {
  n_nonzero <- fit$nzero[best_lam_ind]
  if (n_nonzero > 0) {
    R2 <- rep(0, n_folds)
    for (j in (1:n_folds)) {
      fold_idxs <- which(cv_fold_ids[,1] == j)
      tss <- sum(expr_adj[fold_idxs]**2)
      rss <- sum((expr_adj[fold_idxs] - fit$fit.preval[fold_idxs, best_lam_ind])**2)
      R2[j] <- 1 - (rss/tss)
    }
    best_fit <- fit$glmnet.fit
    expr_adj_pred <- predict(best_fit, as.matrix(cis_gt), s = best_lambda)
    tss_best_fit <- sum(expr_adj**2)
    rss_best_fit <- sum((expr_adj - expr_adj_pred)**2)
	R2_mean <- mean(R2)
    R2_sd <- sd(R2)
    inR2 <- 1 - (rss_best_fit/tss_best_fit) # All sample, in sample R2
    
    n_samp <- length(expr_adj)
    weights <- best_fit$beta[which(best_fit$beta[,best_lam_ind] != 0), best_lam_ind]
    weighted_snps <- names(best_fit$beta[,best_lam_ind])[which(best_fit$beta[,best_lam_ind] != 0)]

    # old way (keep to trace errors)
    pred_perf <- summary(lm(expr_adj ~ fit$fit.preval[,best_lam_ind]))
    pred_perf_rsq <- pred_perf$r.squared
    pred_perf_pval <- pred_perf$coef[2,4]
	
	# up to date way to calcuate coefficient of determination
	y    <- expr_adj
	yhat <- fit$fit.preval[, best_lam_ind]
	tss_global <- sum((y - mean(y))^2)
	rss_global <- sum((y - yhat)^2)
	R2_global  <- 1 - rss_global/tss_global
	r_global <- cor(y, yhat)
	
    out <- list(weights = weights, n_weights = n_nonzero, weighted_snps = weighted_snps, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pred_perf_rsq = pred_perf_rsq, pred_perf_pval = pred_perf_pval, R2_global=R2_global, r_global=r_global)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, weighted_snps = NA, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pred_perf_rsq = NA, pred_perf_pval = NA, R2_global=NA, r_global=NA)
  }
  out
}

###### IV. run analysis ######
# parameter for analysis 
maf=0.05
n_times=3
n_k_folds=10
cis_window=1000000
alpha=0.5

# 1. Inite analysis
snp_annot=fread(snp_annot_file)

# 2. Expression cpG check
meth_df_all=as.data.frame(fread(meth_file))
meth_df_perchr=filter(meth_df_all, CHR==chrom)
meth_df_perchr_cpGs=as.vector(meth_df_perchr$cpg)
meth_df_perchr_cpGs_pos =as.vector(meth_df_perchr$MAPINFO_hg38)

meth_df_perchr=meth_df_perchr[,-(1:3)]
meth_df_perchr_t=t(meth_df_perchr)
colnames(meth_df_perchr_t) = meth_df_perchr_cpGs
rownames(meth_df_perchr_t) = as.vector(colnames(meth_df_all))[4:length(colnames(meth_df_all))]

if(!is.null(meth_df_perchr_t) & (dim(meth_df_perchr_t)[2]!=0)){
  print(paste0("Number of meth cpG to be analyzed ", as.character(dim(meth_df_perchr_t)[2])))
}else{
  print("No cpG remain, exit")
  q()
}
cpGs <- as.vector(colnames(meth_df_perchr_t))
n_cpGs <- length(cpGs)

# 3. Organize gentoyep and expression to contain only common subjects
samples<-intersect(vcfhead[6:length(vcfhead)], rownames(meth_df_perchr_t))
n_samples <- length(samples)
print(paste0("Samples for both genotype and phenotype ", as.character(length(samples))))
if(length(samples)<100){
  print("Errors in matching samples in genotype and expression. Please check your gentoype files and expression files!Exit")
  q()}

gt_df <- get_maf_filtered_genotype(genotype_file, vcfhead, maf, rownames(meth_df_perchr_t))
#Order gene expression rows by gt sampels
indices <- match(row.names(gt_df), rownames(meth_df_perchr_t))
meth_df_perchr_t_match_gt <- meth_df_perchr_t[indices, , drop = FALSE]
if(dim(gt_df)[1] != dim(meth_df_perchr_t_match_gt)[1]){
  print("Selected genotype and expression do not have the same number of samples! Exit")
  q()}

# 4. Output model res
seed <- 20240806
set.seed(seed)
cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
for (j in 1:n_times){
cv_fold_ids[,j] <- sample(1:n_k_folds, n_samples, replace = TRUE)
}

# output file
model_summary_file <- paste0(prefix,'_chr',chrom,'_model_summaries.txt')
model_summary_cols <- c('cpG', 'cpG', 'alpha', 'cv_mse', 'lambda_iteration', 'lambda_min', 'n_snps_in_model',
					  'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2', 'pred_perf_R2', 'pred_perf_pval','R2_global','r_global')
write(model_summary_cols, file = model_summary_file, ncol = 14, sep = '\t')

weights_file <- paste0(prefix,'_chr',chrom,'_weights.txt')
weights_col <- c('cpG', 'rsid', 'varID', 'ref', 'alt', 'beta')
write(weights_col, file = weights_file, ncol = 6, sep = '\t')

covariance_file <- paste0(prefix,'_chr',chrom,'_covariances.txt')
covariance_col <- c('cpG', 'rsid1', 'rsid2', 'corvarianceValues')
write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')

for (i in 1:n_cpGs) { #
	cat(i, "/", n_cpGs, "\n")
	cpG <- cpGs[i]
	model_summary <- c(cpG, cpG, alpha, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA, NA)
	coords <- c(meth_df_perchr_cpGs_pos[i], meth_df_perchr_cpGs_pos[i])
	cis_gt <- get_cis_genotype_varID(gt_df, snp_annot, "pos", coords, cis_window)
	print(dim(cis_gt))
	tryCatch({    
		if (ncol(cis_gt) >= 2) {
			adj_expression <- meth_df_perchr_t_match_gt[,i]
			elnet_out <- do_elastic_net(cis_gt, adj_expression, n_k_folds, cv_fold_ids, n_times, alpha)
			if (length(elnet_out) > 0) {
				eval <- evaluate_performance(cis_gt, adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind, elnet_out$best_lambda, cv_fold_ids, n_k_folds)
				model_summary <- c(cpG, cpG, alpha, elnet_out$min_avg_cvm, elnet_out$best_lam_ind,
								   elnet_out$best_lambda, eval$n_weights, eval$R2_mean, eval$R2_sd, eval$inR2,
								   eval$pred_perf_rsq, eval$pred_perf_pval, eval$R2_global, eval$r_global)
				if (eval$n_weights > 0) {
				  weighted_snps_info <- snp_annot %>% filter(varID %in% eval$weighted_snps) %>% select(SNP, varID, ref , effect)
				  if (nrow(weighted_snps_info) == 0)
					browser()
				  weighted_snps_info$gene <- cpG
				  weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights = eval$weights, varID=eval$weighted_snps), by = 'varID') %>% select(gene, SNP, varID, ref, effect, weights)
				  write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
				  do_covariance(cpG, cis_gt, weighted_snps_info$SNP, weighted_snps_info$varID, covariance_file)
				}
			}
			write(model_summary, file = model_summary_file, append = TRUE, ncol = 14, sep = '\t')
			print(model_summary)
		}else{
			print("Less than 2 genetic variants, not able to conduct ent")
		}
	},error=function(e){})
}