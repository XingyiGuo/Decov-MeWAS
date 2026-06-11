library("methods")
library("glmnet")
library("stringr")
library("data.table")
library("dplyr")
args <- commandArgs(trailingOnly=TRUE)

###### I. paramter ###### 
CTlist=c("ABS","CT","EE","GOB","STM","TAC","TUF")
CT <- CTlist[as.numeric(args[1])]
chrom <- as.character(args[2])
print(CT)
print(chrom)

"%&%" <- function(a,b) paste(a,b, sep = "")
prefix <- paste0("CRC_EUR_",CT)

###### II. input files ######
snp_annot_file <- paste0("/data/l2_bioinfo1/liq17/TF_SNPs_USED/UVA/chr",chrom,"_CRC_TF_used.bed.hg38.in.UK_US.GWAS")

genotype_file <- paste0("/data/l2_bioinfo1/liq17/MethExp/Genotype_CRC_EUR/plink/chr", chrom,".gt.gt")

gene_annot_file  <- "/data/l2_bioinfo1/liq17/ref/gencode/gencode.v26.annotation.gtf"

expression_file <-  paste0("/data/l2_bioinfo1/liq17/MethExp/Expr_CRC_EUR/hires54k_crc_eur_genes10k_SingleCellCount_bulkCPM/",CT,".10k.tpm.qn.invern.afterpeer.invern.csv") # qn -> inverse -> peer -> inverse

# vcf head
vcfhead <- read.table("/data/l2_bioinfo1/liq17/MethExp/Genotype_CRC_EUR/vcf.head",header = F, stringsAsFactors = F)
vcfhead <- gsub("(VM\\d+)_\\S+","\\1",vcfhead,perl=T)

###### III. Set function  ######
source("/data/l2_bioinfo1/liq17/GenoPhenoBasic.R")

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
	
	r_squared_global <- cor(y, yhat)
	
    out <- list(weights = weights, n_weights = n_nonzero, weighted_snps = weighted_snps, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pred_perf_rsq = pred_perf_rsq, pred_perf_pval = pred_perf_pval, R2_global=R2_global, r_squared_global=r_squared_global)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, weighted_snps = NA, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pred_perf_rsq = NA, pred_perf_pval = NA, R2_global=NA, r_squared_global=NA)
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
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
snp_annot <- fread(snp_annot_file)

# 2. Expression genes check
expr_df <- get_gene_expression(expression_file, gene_annot)
if(!is.null(expr_df) & (dim(expr_df)[2]!=0)){
  print(paste0("Number of genes to be analyzed ", as.character(dim(expr_df)[2])))
}else{
  print("No genes remain, exit")
  q()
}
genes <- colnames(expr_df)
n_genes <- length(expr_df)

# 3. Organize gentoyep and expression to contain only common subjects
samples<-intersect(vcfhead[6:length(vcfhead)], rownames(expr_df))
n_samples <- length(samples)
print(paste0("Samples for both genotype and phenotype ", as.character(length(samples))))
if(length(samples)<100){
  print("Errors in matching samples in genotype and expression. Please check your gentoype files and expression files!Exit")
  q()}

gt_df <- get_maf_filtered_genotype_ID(genotype_file, vcfhead, maf, rownames(expr_df))
#Order gene expression rows by gt sampels
indices <- match(row.names(gt_df), rownames(expr_df))
expr_df_match_gt <- expr_df[indices, , drop = FALSE]
if(dim(gt_df)[1] != dim(expr_df_match_gt)[1]){
  print("Selected genotype and expression do not have the same number of samples! Exit")
  q()}


# 4. Output model res
seed <- 20240528
set.seed(seed)
cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
for (j in 1:n_times)
cv_fold_ids[,j] <- sample(1:n_k_folds, n_samples, replace = TRUE)

# output file
model_summary_file <- paste0(prefix,'_chr',chrom,'_model_summaries.txt')
model_summary_cols <- c('gene', 'gene_name', 'alpha', 'cv_mse', 'lambda_iteration', 'lambda_min', 'n_snps_in_model',
					  'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2', 'pred_perf_R2', 'pred_perf_pval','R2_global','r_squared_global')
write(model_summary_cols, file = model_summary_file, ncol = 14, sep = '\t')

weights_file <- paste0(prefix,'_chr',chrom,'_weights.txt')
weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
write(weights_col, file = weights_file, ncol = 6, sep = '\t')

covariance_file <- paste0(prefix,'_chr',chrom,'_covariances.txt')
covariance_col <- c('gene_id', 'rsid1', 'rsid2', 'corvarianceValues')
write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')

for (i in 1:n_genes) {
	cat(i, "/", n_genes, "\n")
	gene <- genes[i]
	gene_name <- as.character(gene_annot$genename[gene_annot$geneid == gene])
	model_summary <- c(gene, gene_name, alpha, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA, NA)
	coords <- get_gene_coords(gene_annot, gene)
	cis_gt <- get_cis_genotype_varID(gt_df, snp_annot, "pos", coords, cis_window)

	tryCatch({    
		if (ncol(cis_gt) >= 2) {
			adj_expression <- expr_df_match_gt[,i]
			elnet_out <- do_elastic_net(cis_gt, adj_expression, n_k_folds, cv_fold_ids, n_times, alpha)
			if (length(elnet_out) > 0) {
				eval <- evaluate_performance(cis_gt, adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind, elnet_out$best_lambda, cv_fold_ids, n_k_folds)
				model_summary <- c(gene, as.character(gene_name), alpha, elnet_out$min_avg_cvm, elnet_out$best_lam_ind,
								   elnet_out$best_lambda, eval$n_weights, eval$R2_mean, eval$R2_sd, eval$inR2,
								   eval$pred_perf_rsq, eval$pred_perf_pval, eval$R2_global, eval$r_squared_global)
				if (eval$n_weights > 0) {
				  weighted_snps_info <- snp_annot %>% filter(varID %in% eval$weighted_snps) %>% select(SNP,varID, ref , effect)
				  if (nrow(weighted_snps_info) == 0)
					browser()
				  weighted_snps_info$gene <- gene
				  weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights = eval$weights, varID=eval$weighted_snps), by = 'varID') %>% select(gene, SNP, varID, ref, effect, weights)
				  write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
				  do_covariance(gene, cis_gt, weighted_snps_info$SNP, weighted_snps_info$varID, covariance_file)
				}
			}
		}
		print(model_summary)
		write(model_summary, file = model_summary_file, append = TRUE, ncol = 14, sep = '\t')
	},error=function(e){})
}
