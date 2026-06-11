#conda activate twasguo
library(sqldf)
library("stringr")
library("data.table")
library("dplyr")
args <- commandArgs(trailingOnly=T)

read_strict_fread <- function(file_path, expected_cols) {
  if (!file.exists(file_path)) {
    cat("Warning: File does not exist ->", file_path, "\n")
    return(NULL)}
  lines <- readLines(file_path)
  split_lines <- strsplit(gsub("\\s+", "\t", lines), "\t")
  valid_lines <- split_lines[sapply(split_lines, length) == expected_cols]
  if (length(valid_lines) == 0) {
    cat("Warning: No valid", expected_cols, "column rows in", file_path, "\n")
    return(NULL)}
  texts <- vapply(valid_lines, function(x) paste(x, collapse = "\t"), character(1))
  dt <- fread(input = paste(texts, collapse = "\n"))
  return(dt)
}

###### I. paramter ###### 
CTlist=c("ABS","CT","EE","GOB","TAC","TUF","STM")
CT <- CTlist[as.numeric(args[1])]
model_prefix <- paste0("CRC_EUR_", CT)
wk="/data/l2_bioinfo1/liq17/MethExp/Expr_CRC_EUR/hires54k_crc_eur_genes10k_SingleCellCount_bulkCPM/TFTWAS_RES_add_corr/"


######II. Collect all model summary files [ensure one ENSG name]
model_summary <- matrix(integer(0), nrow = 0, ncol = 7) %>% as.data.frame()
names(model_summary) <- c("gene", "genename" ,"pred.perf.r2", "pred.perf.pval", "pred.perf.qval", "n.snps.in.model", "r.global")

model_weights <- matrix(integer(0), nrow = 0, ncol = 5) %>% as.data.frame()
names(model_weights) <- c("rsid","gene_id","beta","ref","alt") #gene_id rsid    varID   ref     alt     beta

model_covars <- matrix(integer(0), nrow = 0, ncol = 4) %>% as.data.frame()
names(model_covars) <- c('gene_id', 'rsid1', 'rsid2', 'corvarianceValues')

####Collect model summary#####
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_model_summaries.txt"))) {
		model_summary_perchr=as.data.frame(fread(paste0(wk, model_prefix, "_chr",chrom , "_model_summaries.txt")))
		if(dim(model_summary_perchr)[2]==14){
			#model_summary_perchr_remain=model_summary_perchr[, c(1,2,11,7,12,12)] #gene， gene， pred_perf_R2， n_snps_in_model， pred_perf_pval， pred_perf_pval
			model_summary_perchr_remain=model_summary_perchr[, c(1,2,11,12,12,7,14)] #gene， gene， pred_perf_r2， pred_perf_pval, n_snps_in_model， r_global
			model_summary = rbind(model_summary, model_summary_perchr_remain)
		}else{
			cat("Error: reading ",chrom," model summary file does not have 14 columns! Exit")
			q()
		}
	} else {
	  cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_model_summaries.txt does not exit")
	}
}
##Remove pred.perf.R2 is NA and keep only one ENSG
names(model_summary) <- c("gene","genename","pred_perf_r2","pred_perf_pval", "pred_perf_qval", "n_snps_in_model", "r_global")
model_summary <- model_summary[!is.na(model_summary$r_global),]
model_summary_unique <- model_summary %>% arrange(r_global) %>% distinct(gene, .keep_all=TRUE)
model_summary_unique <- model_summary_unique %>% distinct(gene, .keep_all=TRUE)
colnames(model_summary_unique) <- c("gene", "genename" ,"pred.perf.r2", "pred.perf.pval", "pred.perf.qval", "n.snps.in.model", "r.global")
fwrite(model_summary_unique, paste0(wk, model_prefix, "_model_summaries.csv"), row.names = FALSE, col.names = TRUE)

####Collect model weights#####
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_weights.txt"))) {
		model_weights_perchr=read_strict_fread(paste0(wk, model_prefix, "_chr",chrom , "_weights.txt"), 6)
		if(dim(model_weights_perchr)[2]==6){
			model_weights_perchr_remain=model_weights_perchr[, c(2,1,6,4,5)] #cpG     rsid    varID   ref     alt     beta 
			model_weights = rbind(model_weights, model_weights_perchr_remain)
		}else{
			cat("Error: reading ",chrom," model weights file does not have 6 columns! Exit")
			q()
		}
	}else{
		cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_weights.txt does not exit")
	}
}
colnames(model_weights) <- c("rsid","gene","weight","ref_allele","eff_allele")
model_weights <- model_weights[!is.na(model_weights$weight),]
model_weights_unique <- model_weights %>% distinct(gene, rsid , .keep_all=TRUE)

####Collect model covariates#####
#sed -i 's/ /\t/g' Meth_UVA_chr*_covariances.txt
for(i in 1:22){
	chrom=as.character(i)
	if (file.exists(paste0(wk, model_prefix, "_chr",chrom , "_covariances.txt"))) {
		model_covars_perchr=read_strict_fread(paste0(wk, model_prefix, "_chr",chrom , "_covariances.txt"), 4)
		if(dim(model_covars_perchr)[2]==4){
			model_covars = rbind(model_covars, model_covars_perchr)
		}else{
			cat("Error: reading ",chrom," model covariances file does not have 6 columns! Exit")
			q()
		}
	}else{
		cat("Warning: ", model_prefix, "_chr", as.character(chrom), "_covariances.txt does not exit")
	}
}
colnames(model_covars) <-c("GENE","RSID1","RSID2","VALUE")
model_covars <- model_covars[!is.na(model_covars$VALUE),]
model_covars_unique <- model_covars %>% distinct(GENE, RSID1, RSID2, .keep_all=TRUE)
fwrite(model_covars_unique, paste0(wk, model_prefix, "_cov.txt"), sep=" ", row.names = FALSE, col.names = TRUE)
system(paste("gzip", paste0(wk, model_prefix, "_cov.txt")))


######III. Generate a db file for SPrediXcan
Summary <- model_summary_unique
Weight <- model_weights_unique
dbfile <- paste0(wk, model_prefix, ".db")

db <- dbConnect(SQLite(), dbname= dbfile)
dbWriteTable(conn = db, name = "extra", value = Summary, row.names = FALSE, header = TRUE)
dbWriteTable(conn = db, name = "weights", value = Weight,row.names = FALSE, header = TRUE)

