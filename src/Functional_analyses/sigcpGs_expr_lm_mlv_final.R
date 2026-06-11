# conda activate tftwas_coloc
library("methods")
library("glmnet")
library("stringr")
library("data.table")
library("dplyr")
library("readxl")
setwd("/data/sbcs/GuoLab/backup/liq17/MethExp/Expr_Methy_sigCPGs_lm/202603")
args <- commandArgs(trailingOnly = TRUE)

set.seed(202603)

###### I. Parameters ######
eur_source <- "GTEx"
prefix     <- paste0("Meth_sigCPGs_lm_", eur_source)

# This script has two modes:
#   Rscript sigcpGs_expr_lm_mlv.R <chrom>         -> per-chromosome LM  (chrom 1-22)
#   Rscript sigcpGs_expr_lm_mlv.R aggregate_mlv   -> aggregate + MLV post-processing
run_mode <- ifelse(length(args) > 0, args[1], stop("Provide chromosome number or 'aggregate_mlv'"))

###### II. Input files ######
sigCPGs_df <- as.data.frame(read.csv(
  "/data/sbcs/GuoLab/backup/liq17/MethExp/Manuscript_Tables_Figures/Methy_Cauchy/CRC_EUR_GTEX_Methylation_META_4cts_updateGenes_reportCPGs_GWASLoci.csv"))

CPGs_beta <- as.data.frame(fread(
  "/data/sbcs/GuoLab/backup/liq17/MethExp/Meth_GTEx/Methy_bulk/GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.qn.afterpeer.ivn.hg38.csv"))
sigCPGs_beta <- CPGs_beta[CPGs_beta$cpG %in% sigCPGs_df$cpg, ]
rownames(sigCPGs_beta) <- NULL

# Gene expression: samples as rows, genes as columns (E x N transposed)
expr_df <- as.data.frame(fread(
  "/data/sbcs/GuoLab/backup/RNA-seq/GTEx_CRC_RNA-seq/processed/expr_as_apa/GTEx_Colon_Transverse_expression.txt_QN_inverse_PEER_inverse_02092023.csv"))
rownames(expr_df)    <- expr_df$V1
expr_df$V1           <- NULL
expr_df_ExN          <- t(expr_df)
rownames(expr_df_ExN) <- sub("\\.\\d+$", "", rownames(expr_df_ExN))

# Align samples
common_samples         <- intersect(colnames(sigCPGs_beta), colnames(expr_df_ExN))
print(paste0(length(common_samples), " common Eur samples in both Expr and Methy"))
common_samples_ordered <- colnames(sigCPGs_beta)[colnames(sigCPGs_beta) %in% common_samples]

sigCPGs_beta_ordered <- sigCPGs_beta[, common_samples_ordered, drop = FALSE]
sigCPGs_beta_ordered <- cbind(sigCPGs_beta[, 1:3], sigCPGs_beta_ordered)
expr_df_ExN_ordered  <- expr_df_ExN[, common_samples_ordered, drop = FALSE]

gene_annot_file <- "/data/sbcs/GuoLab/backup/liq17/ref/gencode/gencode.v26.annotation.gtf"

###### III. Helper functions ######
get_gene_annotation <- function(gene_annot_file, chrom, bp, flank) {
  gene_df  <- read.table(gene_annot_file, header = FALSE, stringsAsFactors = FALSE,
                         sep = "\t", fill = TRUE)
  gene_df1 <- filter(gene_df, V3 %in% "gene")
  geneid   <- str_extract(gene_df1[, 9], "ENSG\\d+.\\d+")
  geneid   <- sub("\\.\\d+$", "", geneid)
  genename <- gsub("gene_name (\\S+);", "\\1",
                   str_extract(gene_df1[, 9], "gene_name (\\S+);"), perl = TRUE)
  gene_used <- as.data.frame(cbind(geneid, genename, gene_df1[, c(1, 4, 5, 3)]))
  colnames(gene_used) <- c("geneid", "genename", "chr", "start", "end", "anno")
  gtf_used <- filter(gene_used, gene_used[, 3] %in% paste0("chr", chrom))
  variant_flank_start <- max(as.numeric(bp) - flank, 0)
  variant_flank_end   <- as.numeric(bp) + flank
  gtf_used2 <- filter(gtf_used,
                      ((gtf_used[, 4] > variant_flank_start) & (gtf_used[, 4] < variant_flank_end)) |
                      ((gtf_used[, 5] > variant_flank_start) & (gtf_used[, 5] < variant_flank_end)))
  gtf_used2
}

format_output <- function(model_coefficients) {
  int_beta   <- model_coefficients["(Intercept)", "Estimate"]
  int_se     <- model_coefficients["(Intercept)", "Std. Error"]
  int_t      <- model_coefficients["(Intercept)", "t value"]
  int_p      <- model_coefficients["(Intercept)", "Pr(>|t|)"]
  slope_beta <- model_coefficients["as.numeric(as.vector(cpG_betas))", "Estimate"]
  slope_se   <- model_coefficients["as.numeric(as.vector(cpG_betas))", "Std. Error"]
  slope_t    <- model_coefficients["as.numeric(as.vector(cpG_betas))", "t value"]
  slope_p    <- model_coefficients["as.numeric(as.vector(cpG_betas))", "Pr(>|t|)"]
  return(c(slope_beta, slope_se, slope_t, slope_p, int_beta, int_se, int_t, int_p))
}

###### IV. Mode A – Per-chromosome simple LM ######
if (run_mode != "aggregate_mlv") {

	chrom <- as.numeric(run_mode)
	maf   <- 0.01
	flank <- 1000000

	model_summary_file <- paste0(prefix, "_", chrom, "_model_summaries.txt")
	model_summary_cols <- c("cpG", "chrom", "bp_hg38", "gene_name",
						  "slope_beta", "slope_se", "slope_t", "slope_p",
						  "int_beta", "int_se", "int_t", "int_p")
	write(model_summary_cols, file = model_summary_file, ncol = 13, sep = "\t")

	one_chrom_sigCPGs_df <- filter(sigCPGs_beta_ordered, sigCPGs_beta_ordered[, "CHR"] == chrom)

	for (row_index in 1:nrow(one_chrom_sigCPGs_df)) {
	if (dim(one_chrom_sigCPGs_df)[1] > 0) {
	  bp        <- one_chrom_sigCPGs_df[row_index, "MAPINFO_b38"]
	  cpG       <- one_chrom_sigCPGs_df[row_index, "cpG"]
	  cat(chrom, " ", bp, " ", cpG)
	  cpG_betas <- one_chrom_sigCPGs_df[row_index, 4:dim(one_chrom_sigCPGs_df)[2]]

	  one_chrom_gene_annot <- get_gene_annotation(gene_annot_file, chrom, bp, flank)
	  cat("\nNumber of annot (ENSG/ENST) is ", dim(one_chrom_gene_annot)[1], "\n")

	  overlapping_genes          <- rownames(expr_df_ExN_ordered)[rownames(expr_df_ExN_ordered) %in% one_chrom_gene_annot$geneid]
	  expr_df_ExN_ordered_one_chr <- expr_df_ExN_ordered[overlapping_genes, ]
	  n_tests <- dim(expr_df_ExN_ordered_one_chr)[1]
	  tests   <- rownames(expr_df_ExN_ordered_one_chr)
	  print(paste0("#Num of test is ", n_tests))
	  print(paste0("#Num of samples is ", dim(expr_df_ExN_ordered_one_chr)[2]))

	  if (n_tests != 0) {
		for (gene_index in 1:n_tests) {
		  cat(gene_index, "/", n_tests, "\n")
		  testID                    <- tests[gene_index]
		  profile_df_one_gene_testID <- expr_df_ExN_ordered_one_chr[testID, ]
		  if (length(profile_df_one_gene_testID) > 0) {
			model <- lm(as.numeric(as.vector(profile_df_one_gene_testID)) ~
						  as.numeric(as.vector(cpG_betas)))
		  }
		  model_summary <- c(cpG, chrom, bp, testID, format_output(summary(model)$coefficients))
		  write(model_summary, file = model_summary_file, append = TRUE, ncol = 13, sep = "\t")
		}
	  } else {
		cat("\nNo genes in flank of tested cpGs")
	  }
	} else {
	  cat("\nNo sig cpGs from chromosome", chrom)
	}
	}

	cat("\nPer-chromosome LM complete for chromosome", chrom, "\n")

###### V. Mode B – Aggregate LM results + Forward-selection MLV ######
} else {
	cat("\n=== Aggregating per-chromosome LM results ===\n")
	# ---- V.1  Parse full GTF for geneid <-> genename lookup ----
	gtf_raw   <- read.table(gene_annot_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
	gtf_genes <- filter(gtf_raw, V3 == "gene")
	gtf_annot <- data.frame(
	geneid   = sub("\\.\\d+$", "", str_extract(gtf_genes[, 9], "ENSG\\d+\\.\\d+")),
	genename = gsub("gene_name (\\S+);", "\\1", str_extract(gtf_genes[, 9], "gene_name (\\S+);"), perl = TRUE),stringsAsFactors = FALSE)
	gtf_annot <- gtf_annot[!duplicated(gtf_annot$geneid), ]
	# ---- V.2  Concatenate per-chromosome LM files ----
	cpgs_reg <- as.data.frame(rbindlist(lapply(1:22, function(chrom) {
		f <- paste0(prefix, "_", chrom, "_model_summaries.txt")
		if (file.exists(f)) fread(f, sep = "\t") else NULL})))

	# Rename columns to downstream convention
	colnames(cpgs_reg)[colnames(cpgs_reg) == "cpG"]        <- "cpg"
	colnames(cpgs_reg)[colnames(cpgs_reg) == "gene_name"]  <- "geneid"
	colnames(cpgs_reg)[colnames(cpgs_reg) == "slope_beta"] <- "lm_beta"
	colnames(cpgs_reg)[colnames(cpgs_reg) == "slope_se"]   <- "lm_se"
	colnames(cpgs_reg)[colnames(cpgs_reg) == "slope_t"]    <- "lm_t"
	colnames(cpgs_reg)[colnames(cpgs_reg) == "slope_p"]    <- "lm_p"

	cat(paste0("Total CpG-gene pairs from LM: ", nrow(cpgs_reg), "\n")) #120948

	# ---- V.3  Build lead-CpG table (one row per CpG, strongest association) ----
	# Sort by lm_t ascending (most negative = strongest negative effect)
	cpgs_reg_sort    <- cpgs_reg[order(cpgs_reg$lm_t, cpgs_reg$lm_p), ]
	cpgs_reg_sort <- merge(cpgs_reg_sort,
							gtf_annot[, c("geneid", "genename")],
							by = "geneid", all.x = TRUE)
	# Only pairs with beta < 0 and p < 0.05 from univariate LM are candidates
	# cpgs_reg_sort$fdr <- p.adjust(cpgs_reg_sort$lm_p, method = "BH")
	cpgs_reg_sig <- cpgs_reg_sort[
	!is.na(cpgs_reg_sort$lm_p)    & cpgs_reg_sort$lm_p    < 0.01 &
	!is.na(cpgs_reg_sort$lm_beta) & cpgs_reg_sort$lm_beta < 0,
	]
	cat(paste0("Significant CpG-gene pairs (beta<0, p<0.05): ", nrow(cpgs_reg_sig), "\n"))

	# ---- V.4  Prepare matrices for MLV ----
	cpg_betas_full <- CPGs_beta
	rownames(cpg_betas_full) <- cpg_betas_full$cpG
	cpg_betas_full$cpG       <- NULL
	sample_cols              <- intersect(colnames(cpg_betas_full), common_samples_ordered)
	cpg_betas_match_T        <- as.data.frame(t(cpg_betas_full[, sample_cols]))  # samples x CpGs
	expr_samples_x_genes     <- as.data.frame(t(expr_df_ExN_ordered))            # samples x genes

	# ---- V.5  Identify genes with >1 significant CpG from Step 1 ----
	cpg_counts_per_gene <- cpgs_reg_sig %>% group_by(geneid) %>% summarise(n_cpg = n_distinct(cpg), .groups = "drop")
	multi_cpg_genes <- cpg_counts_per_gene$geneid[cpg_counts_per_gene$n_cpg > 1]
	single_cpg_genes <- cpg_counts_per_gene$geneid[cpg_counts_per_gene$n_cpg == 1]
	cat(paste0("Genes with >1 significant CpG: ", length(multi_cpg_genes), "\n"))
	cat(paste0("Genes with 1 significant CpG: ", length(single_cpg_genes), "\n"))

	# ---- V.6  Forward-selection MLV loop ----
	cat("=== Running forward-selection MLV ===\n")
	Allensg_mlv <- data.frame()
	for (ensg_index in 1:length(multi_cpg_genes)) {
		if (ensg_index %% 100 == 0)
		  cat(paste0(ensg_index, " / ", length(multi_cpg_genes), "\n"))
		ensg <- multi_cpg_genes[ensg_index]
		# Candidate CpGs: only significant pairs from Step 1, ordered by lm_t (strongest first)
		gene_sig_pairs <- cpgs_reg_sig[cpgs_reg_sig$geneid == ensg, ]
		gene_sig_pairs <- gene_sig_pairs[order(gene_sig_pairs$lm_t), ]
		cpgs_avail     <- unique(gene_sig_pairs$cpg)
		if (!(ensg %in% colnames(expr_samples_x_genes))) next
		cpgs_avail <- cpgs_avail[cpgs_avail %in% colnames(cpg_betas_match_T)]
		if (length(cpgs_avail) < 2) {
			cat(paste0("  [info] ", ensg, ": 1 CpG\n"))
			next
		}
		expr_vec <- as.numeric(expr_samples_x_genes[, ensg])
		beta_mat <- as.matrix(cpg_betas_match_T[, cpgs_avail])
		# ── Forward selection: add CpGs one at a time while conditional p < 0.05 ──
		selected_cpgs  <- c()
		remaining_cpgs <- cpgs_avail
		p_entry        <- 0.01
		repeat {
			if (length(remaining_cpgs) == 0) break
			# Test each remaining CpG added conditionally on already-selected CpGs
			candidate_p <- sapply(remaining_cpgs, function(cpg_candidate) {
				x_current <- if (length(selected_cpgs) == 0) {
				  matrix(beta_mat[, cpg_candidate], ncol = 1,
						 dimnames = list(NULL, cpg_candidate))
				} else {
				  cbind(beta_mat[, selected_cpgs, drop = FALSE],
						beta_mat[, cpg_candidate])
				}
				fit <- tryCatch(lm(expr_vec ~ x_current), error = function(e) NULL)
				if (is.null(fit)) return(1.0)
				coef_s <- summary(fit)$coefficients
				coef_s[nrow(coef_s), "Pr(>|t|)"]
			})
			best_cpg <- remaining_cpgs[which.min(candidate_p)]
			best_p   <- min(candidate_p)
			if (best_p >= p_entry) break
			selected_cpgs  <- c(selected_cpgs, best_cpg)
			remaining_cpgs <- remaining_cpgs[remaining_cpgs != best_cpg]
		}
		# ── Determine method tag based on how many CpGs were selected ────────────
		if (length(selected_cpgs) >= 2) {
		  mlv_method <- "ForwardOLS_multi"
		  cat(paste0("  [info] ", ensg, ": ",length(cpgs_avail)," CpGs\n"))
		} else if (length(selected_cpgs) == 1) {
		  # Only lead CpG has independent effect; others are redundant
		  mlv_method <- "ForwardOLS_single"
		  cat(paste0("  [info] ", ensg, ": 1 independent CpG from ", length(cpgs_avail), " candidates\n"))
		} else {
		  # No CpG passed entry — edge case guard
		  cat(paste0("  [warn] ", ensg, ": no CpG passed forward selection\n"))
		  next
		}
		# ── Fit final joint model on selected CpGs only ───────────────────────────
		final_mat <- beta_mat[, selected_cpgs, drop = FALSE]
		colnames(final_mat) <- selected_cpgs
		model_final <- tryCatch(lm(expr_vec ~ final_mat), error = function(e) NULL)
		if (is.null(model_final)) next
		coef_summary   <- summary(model_final)$coefficients
		coef_rows      <- rownames(coef_summary)[rownames(coef_summary) != "(Intercept)"]
		coef_cpg_names <- if (length(coef_rows) > 1) sub("^final_mat", "", coef_rows) else selected_cpgs
		fitted_cpgs    <- intersect(selected_cpgs, coef_cpg_names)
		if (length(fitted_cpgs) == 0) next
		matched_rows   <- if (length(fitted_cpgs) > 1) paste0("final_mat", fitted_cpgs) else "final_mat"
		mlv_betas    <- as.numeric(coef_summary[matched_rows, "Estimate"])
		mlv_pvals    <- as.numeric(coef_summary[matched_rows, "Pr(>|t|)"])
		ensg_mlv <- data.frame(
		  geneid            = ensg,
		  cpg               = fitted_cpgs,
		  mlv_beta          = mlv_betas,
		  mlv_p             = mlv_pvals,
		  mlv_method        = mlv_method,
		  cpgs_from_lm       = length(cpgs_avail),   # sig CpGs from Step 1
		  cpgs_selected_mlv     = length(fitted_cpgs),  # independently significant CpGs
		  stringsAsFactors  = FALSE,
		  row.names         = NULL
		)
		# Rank CpGs within this gene by conditional p-value (rank 1 = lead regulatory CpG)
		ensg_mlv          <- ensg_mlv[order(ensg_mlv$mlv_p), ]
		ensg_mlv$cpg_rank <- seq_len(nrow(ensg_mlv))
		Allensg_mlv <- rbind(Allensg_mlv, ensg_mlv)
	}
	cat(paste0("MLV complete. Rows: ", nrow(Allensg_mlv), "\n"))
	cat(paste0("Genes with multiple independent CpGs (ForwardOLS_multi): ",
			 sum(Allensg_mlv$mlv_method == "ForwardOLS_multi" & Allensg_mlv$cpg_rank == 1), "\n"))
	cat("Method breakdown:\n")
	print(table(Allensg_mlv$mlv_method, useNA = "ifany"))

	cpgs_reg_sig_mlv <- cpgs_reg_sig %>%
	mutate(cpg_gene_type = case_when(
	geneid %in% single_cpg_genes ~ "gene_single_cpg",
	geneid %in% multi_cpg_genes  ~ "gene_multiple_cpgs",
	TRUE                         ~ NA_character_
	)) %>% full_join(Allensg_mlv, by = c("geneid", "cpg"))

	out_file <- paste0(prefix, "_lm_mlv_combined.csv")
	fwrite(cpgs_reg_sig_mlv, file = out_file, sep = ",", quote = FALSE, na = "")
}