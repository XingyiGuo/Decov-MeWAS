suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library("coloc"))
suppressMessages(library(purrr))
suppressMessages(library(readxl))
suppressMessages(library(arrow))
args <- commandArgs(trailingOnly = TRUE)

################################
####Step 1: Define function####
################################

coloc_TwoSamples_SS <- function(cell_qtl_df_pvar_df, gwas_df, cell, out_root, out_folder, flank_length=500000){  # ✓ 去掉末尾逗号
    stash <- matrix(integer(0), nrow = 0, ncol = 8) %>% as.data.frame()
    names(stash) <- c(
        "cell", "phenotype_id",
        "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    total_items <- unique(cell_qtl_df_pvar_df[["phenotype_id"]])
    for(i in 1:length(total_items)){
        print(paste0(i, "/", length(total_items)))
        # get the i-th item
        pheno_item <- total_items[i]                                                        # ✓ 变量名统一
        pheno_item_qtl_df <- cell_qtl_df_pvar_df[cell_qtl_df_pvar_df[["phenotype_id"]] == pheno_item, ]
        # get the overlap
        common_rs <- intersect(pheno_item_qtl_df[["variant_id"]], gwas_df[["SNP"]])
        pheno_item_qtl_df_common <- pheno_item_qtl_df[pheno_item_qtl_df[["variant_id"]] %in% common_rs, ]
        gwas_df_common           <- gwas_df[gwas_df[["SNP"]] %in% common_rs, ]
        # rename columns
        colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "variant_id"] <- "SNP"
        colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "REF"]        <- "A2"
        colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "ALT"]        <- "A1"
		colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "slope"]        <- "BETA"
		colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "slope_se"]        <- "SE"
		colnames(pheno_item_qtl_df_common)[colnames(pheno_item_qtl_df_common) == "pval_nominal"]        <- "P"
        # Match xqtls and GWAS betas
        match_xqtl_gwas <- matching_xqtls_and_gwas(gwas_df_common, pheno_item_qtl_df_common)
		match_xqtl_gwas <- match_xqtl_gwas[order(match_xqtl_gwas$P_g), ] 
		match_xqtl_gwas <- match_xqtl_gwas[!duplicated(match_xqtl_gwas$SNP), ] # remove dup
		match_xqtl_gwas <- match_xqtl_gwas[!is.na(match_xqtl_gwas$BETA_x) & # remove na
                                    !is.na(match_xqtl_gwas$SE_x)   &
                                    !is.na(match_xqtl_gwas$P_x)    &
                                    !is.na(match_xqtl_gwas$maf),   ]
		match_xqtl_gwas <- match_xqtl_gwas[!is.na(match_xqtl_gwas$maf) & # maf (0,1)
                                     match_xqtl_gwas$maf > 0    &
                                     match_xqtl_gwas$maf < 1,   ]
        # Conduct coloc when >= 50 variants overlap
        coloc_res <- rep(NA, 6)                                                             # ✓ 与 result$summary 长度一致
        if(nrow(match_xqtl_gwas) >= 50){
            ds1 <- list(beta     = match_xqtl_gwas$BETA_g,
                        varbeta  = (match_xqtl_gwas$SE_g)^2,
                        pvalues  = match_xqtl_gwas$P_g,
                        snp      = match_xqtl_gwas$SNP,
                        type     = "cc",
                        N        = 180000)
            ds2 <- list(beta     = match_xqtl_gwas$BETA_x,
                        varbeta  = (match_xqtl_gwas$SE_x)^2,
                        pvalues  = match_xqtl_gwas$P_x,
                        snp      = match_xqtl_gwas$SNP,
                        type     = "quant",
                        N        = qtl_sample_size)                                         # ✓ 去掉多余括号
            result    <- coloc.abf(ds1, ds2, MAF = match_xqtl_gwas$maf)
            coloc_res <- unlist(result$summary)
            if(coloc_res["PP.H4.abf"] > 0.5){                                              # ✓ 用名字索引更安全
                fwrite(result$result,
                       paste0(out_root, "/", out_folder, "/", pheno_item, ".coloc.details.csv"))  # ✓ 去掉空格
            }
        } else {
            cat("Skipping", pheno_item, ": less than 50 overlapping SNPs.\n")              # ✓ 用 pheno_item
        }
        stash[nrow(stash) + 1, ] <- c(cell, pheno_item, coloc_res)                               # ✓ 列数匹配 (1+6=7)
    }
    fwrite(stash, paste0(out_root, "/", out_folder, "/", cell,"coloc.", as.character(flank_length / 1000), "K.sum.csv"))                  # ✓ 去掉空格
}


################################
####Step 2: Define variables####
################################
source("/MethExp/src/coloc/coloc_gwas_xqtl_support.R")

flank_length <- 500000
qtls_list <- c("mQTLs", "eQTLs")
resources_list <- c("UVA", "GTEx")
mqtls_cells_list <- c("Enteriendocrine","Enterocyte","Goblet","Progenitor")
eqtls_cells_list <- c("ABS", "GOB", "STM")

qtl <- qtls_list[as.numeric(args[1])]
resource <- resources_list[as.numeric(args[2])]

###load GWAS
gwas_folder = "/data/l2_bioinfo1/liq17/Cancer_GWAS_SS/CRC_GWAS_SS/MetaXcan-UK-US-180K/"
gwas_list <- vector("list", 22)
for(chrom in 1:22){
    gwas_list[[chrom]] <- fread(paste0(gwas_folder, "chr", chrom, ".assoc.dosage"))
}
gwas_df <- rbindlist(gwas_list) 
gwas_df$V1 <- NULL

###load qtls path
qtl_path=""
qtl_sample_size=NA
if(qtl=="mQTLs"){
	if(resource=="UVA"){
		qtl_path="/MethExp/Meth_UVA/Methy_CellTypes/tensor_sigCpGs_cis500kqtls/"
		qtl_sample_size=132
	}else if(resource=="GTEx"){
		qtl_path="/MethExp/Meth_GTEx/Methy_CellTypes/tensor_sigCpGs_cis500kqtls/"
		qtl_sample_size=161
	}else{
		print("Invalide inputs for qtl or resource type")
	}
}else if(qtl=="eQTLs"){
	if(resource=="UVA"){
		qtl_path="/MethExp/Expr_CRC_EUR/hires54k_crc_eur_genes10k_SingleCellCount_bulkCPM/tensor_sigCpGs_cis500kqtls/"
		qtl_sample_size=423
	}else if(resource=="GTEx"){
		qtl_path="/MethExp/Expr_GTEx_SC_Bulk/hires54k_gtex_colontrans_genes10k_SingleCellCount_BulkCPM/tensor_sigCpGs_cis500kqtls/"
		qtl_sample_size=284
	}else{
		print("Invalide inputs for qtl or resource type")
	}
}

###load pvar, pvar contains REF(A2), ALT(A1), which can be used to match A1,A2 with GWAS SS
UVA_meth_pvar_list <- vector("list", 22)
for(chrom in 1:22){ UVA_meth_pvar_list[[chrom]] <- fread(paste0("/MethExp/Meth_UVA/SNPsImputation2019/plink/hg38/hg19_annote_hg19_liftover_hg38/chr",chrom,".dose.filter.recode.R03.bgzip.hg38.rs.pvar"))}
UVA_meth_pvar_df <- rbindlist(UVA_meth_pvar_list) 

UVA_expr_pvar_list <- vector("list", 22)
for(chrom in 1:22){ UVA_expr_pvar_list[[chrom]] <- fread(paste0("/MethExp/Genotype_CRC_EUR/R03_annot_with_hg38/chr",chrom,".dose.filter.R03.hg38.vcf.rs.pvar"))}
UVA_expr_pvar_df <- rbindlist(UVA_expr_pvar_list) 

GTEx_pvar_list <- vector("list", 22)
for(chrom in 1:22){ GTEx_pvar_list[[chrom]] <- fread(paste0("/MethExp/Genotype_GTEx/genotype_wb/vcf/annot_with_hg38/gtexV8_WGS_EUR_blood_genotype_PASS.chr",chrom,".rs.hg38.pvar"))}
GTEx_pvar_df <- rbindlist(GTEx_pvar_list) 

pvar_df=data.frame() # pvar_df has five columns: #CHROM   POS ID    REF    ALT
if((qtl=="mQTLs") & (resource=="UVA")){
	pvar_df=UVA_meth_pvar_df
}else if ((qtl=="eQTLs") & (resource=="UVA")){
	pvar_df=UVA_expr_pvar_df
}else if(resource=="GTEx"){
	pvar_df=GTEx_pvar_df
}

### For each cell load parquet files
target_cell_list<-c()
if(qtl=="mQTLs"){
	target_cell_list=mqtls_cells_list
}else if(qtl=="eQTLs"){
	target_cell_list=eqtls_cells_list
}


for(cell in target_cell_list){
    cell_qtl_list <- vector("list", 22)
    # load qtl from each chrom
    for(chrom in 1:22){
        parquet_file <- paste0(qtl_path, cell, "/", resource, ".cis_qtl_pairs.", chrom, ".parquet")
        if(file.exists(parquet_file)){                          # ✓ R 语法
            chr_qtl_df <- as.data.frame(read_parquet(parquet_file))
            cell_qtl_list[[chrom]] <- chr_qtl_df
        }
    }
    cell_qtl_df <- rbindlist(cell_qtl_list, fill = TRUE)
    # Add rs, REF, ALT from pvar
    cell_qtl_df_pvar_df <- merge(cell_qtl_df, 
                              pvar_df[, c("ID", "REF", "ALT")],
                              by.x = "variant_id",   # cell_qtl_df 里的列名
                              by.y = "ID",           # pvar_df 里的列名
                              all.x = TRUE)
    # Conduct coloc
	print(paste0("###coloc analyses ", cell, " items  ", length(unique(cell_qtl_df_pvar_df["phenotype_id"]))))
    coloc_TwoSamples_SS(cell_qtl_df_pvar_df, gwas_df, cell, out_root="/MethExp/coloc_mqtl_eqtl_gwas/", out_folder=paste0(qtl, "_", resource))
}

