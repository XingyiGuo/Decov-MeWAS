#conda activate tftwas
library(ACAT)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

ts="colontrans"
CTlist=c("ABS","CT","EE","GOB","STM", "TAC","TUF")
CT=CTlist[as.numeric(args[1])]

x1_filename=paste0("/data/l2_bioinfo1/liq17/MethExp/Expr_CRC_EUR/hires54k_crc_eur_genes10k_SingleCellCount_bulkCPM/TFTWAS_RES_add_corr/CRC_EUR_",CT,".TWAS")
x2_filename=paste0("/data/l2_bioinfo1/liq17/MethExp/Expr_GTEx_SC_Bulk/hires54k_gtex_colontrans_genes10k_SingleCellCount_BulkCPM/TFTWAS_RES_add_corr/GTEX_colontrans_",CT,".TWAS")
#x3_filename=paste0("/data/sbcs/GuoLab/backup/liq17/MethExp/Expr_CRC_ASIAN/hires54k_crc_asian_genes10k_SingleCellCount_BulkCPM/TFTWAS_RES/CRC_ASIAN_",CT,".TWAS")
sum1_filename=paste0("/data/l2_bioinfo1/liq17/MethExp/Expr_CRC_EUR/hires54k_crc_eur_genes10k_SingleCellCount_bulkCPM/TFTWAS_RES_add_corr/CRC_EUR_",CT,"_model_summaries.csv")
sum2_filename=paste0("/data/l2_bioinfo1/liq17/MethExp/Expr_GTEx_SC_Bulk/hires54k_gtex_colontrans_genes10k_SingleCellCount_BulkCPM/TFTWAS_RES_add_corr/GTEX_colontrans_",CT,"_model_summaries.csv")

x1_df_raw=as.data.frame(fread(x1_filename))
x2_df_raw=as.data.frame(fread(x2_filename))
#x3_df=as.data.frame(fread(x3_filename))
sum1_df=as.data.frame(fread(sum1_filename))
sum2_df=as.data.frame(fread(sum2_filename))

#keep genes with r>0.1
x1_df <- x1_df_raw %>% left_join(sum1_df %>% select(gene, r.global), by = c("gene" = "gene"))
x2_df <- x2_df_raw %>% left_join(sum2_df %>% select(gene, r.global), by = c("gene" = "gene"))
x1_df = filter(x1_df, r.global>0.1)
x2_df = filter(x2_df, r.global>0.1)
#x3_df = filter(x3_df, pred_perf_r2>0.01)

#add suffix
colnames(x1_df)=paste0(colnames(x1_df), "_UVA")
colnames(x2_df)=paste0(colnames(x2_df), "_GTEx") 
#colnames(x3_df)=paste0(colnames(x3_df), "_ASIAN")

#rename gene names
names(x1_df)[names(x1_df) == "gene_UVA"] <- "gene"
names(x2_df)[names(x2_df) == "gene_GTEx"] <- "gene"
#names(x3_df)[names(x3_df) == "gene_name_ASIAN"] <- "gene"

##merge two EUR dataframe p values toether
x1_x2_df = merge(x1_df, x2_df, by="gene", all = TRUE)
  
#Cannot have NAs in the p-values!, 
x1_x2_df <- x1_x2_df %>%
  rowwise() %>%
  mutate(META_EUR_P = ifelse(is.na(pvalue_UVA), pvalue_GTEx,
         ifelse(is.na(pvalue_GTEx), pvalue_UVA,
                ACAT(c(pvalue_UVA, pvalue_GTEx))))
  )

# # merge EURs with Asian
# x1_x2_x3_df = merge(x1_x2_df, x3_df, by="gene_name", all = TRUE)

# x1_x2_x3_df <- x1_x2_x3_df %>%
  # rowwise() %>%
  # mutate(META_EUR_ASIAN_P = ifelse(is.na(META_EUR_P), pvalue_ASIAN,
                             # ifelse(is.na(pvalue_ASIAN), META_EUR_P,
                                    # ACAT(c(pvalue_ASIAN, META_EUR_P))))
  # )

fwrite(x1_x2_df, paste0("/data/l2_bioinfo1/liq17/MethExp/Manuscript_Tables_Figures/Expr_Cauchy/CRC_EUR_ASIAN_GTEX_CountCPM_META_",CT,".csv"))
