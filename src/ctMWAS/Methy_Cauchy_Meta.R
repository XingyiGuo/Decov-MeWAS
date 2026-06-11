#conda activate tftwas
library(ACAT)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

cts <- c("Enteriendocrine", "Enterocyte","Goblet","Progenitor")
ct <- cts[as.numeric(args[1])]
print(ct)

x1_filename=paste0("/MethExp/Meth_UVA/Methy_CellTypes/TFTWAS_RES_add_corr/",ct,".TWAS.csv")
x2_filename=paste0("/MethExp/Meth_GTEx/Methy_CellTypes/TFTWAS_RES_add_corr/",ct,".TWAS.csv")
sum1_filename=paste0("/MethExp/Meth_UVA/Methy_CellTypes/TFTWAS_RES_add_corr/",ct,"_model_summaries.csv")
sum2_filename=paste0("/MethExp/Meth_GTEx/Methy_CellTypes/TFTWAS_RES_add_corr/",ct,"_model_summaries.csv")

x1_df_raw=as.data.frame(fread(x1_filename))
x2_df_raw=as.data.frame(fread(x2_filename))
sum1_df=as.data.frame(fread(sum1_filename))
sum2_df=as.data.frame(fread(sum2_filename))

#keep genes with r>0.1
x1_df <- x1_df_raw %>% left_join(sum1_df %>% select(gene, r.global), by = c("gene" = "gene"))
x2_df <- x2_df_raw %>% left_join(sum2_df %>% select(gene, r.global), by = c("gene" = "gene"))

x1_df = filter(x1_df, r.global>0.1)
x2_df = filter(x2_df, r.global>0.1)

#add suffix
colnames(x1_df)=paste0(colnames(x1_df), "_UVA")
colnames(x2_df)=paste0(colnames(x2_df), "_GTEx") 

#rename gene names
names(x1_df)[names(x1_df) == "gene_UVA"] <- "gene"
names(x2_df)[names(x2_df) == "gene_GTEx"] <- "gene"

##merge two EUR dataframe p values toether
x1_x2_df = merge(x1_df, x2_df, by="gene", all = TRUE)
  
#Cannot have NAs in the p-values!, 
x1_x2_df <- x1_x2_df %>%
  rowwise() %>%
  mutate(META_EUR_P = ifelse(is.na(pvalue_UVA), pvalue_GTEx,
         ifelse(is.na(pvalue_GTEx), pvalue_UVA,
                ACAT(c(pvalue_UVA, pvalue_GTEx))))
  )

fwrite(x1_x2_df, paste0("/MethExp/Manuscript_Tables_Figures/Methy_Cauchy/CRC_EUR_GTEX_Methylation_META_",ct,".csv"))
