#install.packages("remotes")
#install.packages("data.table")
#remotes::install_github("wjawaid/enrichR")

library(data.table)
library("remotes")
library("readxl")
library("enrichR")
library("dplyr")
library(stringr)
listEnrichrSites()
setEnrichrSite("Enrichr")

format_output <- function(enriched, output_dir ,ouptut_prefix){
  header=c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
  # Initialize empty list to collect filtered data frames
  filtered_list <- list()
  # Loop through each enrichment result
  for (sheet_name in names(enriched)) {
    df <- enriched[[sheet_name]]
    # Check if required column exists
    if ("Adjusted.P.value" %in% colnames(df)) {
      df_sig <- df[df$Adjusted.P.value < 0.1, ]
      
      # Proceed only if non-empty
      if (nrow(df_sig) > 0) {
        df_sig$Source <- sheet_name
        filtered_list[[sheet_name]] <- df_sig
      }
    }
  }
  # Combine all filtered data frames into one
  if (length(filtered_list) > 0) {
    combined_df <- do.call(rbind, filtered_list)
    # Reorder columns if needed
    combined_df <- combined_df[, c(header, "Source")]
    # Save to Excel
    fwrite(combined_df, 
           file = paste0(output_dir, "/",ouptut_prefix,"_Pathway_enriched_genes.tsv"), sep="\t")
  } else {
    message("No significant results (Adjusted.P.value < 0.1) found in any enrichment set.")
  }
}


dbs_lib <- c("Reactome_Pathways_2024", 
             "KEGG_2021_Human", 
             "MSigDB_Hallmark_2020",
             "GO_Biological_Process_2025",
             "GO_Cellular_Component_2025",
             "GO_Molecular_Function_2025")

################################
#### ctMWAS results EnrichR#####
################################

wk="/Decov-MeWAS/data/"

genes_df <- read_excel(paste0(wk, "Supplementary_Tables.xlsx"), sheet="S6", skip=2)
genes_df <- genes_df[genes_df$'Regression results'!="", ]
input <- unique(genes_df$Genes)
print(length(input))
enriched <- enrichr(input, dbs_lib)
format_output(enriched, paste0(wk, "/EnrichR/"), "MWAS_genes")

##Group by cell type
genes_df <- read_excel(paste0(wk, "Supplementary_Tables.xlsx"), sheet="S6", skip=2)
for(model in unique(genes_df$CellTypes)){
  model_df = genes_df[genes_df$CellTypes==model,]
  input <- unique(model_df$Genes)
  print(length(input))
  enriched <- enrichr(input, dbs_lib)
  format_output(enriched,paste0(wk, "/EnrichR/"), paste0("MWAS", model))
}

################################
#### ctTWAS results EnrichR#####
################################
dbs_lib <- c("Reactome_Pathways_2024", "KEGG_2021_Human", "MSigDB_Hallmark_2020")
genes_df <- read_excel(paste0(wk, "Supplementary_Tables.xlsx"), sheet="S7", skip=2)
genes_df <- dplyr::filter(genes_df, PP.H4.abf > 0.8)

input <- unique(genes_df$Genes)
print(length(input))
enriched <- enrichr(input, dbs_lib)
format_output(enriched,paste0(wk, "/EnrichR/"), "TWAS_genes")

##Group by cell type
genes_df <- read_excel(paste0(wk, "Supplementary_Tables.xlsx"), sheet="S7", skip=2)
genes_df <- dplyr::filter(genes_df, PP.H4.abf > 0.8)
for(model in unique(genes_df$CellTypes)){
  model_df = genes_df[genes_df$CellTypes==model,]
  input <- unique(model_df$Genes)
  print(length(input))
  enriched <- enrichr(input, dbs_lib)
  format_output(enriched,paste0(wk, "/EnrichR/"), paste0("TWAS", model))
}


# 示例数据
enrich_df <- enrich_df_raw %>% filter(str_detect(Model, "NL"))

# 添加一个列：-log10 FDR
enrich_df$logFDR <- -log10(enrich_df$`Adjusted.P.value (FDR)`)

# 绘制气泡图
p <- ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
  geom_point(aes(size = Count, color = logFDR)) +
  scale_size_continuous(range = c(1, 6)) +   # reduced bubble size
  scale_color_gradient(low = "navy", high = "magenta") +
  labs(
    x = "Gene Ratio",
    y = "",
    size = "Gene Count",
    color = "-log10(FDR)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))



enrich_df$Term <- factor(enrich_df$Term,
                         levels = rev(unique(enrich_df$Term)))
p <- ggplot(enrich_df, aes(x = GeneRatio, y = Term)) +
  geom_point(aes(size = Count, color = logFDR)) +
  scale_size_continuous(range = c(1, 6)) +
  scale_color_gradient(low = "navy", high = "magenta") +
  labs(
    x = "Gene Ratio",
    y = "",
    size = "Gene Count",
    color = "-log10(FDR)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8))

ggsave("/ms_figures_tables/Figures/figures_before_merge/enrich_bubble_NL.pdf", plot = p, width = 8, height = 12)



