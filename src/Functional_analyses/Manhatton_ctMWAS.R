library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(data.table)

setwd("D:/GitHub/XingyiGuo/Decov-MeWAS/data/Methy_Cauchy/")

# ===============================
# assign Meta_Tag
# ===============================
assign_meta_tag <- function(ct) {
  
  file <- paste0(root, "CRC_EUR_GTEX_Methylation_META_", ct, "_corr_global.csv")
  df <- fread(file)
  
  meta_nonNA <- df[!is.na(META_EUR_P), ]
  meta_threshold <- 0.05 / nrow(meta_nonNA)
  
  meta_sig <- df %>%
    filter(META_EUR_P < meta_threshold)
  
  df <- df %>%
    mutate(
      Meta_Tag = ifelse(gene %in% meta_sig$gene, "Meta_sig", "Not_Sig"),
      CellTypes = ct
    )
  
  cat(sprintf("CellType: %s | Meta_sig: %d (threshold: %.3e)\n",
              ct, nrow(meta_sig), meta_threshold))
  
  return(df)
}

# ===============================
# Manhattan plot
# ===============================
plot_one <- function(celltype, alltraits, pcut) {
  
  df <- alltraits %>%
    filter(CellTypes == celltype) %>%
    mutate(
      chrom = as.numeric(gsub("chr", "", CHR)),
      bp    = as.numeric(MAPINFO_hg38),
      P     = as.numeric(META_EUR_P)
    ) %>%
    filter(!is.na(chrom) & !is.na(bp) & !is.na(P) & P > 0)
  
  th <- 0.05 / nrow(df)
  
  # ===============================
  # color logic (fixed)
  # ===============================
  df <- df %>%
    mutate(
      sig_type = ifelse(Meta_Tag == "Meta_sig", "Meta_sig", "None"),
      P = pmax(P, 1e-20),
      
      plot_color = case_when(
        sig_type == "Meta_sig" ~ paste0("Meta_", CellTypes),
        chrom %% 2 == 0 ~ "even",
        TRUE ~ "odd"
      )
    )
  
  # cumulative position
  chr_info <- df %>%
    group_by(chrom) %>%
    summarise(chr_len = max(bp), .groups = "drop") %>%
    mutate(tot = cumsum(chr_len) - chr_len)
  
  df <- df %>%
    left_join(chr_info, by = "chrom") %>%
    arrange(chrom, bp) %>%
    mutate(BPcum = bp + tot)
  
  axisdf <- df %>%
    group_by(chrom) %>%
    summarise(center = mean(BPcum))
  
  label_df <- df %>%
    filter(P < pcut, !is.na(Genes)) %>%
    group_by(Genes) %>% 
    slice_min(P, n = 1) %>%   # keep top CpG per gene (avoid clutter)
    ungroup()
  
  # ===============================
  # plot
  # ===============================
  p <- ggplot(df, aes(x = BPcum, y = -log10(P))) +
    
    geom_point(aes(color = plot_color), size = 0.6, alpha = 0.8) +
    
    scale_color_manual(values = c(
      "even" = "#D4D4D4", 
      "odd"  = "#999999", 
      
      "Meta_Enteriendocrine" = "#F27684",
      "Meta_Enterocyte"      = "#cc8de2",
      "Meta_Goblet"          = "#3869af",
      "Meta_Progenitor"      = "#99b2d6"
    )) +
    
    scale_x_continuous(breaks = axisdf$center, labels = axisdf$chrom) +
    
    geom_hline(yintercept = -log10(th), linetype = "dashed") +
    
    coord_cartesian(ylim = c(0, 20)) +
    
    theme_bw() +
    theme(
      legend.position = "none",   # ✅ remove legend
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    
    xlab("Chromosome") +
    ylab("-log10(P value)") +
    ggtitle(paste0(celltype, " cell type")) +
    theme(plot.title = element_text(size = 8, hjust = 0.5))+
    geom_text_repel(
      data = label_df,
      aes(label = Genes),
      size = 2,
      max.overlaps = 50,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.color = "black",   # ✅ black line
      segment.size = 0.4,        # ✅ thickness
      segment.alpha = 0.8,       # optional transparency
      min.segment.length = 0,    # ✅ ALWAYS draw line (important)
      force = 1                  # adjust repulsion strength
    )
  
  # save
  outfile <- file.path(outdir, paste0(celltype, "_CpG_Manhattan.tiff"))
  ggsave(outfile, p, width = 8, height = 5, dpi = 300)
  
  return(p)
}

# ===============================
# Main analyses
# ===============================

root = "D:/GitHub/XingyiGuo/Decov-MeWAS/data/Methy_Cauchy/"
outdir = "D:/GitHub/XingyiGuo/Decov-MeWAS/data/Methy_Cauchy/"

cts = c("Enteriendocrine", "Enterocyte", "Goblet", "Progenitor")
cts_methy <- lapply(cts, assign_meta_tag) %>% bind_rows()

# ===============================
# CpG positions
# ===============================
gtex_betas <- fread(paste0(root, "GSE213478_methylation_DNAm_noob_final_BMIQ_colon_224_BETA.hg38.csv"))
uva_betas  <- fread(paste0(root, "CLX_methylation_Normal_Mucosa_BETA_BMIQ.hg38.csv"))

cpg_pos <- bind_rows(
  uva_betas[, .(cpg, CHR, MAPINFO_hg38)],
  gtex_betas[, .(cpg, CHR, MAPINFO_hg38)]
) %>% distinct()

alltraits <- merge(
  cpg_pos, cts_methy,
  by.x = "cpg", by.y = "gene",
  all.y = TRUE
) %>% distinct()

# ===============================
# CpG-gene pairs
# ===============================
cpg_gene_pairs <- read_excel("D:/GitHub/XingyiGuo/Decov-MeWAS/data/Supplementary_Tables.xlsx", sheet ="S6", skip=2)

alltraits$cpg_celltypes <- paste(alltraits$cpg, alltraits$CellTypes, sep = "_")
cpg_gene_pairs$cpg_celltypes <- paste(cpg_gene_pairs$CpGs, cpg_gene_pairs$CellTypes, sep = "_")

alltraits_gene_annot <- merge(
  alltraits,
  cpg_gene_pairs[, c("cpg_celltypes", "Genes")],
  by = "cpg_celltypes",
  all.x = TRUE
) %>% distinct()

# fwrite(alltraits_gene_annot, file = "alltraits_gene_annot.csv")

# ===============================
# plot
# ===============================
pcut=1

plots <- lapply(cts, function(ct) {
  plot_one(ct, alltraits_gene_annot, pcut)
})


