# Deconv-MeWAS

## Overview
This project provides a deconvolution-based framework for identifying cell-type-specific epigenetic and transcriptional mechanisms underlying colorectal cancer (CRC) susceptibility. By integrating bulk colon DNA methylation and transcriptomic datasets with single-cell reference profiles, the pipeline infers cell-type-specific molecular traits and links them to CRC genetic risk through cell-type-specific methylation-wide association studies (ctMWAS) and transcriptome-wide association studies (ctTWAS).
Using normal colon methylation (n = 293) and gene expression (n = 707) datasets together with CRC GWAS summary statistics (78,473 cases and 107,143 controls), the framework prioritizes risk-associated CpG sites, genes, and regulatory pathways at cellular resolution. Integrative analyses combining ctMWAS, ctTWAS, colocalization, and methylation–gene mapping identified high-confidence CRC susceptibility loci, candidate causal genes, and potential therapeutic targets.

![ctMWAS and ctTWAS workflow](./Figures/ctMWAS_ctTWAS.png)

**Step1:** Cell-type-specific deconvolution of bulk DNA methylation and gene expression data

**Step2:** ctMWAS and ctTWAS analyses using genetically predicted molecular traits

**Step3:** Colocalization-based prioritization of candidate causal CpGs and genes

**Step4:** Integrative CpG-to-gene mapping framework

**Step5:** Multi-omics validation using CRC progression datasets and Drug target prioritization analyses

## Methods
### 1. Cell-type-specific deconvolution of bulk DNA methylation and gene expression data
```text
Bulk DNAm (COLONOMICS, GTEx) ──┐
                               ├── EpiSCORE ──> Cell-type proportions
scDNAm Reference (EpiSCORE) ───┘                      │
													  ├─ EpiSCORE + TCA
													  ▼
                                   Cell-type-specific DNAm Profiles


Bulk RNA-seq (BarcUVa-Seq, GTEx) ─┐
                                  ├── CIBERSORTx ──> Cell-type proportions
scRNA-seq Reference (COLON MAP) ──┘                    │
													   ├─ CIBERSORTx (https://cibersortx.stanford.edu/)
													   ▼
                                  Cell-type-specific Expression Profiles
```
**Output:** Cell-type-specific DNA methylation and gene expression matrices for downstream ctMWAS and ctTWAS analyses.
| Analysis | Script |
|----------|---------|
| DNA methylation deconvolution (EpiSCORE + TCA) | `src/Cell_type_specific_expression/EpiSCORE_TCA.R` |


### 2. ctMWAS and ctTWAS analyses using genetically predicted molecular traits

```text
Cell-type-specific DNAm / Expression
                  │
                  ▼
      Molecular Trait Normalization
                  │
                  ▼
      Elastic Net Prediction Models
                  │
                  │  TF-occupied regulatory variants
                  │  (±1 Mb of 51 CRC-associated TFs)
                  ▼
   Genetically Predicted DNAm / Expression
                  │
                  ▼
      Integration with CRC GWAS
      (78,473 cases, 107,143 controls)
                  │
                  ▼
            S-PrediXcan
                  │
          ┌───────┴────────┐
          ▼                ▼
       ctMWAS           ctTWAS
          │                │
          ▼                ▼
     Risk CpGs        Risk Genes
```
| Step | ctMWAS | ctTWAS |
|------|--------|--------|
| Covariate adjustment and PEER correction | `src/ctMWAS/Methy_PEER.R` | `src/ctTWAS/Expr_PEER.R` |
| Prediction model construction | `src/ctMWAS/TWAS_modelbuilding_MethUVA.R` | `src/ctTWAS/TWAS_modelbuildingCRC_EUR.R` |
| Prediction database generation | `src/ctMWAS/WeightsSummaryToDB.R` | `src/ctTWAS/WeightsSummaryToDB.R` |
| Association analysis | `src/ctMWAS/asso.sh` | `src/ctTWAS/asso.sh` |
| Meta-analysis | `src/ctMWAS/Methy_Cauchy_Meta.R` | `src/ctTWAS/Genes_Cauchy_Meta.R` |


### 3. Colocalization-based prioritization of candidate causal CpGs and genes
```text
┌─────────────────────────────┐      ┌─────────────────────────────┐
│            ctMWAS           │      │            ctTWAS           │
├─────────────────────────────┤      ├─────────────────────────────┤
│ Significant CpG Sites       │      │ Significant Genes           │
│           │                 │      │           │                 │
│           ▼                 │      │           ▼                 │
│       mQTL Data             │      │       eQTL Data             │
│           │                 │      │           │                 │
│           ├──────┐          │      │           ├──────┐          │
│           ▼      │          │      │           ▼      │          │
│       coloc.abf  │          │      │       coloc.abf  │          │
│           ▲      │          │      │           ▲      │          │
│           └──────┤          │      │           └──────┤          │
│                  ▼          │      │                  ▼          │
│            CRC GWAS         │      │            CRC GWAS         │
│                  │          │      │                  │          │
│                  ▼          │      │                  ▼          │
│      Colocalized CpGs       │      │      Colocalized Genes      │
│        (PP.H4 > 0.8)        │      │        (PP.H4 > 0.8)        │
└─────────────────────────────┘      └─────────────────────────────┘
```
**Output:** High-confidence CpG sites and genes sharing a causal variant with CRC GWAS signals.
| Coloc | `Functional_analyses/coloc_gwas_xqtl.R` |


### 4. Integrative CpG-to-gene mapping framework
```text
Risk CpGs
    ──► CpG–Gene Association (±1 Mb)
    ──► Independent CpG Selection
    ──► Inverse Regulation Filter (β < 0)
    ──► Candidate Risk Genes
    ──► Oncogene / Tumor Suppressor Classification
```
**Output:** High-confidence CpG–gene regulatory pairs & Putative oncogenes and tumor suppressors
| Genes | `Functional_analyses/sigcpGs_expr_lm_mlv_final.R` |


### 5. Multi-omics validation using CRC progression datasets and Drug target prioritization analyses
```text
Candidate Risk Genes
        │
        ├── TCGA / CPTAC
        │      └── Tumor vs Normal Expression
        │
        ├── DepMap CRISPR Screens
        │      └── Gene Essentiality
        │
        ├── COLON MAP scRNA-seq
        │      └── CRC Progression Analysis
        │
        └── Enrichr
               └── Pathway Enrichment
                    │
                    ▼
         Functionally Supported Risk Genes
         + Biological Pathways
```
**Output:** Functionally supported CRC risk genes; Putative oncogenes and tumor suppressors; Essential genes for CRC cell survival; Cell-type-specific biological pathways
| Genes | `Functional_analyses/enrichR_disease_functional_enrichment_scTWAS_genes.R` |

## Contact
Qing Li: liqingbioinfo@gmail.com \
Xingyi Guo: xingyi.guo@vumc.org
