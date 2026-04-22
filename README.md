# CITE-seq Data Analysis

## Author
Oliver Abinader

---

# Overview

This repository contains a reproducible pipeline for processing **CITE-seq data** using **Cell Ranger**.

The workflow handles:

- Gene Expression (GEX) libraries  
- Feature Barcode (Antibody / Protein) libraries  
- Multiple samples processed in batch  

It generates **filtered feature-barcode matrices** for downstream single-cell analysis (e.g., Seurat).


# Repository Structure

├── README.md
├── scripts/
│ └── citeseq_cellranger_count.sh
├── data/ # (not tracked in Git; example structure below)
│ ├── GEX_fastqs/
│ ├── Protein_Antibody_fastqs/
│ ├── Library_info/
│ └── Feature_ref.csv


# Input Data Structure

- **Lauren_gex_fastqs/** → Gene Expression FASTQ files (GEX)
- **Lauren_protein_antibody_fastqs/** → Antibody (CITE-seq) FASTQs  
- **library_info/** → Library CSV files (one per sample)  
- **Lauren_feature_ref.csv** → Feature barcode reference  


# Example Data

## FASTQ files (GEX)
GEX_fastqs/Folder_C1/
├── Folder_C1_3_gene_06_28_24_S10_L001_R1_001.fastq.gz
├── Folder_C1_3_gene_06_28_24_S10_L001_R2_001.fastq.gz
├── Folder_C1_3_gene_06_28_24_S10_L002_R1_001.fastq.gz
└── Folder_C1_3_gene_06_28_24_S10_L002_R2_001.fastq.gz

## FASTQ files (Antibody / Feature Barcode)
Protein_Antibody_fastqs/Folder_C1/
├── Folder_C1_F_06_28_24_S26_L001_R1_001.fastq.gz
├── Folder_C1_F_06_28_24_S26_L001_R2_001.fastq.gz
├── Folder_C1_F_06_28_24_S26_L002_R1_001.fastq.gz
└── Folder_C1_F_06_28_24_S26_L002_R2_001.fastq.gz

## Library CSV example
library_info/Folder_C1.csv

```csv
fastqs,sample,library_type
GEX_fastqs/Folder_C1,Folder_C1_3_gene,Gene Expression
Protein_Antibody_fastqs/Folder_C1,Folder_C1_F,Antibody Capture
```

## Feature Reference example
Feature_ref.csv
id,name,read,pattern,sequence,feature_type
Ly6G,Ly6G_TotalSeqB,R2,5PNNNNNNNNNN(BC),ACATTGACGCAACTA,Antibody Capture
CD8a,CD8a_TotalSeqB,R2,5PNNNNNNNNNN(BC),TACCCGTAATAGCGT,Antibody Capture
MHCII,MHCII_TotalSeqB,R2,5PNNNNNNNNNN(BC),GGTCACCAGTATGAT,Antibody Capture


# Script

• scripts/citeseq_cellranger_count.sh: Runs Cell Ranger count for each sample using library CSV files.
• Usage: bash scripts/citeseq_cellranger_count.sh OUTPUT_DIR TRANSCRIPTOME FEATURE_REF LIBRARY_DIR
• Arguments:
  - OUTPUT_DIR → base directory for results
  - TRANSCRIPTOME → path to 10x reference genome
  - FEATURE_REF → feature reference CSV
  - LIBRARY_DIR → directory containing library CSV files

• What the script does? For each sample:
  - Reads library CSV
  - Extracts sample ID
  - Creates output directory
  - Runs cellranger count
  - Generates feature-barcode matrices
