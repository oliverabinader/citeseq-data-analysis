# CITE-seq Data Analysis Pipeline

**Author:** Oliver Abinader


## Overview

This repository provides a **reproducible pipeline for processing CITE-seq data** using **10x Genomics Cell Ranger**, followed by downstream multimodal analysis in **Seurat**.

The workflow supports:

* Gene Expression (GEX) libraries
* Feature Barcode (Antibody-derived tags / Protein) libraries
* Batch processing of multiple samples

It generates **filtered feature-barcode matrices** suitable for downstream single-cell analysis (e.g., Seurat, WNN integration).


## Repository Structure

```
.
├── README.md
├── scripts/
│   ├── citeseq_cellranger_count.sh
│   ├── citeseq_rna_adt_integration_pipeline.R
│   └── PI_CD4_Tcell_reclustered_between_conditions_heatmap.R
├── data/
│   ├── GEX_fastqs/
│   ├── Protein_Antibody_fastqs/
│   ├── library_info/
│   └── Feature_ref.csv
```


## Input Data

### 1. FASTQ Files

#### Gene Expression (GEX)

```
data/GEX_fastqs/Folder_C1/
```

Example:

* Folder_C1_3_gene_S10_L001_R1_001.fastq.gz
* Folder_C1_3_gene_S10_L001_R2_001.fastq.gz

#### Antibody / Feature Barcode (ADT)

```
data/Protein_Antibody_fastqs/Folder_C1/
```

Example:

* Folder_C1_F_S26_L001_R1_001.fastq.gz
* Folder_C1_F_S26_L001_R2_001.fastq.gz


### 2. Library Sheet

Each sample is defined using a CSV file:

```
data/library_info/Folder_C1.csv
```

Example:

```csv
fastqs,sample,library_type
GEX_fastqs/Folder_C1,Folder_C1_3_gene,Gene Expression
Protein_Antibody_fastqs/Folder_C1,Folder_C1_F,Antibody Capture
```


### 3. Feature Reference

```
data/Feature_ref.csv
```

Example:

```csv
id,name,read,pattern,sequence,feature_type
Ly6G,Ly6G_TotalSeqB,R2,5PNNNNNNNNNN(BC),ACATTGACGCAACTA,Antibody Capture
CD8a,CD8a_TotalSeqB,R2,5PNNNNNNNNNN(BC),TACCCGTAATAGCGT,Antibody Capture
MHCII,MHCII_TotalSeqB,R2,5PNNNNNNNNNN(BC),GGTCACCAGTATGAT,Antibody Capture
```


## Cell Ranger Processing

### Script

```
scripts/citeseq_cellranger_count.sh
```

### Usage

```bash
bash citeseq_cellranger_count.sh \
  OUTPUT_DIR \
  TRANSCRIPTOME \
  FEATURE_REF \
  LIBRARY_DIR
```

### Arguments

| Argument      | Description                                                |
| ------------- | ---------------------------------------------------------- |
| OUTPUT_DIR    | Base directory for Cell Ranger outputs                     |
| TRANSCRIPTOME | 10x reference transcriptome (e.g. mm10 + custom additions) |
| FEATURE_REF   | Feature barcode reference CSV                              |
| LIBRARY_DIR   | Directory containing sample library CSV files              |

### Workflow

For each sample, the script:

1. Reads library CSV
2. Extracts sample ID
3. Creates output directory
4. Runs `cellranger count`
5. Produces gene expression + ADT matrices


## Output Structure

```
OUTPUT_DIR/
├── Folder_C1/
│   ├── filtered_feature_bc_matrix/
│   ├── raw_feature_bc_matrix/
│   ├── outs/
│   └── BAM files
├── Folder_C2/
└── ...
```


## Requirements

* Cell Ranger (10x Genomics)
* Linux / HPC environment
* Reference transcriptome (e.g. mm10 or custom mm10 + viral genome)
* Feature reference CSV
* R / RStudio (for downstream analysis)


## Downstream Analysis (Seurat + WNN)

### Script

```
scripts/citeseq_rna_adt_integration_pipeline.R
```

### Workflow Summary

#### 1. Seurat Object Creation

* Load RNA (GEX) and ADT matrices
* Create Seurat objects per sample

#### 2. Quality Control

Filters based on:

* `nFeature_RNA`
* `nCount_RNA`
* `nCount_ADT`
* `percent.mt`

#### 3. RNA Processing

* SCTransform normalization
* Integration using SCT anchors

#### 4. ADT Processing

* CLR normalization
* PCA-based dimensional reduction
* Cross-sample integration

#### 5. Multimodal Integration (WNN)

* Weighted Nearest Neighbor (WNN) graph construction
* UMAP embedding
* Clustering (e.g. `wsnn_res.0.3`)

#### 6. Cell Type Annotation

* SingleR classification using ImmGen reference

#### 7. Differential Expression

* `PrepSCTFindMarkers()`
* `FindAllMarkers()`
* Cluster-specific markers

#### 8. Visualization

* UMAP plots
* ADT dot plots


## Additional Analysis: CD4 T Cell Heatmap

### Script

```
scripts/PI_CD4_Tcell_reclustered_between_conditions_heatmap.R
```

### Description

This analysis uses PI-provided CD4 T cell annotations to compare transcriptional programs between conditions:

* **B6 (wild-type)**
* **BLT (knockout)**

### Workflow

#### 1. Data Input

* Seurat integrated object
* External barcode + metadata CSVs

#### 2. Subsetting

* Keeps only annotated CD4 T cells

#### 3. Metadata Mapping

* Assigns condition (B6 vs BLT)
* Assigns cluster labels

#### 4. Expression Aggregation

* Computes average expression per cluster

#### 5. Visualization

* ComplexHeatmap
* Row-scaled (z-score) expression
* Side-by-side comparison (WT vs KO)

### Output

* Cluster-level gene expression heatmap
* Comparative visualization between conditions
* Hierarchical clustering of gene expression patterns


## Notes

* Ensure feature reference matches antibody panel exactly
* Maintain consistent naming between FASTQ, library CSV, and feature reference
* Recommended to validate MD5 checksums for FASTQ integrity prior to processing
