# STARnet
**STARnet** is a spatiotemporal analytical framework designed to systematically uncover drug-induced functional recovery within complex disease networks by integrating single-cell and spatial transcriptomic data with network pharmacology.
<img width="2073" height="1678" alt="figure--1" src="https://github.com/user-attachments/assets/57d3b751-1913-4656-9f9c-e4828fcca509" />
# Introduction to STARnet
STARnet is a spatiotemporal analytical framework for systematically uncovering drug-induced functional recovery in complex diseases. It integrates single-cell and spatial transcriptomic data with network pharmacology modeling to enable unified, multi-scale, and multi-modal analysis.
At spatiotemporal resolution, STARnet quantitatively characterizes drug-driven recovery patterns in biological networks, identifying key target cell populations and gene sets. In addition, it reveals synergistic relationships among cell types and target genes, highlighting critical genes, ligand–receptor pairs, and cell populations that cooperate within recovery-associated networks and pathways.
STARnet has been validated on multiple real-world biological datasets, demonstrating high consistency with previous findings.
# Data Input
Single-cell transcriptomic data, or single-cell data integrated with spatial transcriptomic data.
## Single-cell transcriptomic data
If the input data consist of only single-cell transcriptomic data, the folder structure should follow the standard 10x Genomics format. For example:
```text
10x
└─GSEXXXXXX
├─control
│ barcodes.tsv.gz
│ features.tsv.gz
│ matrix.mtx.gz
│
├─model  
│ barcodes.tsv.gz
│ features.tsv.gz
│ matrix.mtx.gz
│
├─drug  
│ barcodes.tsv.gz
│ features.tsv.gz
│ matrix.mtx.gz
```
## Single-cell + Spatial Transcriptomics data
If the input data include both single-cell transcriptomic data and spatial transcriptomic data, STARnet requires that the single-cell data first be mapped to spatial coordinates using CellTreck.  
The input folder should contain the following files:
```
├─control
│ coord1_.csv
│ spatial__adata.csv
│
├─model  
│ coord1_.csv
│ spatial__adata.csv
│
├─drug  
│ coord1_.csv
│ spatial__adata.csv
```
# Requirements
```
python==3.10.16
adjustText==1.3.0
adjustText==1.3.0
anndata==0.12.6
bokeh==3.7.3
gseapy==1.0.3
harmonypy==0.0.10
holoviews==1.22.1
liana==1.6.1
matplotlib==3.7.2
networkx==3.4.2
numpy==2.3.5
pandas==2.3.3
panel==1.8.4
Requests==2.32.5
scanpy==1.11.5
scikit_learn==1.7.2
scipy==1.16.3
seaborn==0.13.2
umap_learn==0.5.7
```
