# ECE204-Project2
Project repository for ECE 204
This report delves into analysis of Head and neck squamous cancer.

## Languages:
R , C++, Python, Excel

## Packages and installation:
library(GenomicFeatures)

if (!requireNamespace("BiocManager", quietly = TRUE)) <br />
    >> install.packages("BiocManager") <br />
BiocManager::install("biomaRt")

library(openxlsx) <br />
install.packages("org.Hs.eg.db") <br />
install.packages("dataset") <br />
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## Description
The code is available in the clean_data.R file. In this we have compared the mutations as obtained on TCGA dataset from the dN/dS package developed by authors of the paper " Universal Patterns of Selection in Cancer and Somatic Tissues" and our analysis which involves calculating the non synonymous and synonymous mutations of the coding mutations and calculating the dN/dS values for the genes from that. The comparisions of the dN/dS values are done through t-test, z-test, and chi-squeared test in the code. <br />
The plot "Chart (1).png" shows the distribution of the dN/dS values for all non-synonymous genes where infinite values have been excluded. 
