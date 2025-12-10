# ANAT40040 Final Project
ANAT40040 Final Project: Gene Expression Analysis of HER2Amplified Breast Cancer

This repository contains the R workflow used to complete Assignment 2: Gene Expression Analysis and Interpretation for the module ANAT40040 Bio Principles and Cell Organisation.

The analysis uses TCGA BRCA PanCancer Atlas 2018 data downloaded from cBioPortal, and compares ERBB2 amplified vs nonamplified breast tumours.

## Repository Contents
1. ANANT40040_Final_BRCA_ERBB2_Pipeline_RScript.R

This script performs the analysis pipeline:

Data Import & Preprocessing

Loads RNAseq, clinical, and CNA data from cBioPortal

Matches patient barcodes across data types

Defines ERBB2 amplification based on CNA (> 0)

Differential Gene Expression (DESeq2)

Normalises counts with DESeq2

Computes differential expression between Amplified vs Not Amplified

Extracts the top 10 DE genes

Pathway Enrichment Analysis

Uses clusterProfiler, org.Hs.eg.db, KEGG, and Reactome to identify enriched:

GO Biological Processes

KEGG pathways

Reactome pathways

Unsupervised Analysis

Variancestabilised transformation (VST)

PCA plot coloured by ERBB2 status

Heatmap of the top 100 DE genes

Survival Modelling (LASSO Cox)

Builds a risk model using glmnet

Computes risk scores, hazard ratios, and Kaplan–Meier curves

## How to Run the Script

Download the TCGA BRCA dataset from cBioPortal:
https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018

Place the .tar.gz file in the working directory specified in the script.

Open the R script in RStudio.

Install required packages if missing:

install.packages(c("ggplot2", "pheatmap", "survival", "glmnet"))
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "clusterProfiler", "ReactomePA"))


Run the script from top to bottom.

## Summary of Outputs

The script produces:

Top 10 DE genes table

GO, KEGG, and Reactome enrichment plots

PCA plot of VST data

Heatmap of top 100 DE genes

LASSO Cox KM curves & risk model metrics

## Author

James McClatchie
MSc in Artificial Intelligence for Medicine & Medical Research (2026)
