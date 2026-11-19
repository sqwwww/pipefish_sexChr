# Turnovers of sex chromosomes and sex differentiation pathways across pipefishes and seahorses

This repository contains the analysis scripts and pipelines associated with the study **“Turnovers of sex chromosomes and sex differentiation pathways across pipefishes and seahorses.”**

The code covers genome assembly and evaluation, sex chromosome analyses, testicular single-cell transcriptomics, chromosome-end GC-content analyses, and comparative tissue- and cell-type-level gene expression analyses.

## Repository structure

| Directory | Description |
|---|---|
| `01.genomeAssembly/` | Genome assembly and quality assessment pipelines for *Microphis brachyurus* (Mbr), *Hippichthys heptagonus* (Hhe), and *Dunckerocampus dactyliophorus* (Dda). |
| `02.sex_chr_analysis/` | Sex chromosome analyses for Mbr, Hhe, *Corythoichthys intestinalis* (Cin), and Dda. |
| `03.testis_SC_analysis/` | Testicular single-cell RNA-seq analyses, including species-level processing, cross-species integration, and RNA velocity analyses. |
| `04.GCspike_analysis/` | Calculation and visualization of GC-content enrichment near chromosome ends. |
| `05.tissue_and_cellType_analysis/` | Comparative analyses of gonadal tissue-expression bias and testicular cell-type expression conservation and turnover, including analyses of selected sex-related genes. |

## Usage

Scripts are organized by analysis and should be run with the corresponding reference genomes, sequencing data, expression matrices, and gene-orthology tables described in the manuscript and its supplementary materials. Several scripts use project-specific input and output paths; these should be adjusted to the user's local directory structure before running.

Required software and R packages vary by analysis. Please consult the individual scripts for tool names and analysis parameters.


