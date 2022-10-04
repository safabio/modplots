# modplots package

[![ISSUES](https://img.shields.io/github/issues/safabio/modplots)](https://github.com/safabio/modplots/issues)
[![Stars](https://img.shields.io/github/stars/safabio/modplots?style=social)](https://github.com/safabio/modplots/)


This package contains currently 3 plotting functions for chicken RNA seq data, as well as an annotated gene id /gene name table.

All three functions allow you to add your own gnames table, making it applicable to all species. The table should contain to columns named: "Gene.stable.ID" and "Gene.name".

## Installation

To install modplots run in R:

`devtools::install_github("safabio/modplots", ref="main")`

## mFeatureplot

Allows you to plot any embedding of your Seurat object, and highlight gene expression from any assay. Currently there is no control over the assay used. Please set the active assay beforehand.
Contrary to Seurats FeaturePlot, you have greater control over the color gradient, and it allows you to specify genes via name and not IDs. Also the title of the plot will be the gene's name, instead of its stable ID.

<img src="man/figures/logo.png" align="right" height="20%" width="20%" />

## mDotPlot2

Modified DotPlot function from Seurat. Automatically adds gene names as labels. Input has to be gene IDs for now.

## mVolcanoplot 

Volcano plot function for Seurat and DESeq2 generated DE tables. The backbone of the function is taken from [here](https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/).

Also any other DE table can be given as input, the column names should contain: "gene" or "id" for the ID column, "2F" for the log2 fold change, and "adj" for the adjusted p values. 
Plot can be generated with plotly to create an interactive plot to inspect directly gene names, p val adj and lo2FC's.

## mPCA

Function using plotPCA from DESeq2. You can specify grouping variables and PCs to plot. Additionally you can use shape and label aesthetics to control those ggplot parameters in geom_point(). Automatically returns pca object (prcomp output) to global env.

## mPheatmapDESEq2

Modified pheatmap function for vst transformed dds objects. 

## loadingPlot

Plots the absolute loadings of a prcomp object. At the same time it adds the names of the top n loading genes (by absolute value).

## gnames

Gg6 gene ID and name table. For details use: `?gnames`.
