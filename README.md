# stCancer

## Introduction

The `stCancer` package focuses on processing and analysis spatial transcriptome data for cancer research. Except basic data processing steps, stCancer takes more considerations to spatial information and cancer-specific features.

The workflow of `stCancer` mainly consists of two modules: `stStatistics` and `stAnnotation`.
* The `stStatistics` performs quality control and basic statistical analyses.
* The `stAnnotation` performs functional data analyses and visualization of single sample, including low-dimensional representation, clustering, gene expression pattern analysis, copy number variants estimation, ligand-receptor interaction analysis, phenotype heterogeneity analysis, etc.

After all the computational analyses finished, detailed and graphical reports were automatically generated in user-friendly HTML format.

## System Requirements
* R version: >= 4.0.0

## Current version
* stCancer 0.1.0

## Installation

```
checkPkg <- function(pkg){
return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")

library(devtools)
if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")
if(!checkPkg("markdown")) install.packages("markdown")
if(!checkPkg("R.utils")) install.packages("R.utils")
if(!checkPkg("NNLM")) install_github("linxihui/NNLM")
if(!checkPkg("copykat")) install_github("navinlabcode/copykat")
if(!checkPkg("Seurat")) install.packages("Seurat")
if(!checkPkg("knitr")) BiocManager::install("knitr")
if(!checkPkg("GSVA")) BiocManager::install("GSVA")
if(!checkPkg("pheatmap")) BiocManager::install("pheatmap")
if(!checkPkg("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
if(!checkPkg("ClusterProfiler")) BiocManager::install("ClusterProfiler")
if(!checkPkg("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if(!checkPkg("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
if(!checkPkg("glmGamPoi")) BiocManager::install('glmGamPoi')

install_github("Miaoyx323/stCancer")
```

## Usage

`stCancer` takes the output of spaceranger as input and outputs the preliminary analysis results of samples. 

```
library(stCancer)

sample.name <- "ST_data"
data.path <- "Path/to/the/outputs/of/Spaceranger"  # including 'filtered_feature_bc_matrix' and 'spatial' folders
save.path <- "the/save/path"

results <- stStatistics(sample.name, 
                        dataPath = data.path, 
                        savePath = save.path, 
                        species = "human")
object <- results$object

results <- stAnnotation(object, 
                        savePath = save.path, 
                        species = "human", 
                        bool.NMF = T, 
                        bool.CellType = T, 
                        bool.CNV = T, 
                        bool.interaction = T, 
                        bool.tumor.feature = T)
object <- results$object
```


## Data
Here, we provide some example data of [HCC](http://lifeome.net/supp/livercancer-st/data.htm) from [Comprehensive analysis of spatial architecture in primary liver cancer](https://www.science.org/doi/10.1126/sciadv.abg3750). 

<!--
You can download them and run scripts in [wiki](https://github.com/Miaoyx323/stCancer/wiki) to understand the workflow of stCancer.
-->

## Citation
Please use the following citation:

Chen Z, Miao Y, Tan Z, et al. scCancer2: data-driven in-depth annotations of the tumor microenvironment at single-level resolution[J]. Bioinformatics, 2024, 40(2): btae028.

RUI WU, WENBO GUO, XINYAO QIU,et al. Comprehensive analysis of spatial architecture in primary liver cancer. SCIENCE ADVANCES, 2021. DOI: 10.1126/sciadv.abg3750
https://www.science.org/doi/10.1126/sciadv.abg3750

