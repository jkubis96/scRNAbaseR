### scRNAbaseR - R library for single cell analyses

<br />



<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="180" />
</p>


### Author: Jakub Kubiś

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Laboratory of Single Cell Analyses<br />

<br />

<p align="left">
<img  src="fig/lsca.png" alt="drawing" width="250" />
</p>
</div>


<br />


## Description

A comprehensive toolkit designed for the basic analysis of single-cell RNA sequencing (scRNA-seq) data without cell barcodes, relying solely on unique molecular identifiers (UMIs). Developed within the Laboratory of Single-Cell Analyses at the Institute of Bioorganic Chemistry, PAS, in Poznań, Poland. 

Core functionalities include:
- Merging matrices from featureCounts.
- Quality control of cell content to identify and exclude low-quality data.
- Principal Component Analysis (PCA) for dimensionality reduction.
- Uniform Manifold Approximation and Projection (UMAP) for visualization.
- Clustering for group identification within datasets.
- Marker selection to identify key features for cell types or states.
- Visualization tools to facilitate exploration and presentation of results.

scRNAbaseR empowers researchers to perform essential preprocessing, analysis, and visualization tasks on scRNA-seq and RNA-seq datasets in a straightforward and efficient manner.



#### Installation

```
install.packages("https://github.com/jkubis96/scRNAbaseR/raw/refs/heads/main/scRNAbaseR_0.1.0.tar.gz", repos = NULL, type = "source")
```


#### Loading

```
library(scRNAbaseR)
```


#### Documentation

* [scRNAbaseR](https://jkubis96.github.io/scRNAbaseR/index.html)

<br />

<br />


#### Have fun JBS©