
<!-- badges: start -->
[![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/APL.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/APL)
[![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/APL.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/APL)
[![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/APL.svg)](http://bioconductor.org/packages/stats/bioc/APL/)
[![Bioc support](https://bioconductor.org/shields/posts/APL.svg)](https://support.bioconductor.org/tag/APL)
[![Bioc history](https://bioconductor.org/shields/years-in-bioc/APL.svg)](https://bioconductor.org/packages/release/bioc/html/APL.html#since)
[![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/APL.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/APL/)
[![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/APL.svg)](https://bioconductor.org/packages/release/bioc/html/APL.html#since)
<!-- badges: end -->

<img src="man/figures/fig_AP.png" width="700">

# APL

`APL` is a package developed for computation of Association Plots, a method for visualization and analysis of single cell transcriptomics data. The main focus of `APL` is the identification of genes characteristic for individual clusters of cells from input data. 

When working with `APL` package please cite:
```
Gralinska, E., Kohl, C., Fadakar, B. S., & Vingron, M. (2022). 
Visualizing Cluster-specific Genes from Single-cell Transcriptomics Data Using Association Plots. 
Journal of Molecular Biology, 434(11), 167525.
```

## Installation

The `APL` can be installed from GitHub:
    
    library(devtools)
    install_github("VingronLab/APL")
    

To additionally build the package vignette, run instead:

    install_github("VingronLab/APL", build_vignettes = TRUE, dependencies = TRUE)


Building the vignette will however take considerable time.

**The vignette can also be found under the link: https://vingronlab.github.io/APL/ (hyperlink in the GitHub repository description).**

To install the `APL` from Bioconductor, run:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("APL")

 
## Pytorch installation

In order to speed up the singular value decomposition, we highly recommend the installation of `pytorch`.
Users can instead also opt to use the slower R native SVD. For this, please set the argument `python = FALSE` wherever applicable in the package vignette.

### Install pytorch with reticulate

    library(reticulate)
    install_miniconda() 
    conda_install(envname = "r-reticulate", packages = "numpy")
    conda_install(envname = "r-reticulate", packages = "pytorch")

### Manually install pytorch with conda

Download the appropriate Miniconda installer for your system from [the conda website](https://docs.conda.io/en/latest/miniconda.html). 
Follow the installation instructions on their website and make sure the R package `reticulate` is also installed before proceeding.
Once installed, list all available conda environments via <br>
`conda info --envs` <br>
One of the environments should have `r-reticulate` in its name. Depending on where
you installed it and your system, the exact path might be different.
Activate the environment and install pytorch into it.

    conda activate ~/.local/share/r-miniconda/envs/r-reticulate # change path accordingly.
    conda install numpy
    conda install pytorch


## Feature overview

Please run 
    
    vignette("APL")

after installation with `build_vignettes = TRUE` for an introduction into the package.
