

# APL

When working with `APL` package please cite:
```
Association Plots: Visualizing associations in high-dimensional correspondence analysis biplots
Elzbieta Gralinska, Martin Vingron
bioRxiv 2020.10.23.352096; doi: https://doi.org/10.1101/2020.10.23.352096
```

## Installation

    library(devtools)
    install_github(" elagralinska/APL-Rpackage")
    
## Pytorch installation

In order to speed up the singular value decomposition, we highly recommend the installation of `pytorch`.
Users can instead also opt to use the slower R native SVD. For this, please turn the argument `python = FALSE` wherever applicable in this vignette.

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

after installation for an introduction into the package.
