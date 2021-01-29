

# APL


## Installation

As of now all Bioconductor dependencies ("SingleCellExperiment") have to be installed manually:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("SingleCellExperiment")

The package is then installed as follows:

    library(devtools)
    install_github("ClemensKohl/APL")

Additionally, python with the libraries torch and numpy have to be installed if you want to use the pytorch SVD implementation.
One way to do this is with the R package reticulate:

    # TODO: DOES NOT WORK CURRENTLY DO NOT RUN! 
    # cant import/find pytorch for some reason after install.
    
    library(reticulate)
    install_miniconda() # Installs Miniconda
    conda_create("APLpy") # Creates a new environment
    conda_install("APLpy", c("numpy", "pytorch")) # Installs numpy and pytorch to the newly created env.

    # Call before running the code below
    use_condaenv("APLpy")

## Feature overview

The package contains functions for correspondence analysis, to display the results in biplots and to explore the dataset with Association Plots and to identify genes of interest.
The starting input can either be a count matrix or alternatively a Seurat or SingleCellExperiment container.

### Basic functionality

Here we demonstrate the basic functionality on a count matrix which consists of a subset of the GTEx data:

    library(APL)
    library(dplyr)
    library(readr)
    
    set.seed(1234)
    
    # GTEx dataset included in package.
    gtex[1:5,1:5]
    
    # Filter out genes without reads
    gtex <- gtex[rowSums(gtex) > 0,]
  

To run Correspondence Analysis we can use the `cacomp()` function. By default it will also calculate the standard coordinates of the rows (genes) and columns (samples/cells), as well as the principal coordinates of the rows. This will be of use later for plotting, but can be turned off if undesired. Here we additionally only keep the 5000 most variable genes as determined by the  variance of the chisquare components matrix in order to filter out genes that do not vary over the samples (here tissues). Currently there are two functions available to perform singular value decomposition: the base R `svd` function and SVD as implemented in pytorch in python. The latter significantly speeds up the computation and it is therefore highly recommended to use it whenever possible (`python = TRUE`), however a working python installation with numpy and torch installed is required.

    # Change the path to your python installation
    reticulate::use_python("/usr/bin/python3", required = TRUE)
    
    ca <- cacomp(obj = gtex,
                 top = 5000,
                 python = TRUE)
                 
    names(ca)


A convenient way to estimate the number of dimensions that we should keep for downstream analysis is to use `pick_dims()`. Several different methods are implemented to estimate the number of interesting dimensions, but here we use a formalization of the elbow rule to find a cutoff. This relies on rerunning `cacomp` on permutations of the data, which can take long for large count matrices, in which case either the number of `reps` should be decreased or one of the other methods can be chosen.

    pd <- pick_dims(obj = ca,
                    mat = gtex,
                    method = "elbow_rule",
                    reps = 3,
                    return_plot = TRUE,
                    python = TRUE)
    
    pd$dims

    
This gives us 29 dimensions to keep. The accompanying scree plot for visual inspection of the results can also be plotted:

    pd$plot
    
We can further explore the data in an assymetric biplot of the first 3 dimensions of the CA results to give us an idea of which samples might be interesting.

    ca_3Dplot(obj = ca)
    
In the interactive 3D plot we can see that samples derived from the Pancreas are quite different from other tissues. In order to further explore which genes are most highly associated with these samples we can use `runAPL()`. This function will calculate the coordinates of the genes and samples in the Association Plot. The further right a gene is located the more highly is it associated with the group of chosen samples. The x-Axis can be interpreted as the direction of the centroid of the chosen samples. The further up a gene is located on the y-Axis, the more it is also associated with other samples, which "pull" on the gene. Therefore a highly condition specific gene would be found far too the right with a very low y-value. By setting `score = TRUE` we can rank the genes on their association with the chosen sample group while also considering that we are not interested in genes that are too highly associated with other conditions too. By default the top 10 most highly ranked genes will be displayed. To mark specific genes, use the `mark_rows` parameter. By also submitting the previously computed `ca` we can save time, as `runAPL` would otherwise have to recompute the values.

    runAPL(obj = gtex,
           caobj = ca,
           dims = 29, 
           group = grep("Pancreas", x = colnames(gtex)),
           score = TRUE) 
           

The wrapper `runAPL` makes it easy to try out different groups and explore the data quickly. If, however, we want a more fine grained control we can also call all functions seperately:

    ca <- cacomp(obj = gtex,
                     coords = TRUE,
                     princ_coords = 1,
                     dims = 29,
                     top = 5000,
                     python = TRUE)
    
    group <-  grep("Pituitary", x = colnames(gtex))
    ca <- apl_coords(caobj = ca, 
                     group = group)
    ca <- apl_score(caobj = ca,
                    mat = gtex,
                    dims = ca$dims,
                    group = ca$group,
                    reps = 3)
    apl(caobj = ca,
        type = "plotly",
        rowlabels = TRUE,
        collabels = TRUE,
        rows_idx = head(ca$APL_score$Row_num,10),
        cols_idx = ca$group)
        
If we do not want to calculate the CA coordinates and subset the dimensions in the beginnings (e.g. because we want to run `pick_dims` first) we can do this later by running:

    ca <- ca_coords(caobj = ca,
                    dims = 29,
                    princ_coords = 1)

A biplot of the data can be generated similarly to the 3D plot by running:

    ca_biplot(obj = ca)
 

### Seurat/SingleCellExperiment integration

`cacomp`, `pick_dims`, `runAPL` and the biplots work equally well with Seurat and SingleCellExperiment containers:

#### Seurat

    library(Seurat)
    
    pbmc_small # included in Seurat package
    
    pbmc_small <- cacomp(obj = pbmc_small,
                         assay = "RNA",
                         coords = TRUE,
                         return_input = TRUE)
                         
    pd <- pick_dims(pbmc_small,
                    assay = "RNA",
                    method = "elbow_rule",
                    reps = 3,
                    return_plot = FALSE)
    
    runAPL(pbmc_small,
           assay = "RNA",
           dims = pd,
           group = grep("g1", pbmc_small@meta.data$groups),
           score = TRUE,
           nrow = 10,
           reps = 5,
           python = TRUE)

#### SingleCellExperiment

    library(scRNAseq)
    library(scater)
    
    # Load data
    sce <- ReprocessedAllenData("tophat_counts")
    
    counts <- assay(sce, "tophat_counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    logcounts(sce) <- log2(t(t(counts)/size.factors) + 1)
    counts(sce) <- assay(sce, "tophat_counts")
    
    sce <- cacomp(obj = sce, assay = "counts", dims = 50, top = 5000, return_input = TRUE)
    pd <- pick_dims(sce, assay = "counts", method = "elbow_rule", return_plot = FALSE)

    grp <- grep("L4 Arf5", colData(sce)$Primary.Type)
    
    runAPL(obj = sce,
            group = grp,
            dims = pd,
            assay = "counts",
            score = TRUE,
            reps = 3,
            nrow = 10)

## TODO

- S3 Constructor function! https://adv-r.hadley.nz/s3.html#s3-constructor
- document datasets
- Different (?) Documentation for S3 methods.
- Implement print.cacomp() method for nicer printing of cacomp objects.

