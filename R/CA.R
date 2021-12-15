#' @include constructor.R
NULL

#' Compute Standardized Residuals
#'
#' @description
#' `comp_std_residuals` computes the standardized residual matrix S,
#' which is the basis for correspondence analysis and serves
#' as input for singular value decomposition (SVD).
#'
#' @details
#' Calculates standardized residual matrix S from the proportion matrix P and
#' the expected values E according to \eqn{S = \frac{(P-E)}{sqrt(E)}}.
#'
#' @param mat A numerical matrix or coercible to one by `as.matrix()`.
#' Should have row and column names.
#' @return
#' A named list with standardized residual matrix "S",
#' grand total of the original matrix "tot"
#' as well as row and column masses "rowm" and "colm" respectively.
#'
comp_std_residuals <- function(mat){

  if (!is(mat, "matrix")){
    mat <- as.matrix(mat)
  }
  stopifnot(
    "Input matrix does not have any rownames!" = !is.null(rownames(mat)))
  stopifnot(
    "Input matrix does not have any colnames!" = !is.null(colnames(mat)))

  tot <- sum(mat)
  P <- mat/tot               # proportions matrix
  rowm <- rowSums(P)          # row masses
  colm <- colSums(P)          # column masses

  E <- rowm %o% colm      # expected proportions
  S <-  (P - E) / sqrt(E)         # standardized residuals
  S[is.nan(S)] <- 0

  out <- list("S"=S, "tot"=tot, "rowm"=rowm, "colm"=colm)
  return(out)
}

#' removes 0-only rows and columns in a matrix.
#'
#' @param obj A matrix.
#' @return Input matrix with rows & columns consisting of only 0 removed.
rm_zeros <- function(obj){
  stopifnot(is(obj, "matrix"))

  no_zeros_rows <- rowSums(obj) > 0
  no_zeros_cols <- colSums(obj) > 0
  if (sum(!no_zeros_rows) != 0){
    ## Delete genes with only zero values across all columns
    warning("Matrix contains rows with only 0s. ",
            "These rows were removed. ",
            "If undesired set rm_zeros = FALSE.")
    obj <- obj[no_zeros_rows,]
  }
  if (sum(!no_zeros_cols) != 0){
    ## Delete cells with only zero values across all genes
    warning("Matrix contains columns with only 0s. ",
            "These columns were removed. ",
            "If undesired set rm_zeros = FALSE.")
    obj <- obj[,no_zeros_cols]
  }

  return(obj)
}
#' Find most variable rows
#'
#' @description
#' Calculates the variance of the chi-square component matrix and selects the 
#' rows with the highest variance, e.g. 5,000.
#'
#' @return
#' Returns a matrix, which consists of the top variable rows of mat.
#'
#' @param mat A numeric matrix. For sequencing a count matrix,
#' gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param top Integer. Number of most variable rows to retain. Default 5000.
#' @export
#' @examples
#' set.seed(1234)
#'
#' # Simulate counts
#'cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'               x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Choose top 5000 most variable genes
#' cnts <- var_rows(mat = cnts, top = 5000)
#'
#'
var_rows <- function(mat, top = 5000){

  res <-  comp_std_residuals(mat=mat)

  if(top>nrow(mat)) {
    warning("Top is larger than the number of rows in matrix. ",
            "Top was set to nrow(mat).")
  }
  
  top <- min(nrow(mat), top)

  chisquare <- res$tot * (res$S^2)		# chi-square components matrix
  variances <- apply(chisquare,1,var) #row-wise variances
  ix_var <- order(-variances)
  mat <- mat[ix_var[seq_len(top)],] # choose top rows
  return(mat)

}


#' Internal function for `cacomp`
#'
#' @description
#' `run_cacomp` performs correspondence analysis on a matrix and returns the 
#' transformed data.
#'
#' @details
#' The calculation is performed according to the work of Michael Greenacre. 
#' Singular value decomposition
#' can be performed either with the base R function `svd` or preferably by the 
#' faster
#' pytorch implementation (python = TRUE). When working with large matrices, 
#' CA coordinates and
#' principal coordinates should only be computed when needed to save 
#' computational time.
#'
#' @return
#' Returns a named list of class "cacomp" with components
#' U, V and D: The results from the SVD.
#' row_masses and col_masses: Row and columns masses.
#' top_rows: How many of the most variable rows/genes were retained for the 
#' analysis.
#' tot_inertia, row_inertia and col_inertia: Only if inertia = TRUE. Total, 
#' row and column inertia respectively.
#' @references
#' Greenacre, M. Correspondence Analysis in Practice, Third Edition, 2017.

#' @param obj A numeric matrix or Seurat/SingleCellExperiment object. For 
#' sequencing a count matrix, gene expression values with genes in rows and 
#' samples/cells in columns.
#' Should contain row and column names.
#' @param coords Logical. Indicates whether CA standard coordinates should be 
#' calculated. Default TRUE
#' @param python A logical value indicating whether to use singular-value 
#' decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @param princ_coords Integer. Number indicating whether principal 
#' coordinates should be calculated for the rows (=1), columns (=2), 
#' both (=3) or none (=0).
#' Default 1.
#' @param dims Integer. Number of CA dimensions to retain. Default NULL 
#' (keeps all dimensions).
#' @param top Integer. Number of most variable rows to retain. Default 5000.
#' @param inertia Logical.. Whether total, row and column inertias should be 
#' calculated and returned. Default TRUE.
#' @param rm_zeros Logical. Whether rows & cols containing only 0s should be 
#' removed. Keeping zero only rows/cols might lead to unexpected results. 
#' Default TRUE.
#' @param ... Arguments forwarded to methods.
run_cacomp <- function(obj,
                   coords = TRUE,
                   princ_coords = 3,
                   python = FALSE,
                   dims = NULL,
                   top = 5000,
                   inertia = TRUE,
                   rm_zeros = TRUE,
                   ...){

  stopifnot("Input matrix does not have any rownames!" =
              !is.null(rownames(obj)))
  stopifnot("Input matrix does not have any colnames!" = 
              !is.null(colnames(obj)))

  if (rm_zeros == TRUE){
    obj <- rm_zeros(obj)
  }


  # Choose only top # of variable genes
  if (is.null(top) || top == nrow(obj)) {
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  } else if (!is.null(top) && top < nrow(obj)){

    obj <- var_rows(mat = obj, top = top)
    res <-  comp_std_residuals(mat=obj)
    toptmp <- top
    
  } else if (top > nrow(obj)) {
    warning("\nParameter top is >nrow(obj) and therefore ignored.")
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  } else {
    warning("\nUnusual input for top, argument ignored.")
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  }

  S <- res$S
  tot <- res$tot
  rowm <- res$rowm
  colm <- res$colm
  rm(res)

  k <- min(dim(S))-1

  if (is.null(dims)) dims <- k
  if (dims > k) dims <- k
  # S <- (diag(1/sqrt(r)))%*%(P-r%*%t(c))%*%(diag(1/sqrt(c)))
  # message("Running singular value decomposition ...")

  if (python == TRUE){
    svd_torch <- NULL
    # require(reticulate)
    # source_python('./python_svd.py')
    reticulate::source_python(system.file("python/python_svd.py", package = "APL"))
    SVD <- svd_torch(S)
    # SVD <- svd_linalg_torch(S)
    names(SVD) <- c("U", "D", "V")
    SVD$D <- as.vector(SVD$D)

  } else {

    SVD <- svd(S, nu = dims, nv = dims)
    names(SVD) <- c("D", "U", "V")
    SVD <- SVD[c(2, 1, 3)]
    if(length(SVD$D) > dims) SVD$D <- SVD$D[seq_len(dims)]
  }

  names(SVD$D) <- paste0("Dim", seq_len(length(SVD$D)))
  dimnames(SVD$V) <- list(colnames(S), paste0("Dim", seq_len(ncol(SVD$V))))
  dimnames(SVD$U) <- list(rownames(S), paste0("Dim", seq_len(ncol(SVD$U))))


  if(inertia == TRUE){
    #calculate inertia
    SVD$tot_inertia <- sum(SVD$D^2)
    SVD$row_inertia <- rowSums(S^2)
    SVD$col_inertia <- colSums(S^2)
  }

  SVD$row_masses <- rowm
  SVD$col_masses <- colm
  SVD$top_rows <- toptmp

  SVD <- do.call(new_cacomp, SVD)
  SVD <- subset_dims(SVD, dims)
  # class(SVD) <- "cacomp"

  if (coords == TRUE){
    # message("Calculating coordinates...")

    SVD <- ca_coords(caobj = SVD,
                     dims = dims,
                     princ_coords = princ_coords,
                     princ_only = FALSE)
  } else {
    if(!is.null(dims)){
      if (dims >= length(SVD@D)){
        if (dims > length(SVD@D)){
          warning("Chosen number of dimensions is larger than the ",
                  "number of dimensions obtained from the singular ",
                  "value decomposition. Argument ignored.")
        }
        SVD@dims <- length(SVD@D)
      } else {
        dims <- min(dims, length(SVD@D))
        SVD@dims <- dims

        dims <- seq(dims)

        # subset to number of dimensions
        SVD@U <- SVD@U[,dims]
        SVD@V <- SVD@V[,dims]
        SVD@D <- SVD@D[dims]
      }
    } else {
      SVD@dims <- length(SVD@D)
    }
  }

  stopifnot(validObject(SVD))
  return(SVD)
}


#' Correspondance Analysis
#'
#' @description
#' `cacomp` performs correspondence analysis on a matrix or
#' Seurat/SingleCellExperiment object and returns the transformed data.
#'
#' @details
#' The calculation is performed according to the work of Michael Greenacre. 
#' Singular value decomposition can be performed either with the base R 
#' function `svd` or preferably by the faster pytorch implementation 
#' (python = TRUE). When working with large matrices, CA coordinates and
#' principal coordinates should only be computed when needed to save 
#' computational time.
#'
#' @return
#' Returns a named list of class "cacomp" with components
#' U, V and D: The results from the SVD.
#' row_masses and col_masses: Row and columns masses.
#' top_rows: How many of the most variable rows were retained for the analysis.
#' tot_inertia, row_inertia and col_inertia: Only if inertia = TRUE.
#' Total, row and column inertia respectively.
#' @references
#' Greenacre, M. Correspondence Analysis in Practice, Third Edition, 2017.

#' @param obj A numeric matrix or Seurat/SingleCellExperiment object.
#' For sequencing a count matrix, gene expression values with genes in rows 
#' and samples/cells in columns.
#' Should contain row and column names.
#' @param coords Logical. Indicates whether CA standard coordinates should be 
#' calculated. Default TRUE
#' @param python A logical value indicating whether to use singular-value 
#' decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @param princ_coords Integer. Number indicating whether principal 
#' coordinates should be calculated for the rows (=1), columns (=2), 
#' both (=3) or none (=0).
#' Default 1.
#' @param dims Integer. Number of CA dimensions to retain. Default NULL 
#' (keeps all dimensions).
#' @param top Integer. Number of most variable rows to retain. 
#' Default 5000. (set NULL to keep all).
#' @param inertia Logical.. Whether total, row and column inertias should be 
#' calculated and returned. Default TRUE.
#' @param rm_zeros Logical. Whether rows & cols containing only 0s should be 
#' removed.
#' Keeping zero only rows/cols might lead to unexpected results. Default TRUE.
#' @param ... Arguments forwarded to methods.
#' @examples
#' # Simulate scRNAseq data.
#' cnts <- data.frame(cell_1 = rpois(10, 5),
#'                    cell_2 = rpois(10, 10),
#'                    cell_3 = rpois(10, 20))
#' rownames(cnts) <- paste0("gene_", 1:10)
#' cnts <- as.matrix(cnts)
#'
#' # Run correspondence analysis.
#' ca <- cacomp(obj = cnts, princ_coords = 3, top = 5)
#' @export
setGeneric("cacomp", function(obj,
                              coords = TRUE,
                              princ_coords = 3,
                              python = FALSE,
                              dims = NULL,
                              top = 5000,
                              inertia = TRUE,
                              rm_zeros = TRUE,
                              ...) {
  standardGeneric("cacomp")
})


#' @rdname cacomp
#' @export
setMethod(f = "cacomp",
          signature=(obj="matrix"),
          function(obj,
                   coords = TRUE,
                   princ_coords = 3,
                   python = FALSE,
                   dims = NULL,
                   top = 5000,
                   inertia = TRUE,
                   rm_zeros = TRUE,
                   ...){

    caobj <- run_cacomp(obj = obj,
                        coords = coords,
                        princ_coords = princ_coords,
                        python = python,
                        dims = dims,
                        top = top,
                        inertia = inertia,
                        rm_zeros = rm_zeros,
                        ...)

    return(caobj)

})


#' Correspondance Analysis for Seurat objects
#'
#' @description
#' `cacomp.seurat` performs correspondence analysis on an assay from a Seurat 
#' container and stores the standardized coordinates of the columns (= cells) 
#' and the principal coordinates of the rows (= genes) as a DimReduc Object in 
#' the Seurat container.
#'
#' @return
#' If return_imput = TRUE with Seurat container: Returns input obj of class 
#' "Seurat" with a new Dimensional Reduction Object named "CA".
#' Standard coordinates of the cells are saved as embeddings,
#' the principal coordinates of the genes as loadings and
#' the singular values (= square root of principal intertias/eigenvalues)
#' are stored as stdev.
#' To recompute a regular "cacomp" object without rerunning cacomp use 
#' `as.cacomp()`.
#' @param assay Character. The assay from which extract the count matrix for 
#' SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#' @param slot character. The slot of the Seurat assay. Default "counts".
#' @param return_input Logical. If TRUE returns the input 
#' (SingleCellExperiment/Seurat object) with the CA results saved in the 
#' reducedDim/DimReduc slot "CA".
#' Otherwise returns a "cacomp". Default FALSE.
#' @param ... Other parameters
#' @rdname cacomp
#' @export
#' @examples
#' 
#' ###########
#' # Seurat  #
#' ###########
#' library(Seurat)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                      x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Create Seurat object
#' seu <- CreateSeuratObject(counts = cnts)
#'
#' # Run CA and save in dim. reduction slot
#' seu <- cacomp(seu, return_input = TRUE, assay = "RNA", slot = "counts")
#'
#' # Run CA and return cacomp object
#' ca <- cacomp(seu, return_input = FALSE, assay = "RNA", slot = "counts")
setMethod(f = "cacomp",
          signature=(obj="Seurat"),
          function(obj,
                   coords = TRUE,
                   princ_coords = 3,
                   python = FALSE,
                   dims = NULL,
                   top = 5000,
                   inertia = TRUE,
                   rm_zeros = TRUE,
                   ...,
                   assay = Seurat::DefaultAssay(obj),
                   slot = "counts",
                   return_input = FALSE){

  stopifnot("obj doesnt belong to class 'Seurat'" = is(obj, "Seurat"))

  stopifnot("Set coords = TRUE when inputting a Seurat object and return_input = TRUE." = coords == TRUE)


  seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = slot)
  seu <- as.matrix(seu)

  caobj <- run_cacomp(obj = seu,
                      coords = coords,
                      top = top,
                      princ_coords = princ_coords,
                      dims = dims,
                      python = python,
                      rm_zeros = rm_zeros,
                      inertia = inertia,
                      ...)

  if (return_input == TRUE){
    colnames(caobj@V) <- paste0("DIM_", seq(ncol(caobj@V)))
    colnames(caobj@U) <- paste0("DIM_", seq(ncol(caobj@U)))

    obj[["CA"]] <- Seurat::CreateDimReducObject(embeddings = caobj@std_coords_cols,
                                               loadings = caobj@prin_coords_rows,
                                               stdev = caobj@D,
                                               key = "Dim_",
                                               assay = assay)

    return(obj)
  } else {
    return(caobj)
  }

})


#' Correspondance Analysis for SingleCellExperiment objects
#'
#' @description
#' `cacomp.SingleCellExperiment` performs correspondence analysis on an assay 
#' from a SingleCellExperiment and stores the standardized coordinates
#'  of the columns (= cells) and the principal coordinates of the rows 
#'  (= genes) as a matrix in the SingleCellExperiment container.
#'
#' @return
#' If return_input =TRUE for SingleCellExperiment input returns a 
#' SingleCellExperiment object with a matrix of standardized coordinates of 
#' the columns in
#' reducedDim(obj, "CA"). Additionally, the matrix contains the following 
#' attributes:
#' "prin_coords_rows": Principal coordinates of the rows.
#' "singval": Singular values. For the explained inertia of each principal 
#' axis calculate singval^2.
#' "percInertia": Percent explained inertia of each principal axis.
#' To recompute a regular "cacomp" object from a SingleCellExperiment without 
#' rerunning cacomp use `as.cacomp()`.
#' @param assay Character. The assay from which extract the count matrix for 
#' SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#' @param return_input Logical. If TRUE returns the input 
#' (SingleCellExperiment/Seurat object) with the CA results saved in the 
#' reducedDim/DimReduc slot "CA".
#'  Otherwise returns a "cacomp". Default FALSE.
#' @rdname cacomp
#' @export
#' @examples
#' 
#' ########################
#' # SingleCellExperiment #
#' ########################
#' library(SingleCellExperiment)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#' logcnts <- log2(cnts + 1)
#'
#' # Create SingleCellExperiment object
#' sce <- SingleCellExperiment(assays=list(counts=cnts, logcounts=logcnts))
#'
#' # run CA and save in dim. reduction slot.
#' sce <- cacomp(sce, return_input = TRUE, assay = "counts") # on counts
#' sce <- cacomp(sce, return_input = TRUE, assay = "logcounts") # on logcounts
#'
#' # run CA and return cacomp object.
#' ca <- cacomp(sce, return_input = FALSE, assay = "counts")
setMethod(f = "cacomp",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   coords = TRUE,
                   princ_coords = 3,
                   python = FALSE,
                   dims = NULL,
                   top = 5000,
                   inertia = TRUE,
                   rm_zeros = TRUE,
                   ...,
                   assay = "counts",
                   return_input = FALSE){

  stopifnot("obj doesnt belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))
  stopifnot("Set coords = TRUE when inputting a SingleCellExperiment object and return_input = TRUE." = coords == TRUE)

  mat <- SummarizedExperiment::assay(obj, assay)
  mat <- as.matrix(mat)

  top <- min(nrow(mat), top)

  caobj <- run_cacomp(obj = mat,
                     coords = coords,
                     top = top,
                     princ_coords = princ_coords,
                     dims = dims,
                     python = python,
                     rm_zeros = rm_zeros,
                     inertia = inertia,
                     ...)

  if (return_input == TRUE){
    prinInertia <- caobj@D^2
    percentInertia <- prinInertia / sum(prinInertia) * 100

    # Saving the results
    ca <- caobj@std_coords_cols
    attr(ca, "prin_coords_rows") <- caobj@prin_coords_rows
    attr(ca, "singval") <- caobj@D
    attr(ca, "percInertia") <- percentInertia

    SingleCellExperiment::reducedDim(obj, "CA") <- ca

    return(obj)

  } else {
    return(caobj)
  }

})


#' Subset dimensions of a caobj
#'
#' @description Subsets the dimensions according to user input.
#'
#' @return Returns caobj.
#'
#' @param caobj A caobj.
#' @param dims Integer. Number of dimensions.
#' @examples
#' # Simulate scRNAseq data.
#' cnts <- data.frame(cell_1 = rpois(10, 5),
#'                    cell_2 = rpois(10, 10),
#'                    cell_3 = rpois(10, 20))
#' rownames(cnts) <- paste0("gene_", 1:10)
#' cnts <- as.matrix(cnts)
#'
#' # Run correspondence analysis.
#' ca <- cacomp(cnts)
#' ca <- subset_dims(ca, 2)
#' @export
subset_dims <- function(caobj, dims){
  
  if (dims == 1){stop("Please choose more than 1 dimension.")}
  
  stopifnot(is(caobj, "cacomp"))
  
  if (is.null(dims)){
    return(caobj)
  }
  
  if(dims > length(caobj@D)){
    warning("dims is larger than the number of available dimensions.",
            " Argument ignored")
  } else if (dims == length(caobj@D)){
    caobj@dims <- dims
    return(caobj)
  }

  dims <- min(dims, length(caobj@D))
  caobj@dims <- dims
  dims <- seq(dims)
  caobj@U <- caobj@U[,dims]
  caobj@V <- caobj@V[,dims]
  caobj@D <- caobj@D[dims]

  if (!is.empty(caobj@std_coords_cols)){
    caobj@std_coords_cols <- caobj@std_coords_cols[,dims]
  }
  if (!is.empty(caobj@prin_coords_cols)){
    caobj@prin_coords_cols <- caobj@prin_coords_cols[,dims]
  }

  if (!is.empty(caobj@std_coords_rows)){
    caobj@std_coords_rows <- caobj@std_coords_rows[,dims]
  }
  if (!is.empty(caobj@prin_coords_rows)){
    caobj@prin_coords_rows <- caobj@prin_coords_rows[,dims]
  }

  stopifnot(validObject(caobj))
  return(caobj)
}


#' Calculate correspondence analysis row and column coordinates.
#'
#' @description `ca_coords` calculates the standardized and principal 
#' coordinates of the rows and columns in CA space.
#'
#' @details
#' Takes a "cacomp" object and calculates standardized and principal 
#' coordinates for the visualization of CA results in a biplot or
#' to subsequently calculate coordinates in an Association Plot.
#'
#' @return
#' Returns input object with coordinates added.
#' std_coords_rows/std_coords_cols: Standardized coordinates of rows/columns.
#' prin_coords_rows/prin_coords_cols: Principal coordinates of rows/columns.
#'
#' @param caobj A "cacomp" object as outputted from `cacomp()`.
#' @param dims Integer indicating the number of dimensions to use for the 
#' calculation of coordinates.
#' All elements of caobj (where applicable) will be reduced to the given 
#' number of dimensions. Default NULL (keeps all dimensions).
#' @param princ_only Logical, whether only principal coordinates should be 
#' calculated.
#' Or, in other words, whether the standardized coordinates are already 
#' calculated and stored in `caobj`. Default `FALSE`.
#' @param princ_coords Integer. Number indicating whether principal 
#' coordinates should be calculated for the rows (=1), columns (=2), both (=3) 
#' or none (=0).
#' Default 3.
#' @examples
#' # Simulate scRNAseq data.
#' cnts <- data.frame(cell_1 = rpois(10, 5),
#'                    cell_2 = rpois(10, 10),
#'                    cell_3 = rpois(10, 20))
#' rownames(cnts) <- paste0("gene_", 1:10)
#' cnts <- as.matrix(cnts)
#'
#' # Run correspondence analysis.
#' ca <- cacomp(obj = cnts, princ_coords = 1)
#' ca <- ca_coords(ca, princ_coords = 3)
#' @export
ca_coords <- function(caobj, dims=NULL, princ_coords = 3, princ_only = FALSE){

  stopifnot(is(caobj, "cacomp"))
  stopifnot(dims <= length(caobj@D))

  if(!is.null(dims)){
    if (dims > length(caobj@D)){
      warning("Chosen dimensions are larger than the number of ",
              "dimensions obtained from the singular value ",
              "decomposition. Argument ignored.")
     }
      caobj <- subset_dims(caobj = caobj, dims = dims)
    }


  if(princ_only == FALSE){

    #standard coordinates
    if(dims == 1 && !is.null(dims)){
      caobj@std_coords_rows <- caobj@U/sqrt(caobj@row_masses)
      caobj@std_coords_cols <- caobj@V/sqrt(caobj@col_masses)
    } else {
      caobj@std_coords_rows <- sweep(x = caobj@U,
                                     MARGIN = 1,
                                     STATS = sqrt(caobj@row_masses),
                                     FUN = "/")
      caobj@std_coords_cols <- sweep(x = caobj@V,
                                     MARGIN = 1,
                                     STATS = sqrt(caobj@col_masses),
                                     FUN = "/")
    }


    # Ensure no NA/Inf after dividing by 0.
    caobj@std_coords_rows[is.na(caobj@std_coords_rows)] <- 0
    caobj@std_coords_cols[is.na(caobj@std_coords_cols)] <- 0
    caobj@std_coords_rows[is.infinite(caobj@std_coords_rows)] <- 0
    caobj@std_coords_cols[is.infinite(caobj@std_coords_cols)] <- 0

  }


  stopifnot("princ_coords must be either 0, 1, 2 or 3" = 
              (princ_coords == 0 || 
               princ_coords == 1 || 
               princ_coords == 2 || 
               princ_coords == 3))

  if(princ_coords != 0){
    stopifnot(!is.empty(caobj@std_coords_rows))
    stopifnot(!is.empty(caobj@std_coords_cols))

      if (princ_coords == 1){
        #principal coordinates for rows
        if (dims == 1 && !is.null(dims)){
          caobj@prin_coords_rows <- caobj@std_coords_rows*caobj@D
        } else {
          caobj@prin_coords_rows <- sweep(caobj@std_coords_rows,
                                          2,
                                          caobj@D,
                                          "*")
        }

      } else if (princ_coords == 2) {
        #principal coordinates for columns
        if (dims == 1 && !is.null(dims)){
          caobj@prin_coords_cols <- caobj@std_coords_cols*caobj@D
        } else {
          caobj@prin_coords_cols <- sweep(caobj@std_coords_cols,
                                          2,
                                          caobj@D,
                                          "*")
        }
      } else if (princ_coords  == 3) {

        if (dims == 1 && !is.null(dims)){
          #principal coordinates for rows
          caobj@prin_coords_rows <- caobj@std_coords_rows*caobj@D
          #principal coordinates for columns
          caobj@prin_coords_cols <- caobj@std_coords_cols*caobj@D
        } else {
          #principal coordinates for rows
          caobj@prin_coords_rows <- sweep(caobj@std_coords_rows,
                                          2,
                                          caobj@D,
                                          "*")
          #principal coordinates for columns
          caobj@prin_coords_cols <- sweep(caobj@std_coords_cols,
                                          2,
                                          caobj@D,
                                          "*")
        }

      }

  }

  stopifnot(validObject(caobj))
  return(caobj)
}


#' Scree Plot
#'
#'@description Plots a scree plot.
#'
#'@return
#'Returns a ggplot object.
#'
#'@param df A data frame with columns "dims" and "inertia".
scree_plot <- function(df){

  stopifnot(c("dims", "inertia") %in% colnames(df))

  avg_inertia <- 100/nrow(df)
  max_num_dims <- nrow(df)

  screeplot <- ggplot2::ggplot(df, ggplot2::aes(x=.data$dims,
                                                y=.data$inertia)) +
    ggplot2::geom_col(fill="#4169E1") +
    ggplot2::geom_line(color="#B22222", size=1) +
    ggplot2::labs(
         title = "Scree plot of explained inertia per dimensions and the average inertia",
         y="Explained inertia [%]",
         x="Dimension") +
    ggplot2::theme_bw()
  return(screeplot)
}

#' Runs elbow method
#' 
#' @description Helper function for pick_dims() to run the elbow method.
#' 
#' @param obj A "cacomp" object as outputted from `cacomp()`
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression 
#' values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param reps Integer. Number of permutations to perform when choosing 
#' "elbow_rule".
#' @param return_plot TRUE/FALSE. Whether a plot should be returned when 
#' choosing "elbow_rule".
#' @param python A logical value indicating whether to use singular value 
#' decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @return
#' `elbow_method` (for `return_plot=TRUE`) returns a list with two elements: 
#' "dims" contains the number of dimensions and "plot" a ggplot. if 
#' `return_plot=TRUE` it just returns the number of picked dimensions.
#' @references 
#' Ciampi, Antonio, González Marcos, Ana and Castejón Limas, Manuel. \cr
#' Correspondence analysis and 2-way clustering. (2005), SORT 29(1).
#' 
#' @examples 
#' 
#' # Get example data from Seurat
#' library(Seurat)
#' set.seed(2358)
#' cnts <- as.matrix(Seurat::GetAssayData(pbmc_small, "data"))
#' # Run correspondence analysis.
#' ca <- cacomp(obj = cnts)
#' 
#' # pick dimensions with the elbow rule. Returns list.
#' pd <- pick_dims(obj = ca,
#'                 mat = cnts,
#'                 method = "elbow_rule",
#'                 return_plot = TRUE,
#'                 reps = 10)
#' pd$plot
#' ca_sub <- subset_dims(ca, dims = pd$dims)
#' 
elbow_method <- function(obj,
                         mat,
                         reps,
                         python,
                         return_plot){
  ev <- obj@D^2
  expl_inertia <- (ev/sum(ev)) *100
  max_num_dims <- length(obj@D)
  
  matrix_expl_inertia_perm <- matrix(0, nrow = max_num_dims , ncol = reps)
  
  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  
  for (k in seq(reps)) {
    
    mat <- as.matrix(mat)
    mat_perm <- apply(mat, 2, FUN=sample)
    colnames(mat_perm) <- colnames(mat)
    rownames(mat_perm) <- seq_len(nrow(mat_perm))
    
    obj_perm <- cacomp(obj=mat_perm,
                       top = obj@top_rows,
                       dims = obj@dims,
                       coords = FALSE,
                       python = python)
    
    ev_perm <- obj_perm@D^2
    expl_inertia_perm <- (ev_perm/sum(ev_perm))*100
    
    matrix_expl_inertia_perm[,k] <- expl_inertia_perm
    colnames(matrix_expl_inertia_perm) <- paste0("perm",seq_len(reps))
    
    setTxtProgressBar(pb, k)
    
  }
  close(pb)
  
  
  if (return_plot == TRUE){
    df <- data.frame(dims = seq_len(max_num_dims),
                     inertia = expl_inertia)
    
    df <- cbind(df, matrix_expl_inertia_perm)
    
    screeplot <- scree_plot(df)
    
    for (k in seq_len(reps)) {
      
      colnm <- ggplot2::sym(paste0("perm",k))
      
      screeplot <- screeplot +
        ggplot2::geom_line(data = df, ggplot2::aes(x=.data$dims,
                                                   y=!!colnm),
                           color="black",
                           alpha=0.8,
                           linetype=2)
      
    }
  }
  
  avg_inertia_perm <- rowMeans(matrix_expl_inertia_perm)
  
  tmp <- as.integer(expl_inertia>avg_inertia_perm)
  if (sum(tmp)==0 || sum(tmp)==max_num_dims){
    dim_number <- max_num_dims								
  } else if (tmp[1] == 0){
    stop("Average inertia of the permutated data is above ",
         "the explained inertia of the data in the first dimension. ",
         "Please either try more permutations or a different method.")
  }else{
    dim_number <- length(tmp[cumsum(tmp == 0)<1 & tmp!=0])		
  }
  
  if (return_plot == FALSE){
    return(dim_number)
  } else {
    return(list("dims" = dim_number, "plot" = screeplot))
  }
}


#' Compute statistics to help choose the number of dimensions
#'
#' @description
#' Allow the user to choose from 4 different methods ("avg_inertia", 
#' "maj_inertia", "scree_plot" and "elbow_rule")
#' to estimate the number of dimensions that best represent the data.
#'
#' @details
#' * "avg_inertia" calculates the number of dimensions in which the inertia is 
#' above the average inertia.
#' * "maj_inertia" calculates the number of dimensions in which cumulatively 
#' explain up to 80% of the total inertia.
#' * "scree_plot" plots a scree plot.
#' * "elbow_rule" formalization of the commonly used elbow rule. Permutes the 
#' rows for each column and reruns `cacomp()` for a total of `reps` times.
#' The number of relevant dimensions is obtained from the point where the 
#' line for the explained inertia of the permuted data intersects with the 
#' actual data.
#'
#' @return
#' For `avg_inertia`, `maj_inertia` and `elbow_rule` (when `return_plot=FALSE`) 
#' returns an integer, indicating the suggested number of dimensions to use.
#' * `scree_plot` returns a ggplot object.
#' * `elbow_rule` (for `return_plot=TRUE`) returns a list with two elements: 
#' "dims" contains the number of dimensions and "plot" a ggplot.
#'
#' @param obj A "cacomp" object as outputted from `cacomp()`,
#' a "Seurat" object with a "CA" DimReduc object stored,
#' or a "SingleCellExperiment" object with a "CA" dim. reduction stored.
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression 
#' values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param method String. Either "scree_plot", "avg_inertia", "maj_inertia" or 
#' "elbow_rule" (see Details section). Default "scree_plot".
#' @param reps Integer. Number of permutations to perform when choosing 
#' "elbow_rule". Default 3.
#' @param return_plot TRUE/FALSE. Whether a plot should be returned when 
#' choosing "elbow_rule". Default FALSE.
#' @param python A logical value indicating whether to use singular value 
#' decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @param ... Arguments forwarded to methods.
#' @examples
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis.
#' ca <- cacomp(obj = cnts)
#'
#' # pick dimensions with the elbow rule. Returns list.
#'
#' set.seed(2358)
#' pd <- pick_dims(obj = ca,
#'                 mat = cnts,
#'                 method = "elbow_rule",
#'                 return_plot = TRUE,
#'                 reps = 10)
#' pd$plot
#' ca_sub <- subset_dims(ca, dims = pd$dims)
#'
#' # pick dimensions which explain cumulatively >80% of total inertia.
#' # Returns vector.
#' pd <- pick_dims(obj = ca,
#'                 method = "maj_inertia")
#' ca_sub <- subset_dims(ca, dims = pd)
#' @export
#' @md
setGeneric("pick_dims", function(obj,
                                 mat = NULL,
                                 method="scree_plot",
                                 reps=3,
                                 python = TRUE,
                                 return_plot = FALSE,
                                 ...) {
  standardGeneric("pick_dims")
})


#' @rdname pick_dims
#' @export
setMethod(f = "pick_dims",
          signature=(obj="cacomp"),
          function(obj,
                   mat = NULL,
                   method="scree_plot",
                   reps=3,
                   python = TRUE,
                   return_plot = FALSE,
                   ...){

  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }

  ev <- obj@D^2
  expl_inertia <- (ev/sum(ev)) *100
  max_num_dims <- length(obj@D)

  if (method == "avg_inertia"){
    # Method 1: Dim's > average inertia
    # percentage of inertia explained by 1 dimension (on average)
    avg_inertia <- 100/max_num_dims
    # result: number of dimensions, all of which explain more than avg_inertia
    dim_num <- sum(expl_inertia > avg_inertia)  
    return(dim_num)

  } else if (method == "maj_inertia"){
    # Method 2: Sum of dim's > 80% of the total inertia
    # the first dimension for which the cumulative sum of inertia (from dim1 
    # up to given dimension) is higher than 80%
    dim_num <- min(which(cumsum(expl_inertia)>80)) 
    return(dim_num)

  } else if (method == "scree_plot"){
    # Method 3: Graphical representation of explained inertia (scree plot)
    # the user can set the threshold based on the scree plot

    df <- data.frame(dims = seq_len(max_num_dims),
                     inertia = expl_inertia)

    screeplot <- scree_plot(df)

    return(screeplot)

  } else if (method == "elbow_rule") {

    if(is.null(mat)){
      cat(paste0("When running method=\"elbow_rule\", ",
                 "please provide the original data matrix (paramater mat) ",
                 "which was earlier submitted to cacomp()!"))
      stop()
    }
    
    pd <- elbow_method(obj = obj,
                       mat = mat,
                       reps = reps,
                       python = python,
                       return_plot = return_plot)
    return(pd)
    
  } else {
    cat("Please pick a valid method!")
    stop()
  }
})



#' @param assay Character. The assay from which extract the count matrix for 
#' SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#' @param slot Character. Data slot of the Seurat assay. 
#' E.g. "data" or "counts". Default "counts".
#'
#' @rdname pick_dims
#' @export
#' @examples
#'
#' ################################
#' # pick_dims for Seurat objects #
#' ################################
#' library(Seurat)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Create Seurat object
#' seu <- CreateSeuratObject(counts = cnts)
#'
#' # run CA and save in dim. reduction slot.
#' seu <- cacomp(seu, return_input = TRUE, assay = "RNA", slot = "counts")
#'
#' # pick dimensions
#' pd <- pick_dims(obj = seu,
#'                 method = "maj_inertia",
#'                 assay = "RNA",
#'                 slot = "counts")
setMethod(f = "pick_dims",
          signature=(obj="Seurat"),
          function(obj,
                   mat = NULL,
                   method="scree_plot",
                   reps=3,
                   python = TRUE,
                   return_plot = FALSE,
                   ...,
                   assay = Seurat::DefaultAssay(obj),
                   slot = "counts"){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  if (method == "elbow_rule") {
    seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = slot)
    seu <- as.matrix(seu)
  } else {
    seu <- NULL
  }

  if ("CA" %in% Seurat::Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  } else {
    stop("No 'CA' dimension reduction object found. ",
         "Please run cacomp(seurat_obj, top, coords = FALSE, ",
         "return_input=TRUE) first.")
  }

  stopifnot(is(caobj, "cacomp"))

  pick_dims(obj = caobj,
            mat = seu,
            method = method,
            reps = reps,
            return_plot = return_plot,
            python = python)
})


#' @param assay Character. The assay from which to extract the count matrix 
#' for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#'
#' @rdname pick_dims
#' @export
#' @examples
#'
#' ##############################################
#' # pick_dims for SingleCellExperiment objects #
#' ##############################################
#' library(SingleCellExperiment)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Create SingleCellExperiment object
#' sce <- SingleCellExperiment(assays=list(counts=cnts))
#'
#' # run CA and save in dim. reduction slot.
#' sce <- cacomp(sce, return_input = TRUE, assay = "counts")
#'
#' # pick dimensions
#' pd <- pick_dims(obj = sce,
#'                 method = "maj_inertia",
#'                 assay = "counts")
setMethod(f = "pick_dims",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   mat = NULL,
                   method="scree_plot",
                   reps=3,
                   python = TRUE,
                   return_plot = FALSE,
                   ...,
                   assay = "counts"){

  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = 
              is(obj, "SingleCellExperiment"))

  if (method == "elbow_rule") {
    mat <- SummarizedExperiment::assay(obj, assay)
  } else {
    mat <- NULL
  }

  if ("CA" %in% SingleCellExperiment::reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  } else {
    stop("No 'CA' dim. reduction object found. ",
         "Please run cacomp(sce, top, coords = FALSE, ",
         "return_input=TRUE) first.")
  }

  stopifnot(is(caobj, "cacomp"))
  pick_dims(obj = caobj,
           mat = mat,
           method = method,
           reps = reps,
           return_plot = return_plot,
           python = python)

})







