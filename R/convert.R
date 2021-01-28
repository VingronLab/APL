
#' Recompute missing values of cacomp object.
#'
#' @description
#' The caobj needs to have the std_coords_cols, the prin_coords_rows and D calculated. From this the remainder will be calculated.
#' Future updates might extend this functionality.
#'
#' @return
#' A cacomp object with additional calculated row_masses, col_masses, std_coords_rows, U and V.
#'
#' @param caobj A cacomp object with std_coords_cols, the prin_coords_rows and D.
#' @param mat A matrix from which the cacomp object derives from.
recompute <- function(caobj, mat){
  stopifnot(is(caobj, "cacomp"))
  stopifnot(is(mat, "matrix"))

  mat <- var_rows(mat = mat,
                  top = caobj$top_rows)
  res <-  comp_std_residuals(mat=mat)

  S <- res$S
  tot <- res$tot
  rowm <- res$rowm
  colm <- res$colm

  ordidx <- match(rownames(caobj$prin_coords_rows), names(rowm))
  caobj$row_masses <- rowm[ordidx]

  ordidx <- match(rownames(caobj$std_coords_cols), names(colm))
  caobj$col_masses <- colm[ordidx]

  caobj$std_coords_rows <- sweep(caobj$prin_coords_rows, 2, caobj$D, "/")
  caobj$U <- sweep(caobj$std_coords_rows, 1, sqrt(caobj$row_masses), "*")
  caobj$V <- sweep(caobj$std_coords_cols, 1, sqrt(caobj$col_masses), "*")

  return(caobj)
}


#' Create cacomp object from Seurat container
#'
#' @description
#' Converts the values stored in the Seurat DimReduc slot "CA" to an cacomp object.
#' If recompute = TRUE additional parameters are recomputed from the saved values without rerunning SVD.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "cacomp", "Seurat" or "SingleCellExperiment" with a DimReduc Object in the "CA" slot.
#' @param assay The name of the assay to use, e.g. "RNA".
#' @param recompute Should all values be recomputed (std_coords_rows, U and V as well as row and column masses)? Default TRUE.
#' @export
as.cacomp <- function(obj, assay, recompute = TRUE) UseMethod("as.cacomp")

#' Create cacomp object from Seurat container
#'
#' @description
#' Converts the values stored in the Seurat DimReduc slot "CA" to an cacomp object.
#' If recompute = TRUE additional parameters are recomputed from the saved values without rerunning SVD.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "Seurat" with a DimReduc Object in the "CA" slot.
#' @param assay The name of the assay to use, e.g. "RNA".
#' @param recompute Should all values be recomputed (std_coords_rows, U and V as well as row and column masses)? Default TRUE.
#' @export
as.cacomp.default <- function(obj, assay, recompute = TRUE){
  stop(paste0("as.cacomp does not know how to handle objects of class ",
              class(obj),
              ". Currently only objects of class 'Seurat' or 'SingleCellExperiment' can be converted to 'cacomp'."))
}

#' Convert to cacomp object
#'
#' @description
#' If an object already is of class cacomp, the function simply returns the input.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "cacomp".
#' @export
as.cacomp.cacomp <- function(obj){
  stopifnot(is(obj, "cacomp"))
  return(obj)
}



#' Create cacomp object from Seurat container
#'
#' @description
#' Converts the values stored in the Seurat DimReduc slot "CA" to an cacomp object.
#' If recompute = TRUE additional parameters are recomputed from the saved values without rerunning SVD.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "Seurat" with a DimReduc Object in the "CA" slot.
#' @param assay The name of the assay to use, e.g. "RNA".
#' @param recompute Should all values be recomputed (std_coords_rows, U and V as well as row and column masses)? Default TRUE.
#' @export
as.cacomp.Seurat <- function(obj, assay=NULL, recompute = TRUE){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))
  stopifnot("obj doesn't contain a DimReduc object named 'CA'. Try running cacomp()." = "CA" %in% Reductions(pbmc_small))

  ca_obj <- list("std_coords_cols" = Seurat::Embeddings(obj, reduction = "CA"),
                 "D" = Seurat::Stdev(obj, reduction = "CA"),
                 "prin_coords_rows" = Seurat::Loadings(obj, reduction = "CA"))
  ca_obj$top_rows <- nrow(ca_obj$prin_coords_rows)
  ca_obj$dims <- length(ca_obj$D)

  class(ca_obj) <- "cacomp"

  if (recompute == TRUE){
    stopifnot("Assay is needed to recompute cacomp." = !is.null(assay))

    seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = "data")
    seu <- as.matrix(seu)
    # res <-  comp_std_residuals(mat=seu)

    ca_obj <- recompute(caobj = ca_obj, mat = seu)
    # seu <- var_rows(mat = seu,
    #                    top = ca_obj$top_rows)
    #
    #
    # res <-  comp_std_residuals(mat=seu)
    # S <- res$S
    # tot <- res$tot
    # rowm <- res$rowm
    # colm <- res$colm
    #
    # ordidx <- match(rownames(ca_obj$prin_coords_rows), names(rowm))
    # ca_obj$row_masses <- rowm[ordidx]
    #
    # ordidx <- match(rownames(ca_obj$std_coords_cols), names(colm))
    # ca_obj$col_masses <- colm[ordidx]
    #
    # ca_obj$std_coords_rows <- sweep(ca_obj$prin_coords_rows, 2, ca_obj$D, "/")
    # ca_obj$U <- sweep(ca_obj$std_coords_rows, 1, sqrt(ca_obj$row_masses), "*")
    # ca_obj$V <- sweep(ca_obj$std_coords_cols, 1, sqrt(ca_obj$col_masses), "*")
  }
  return(ca_obj)
}








#' Create cacomp object from SingleCellExperiment object.
#'
#' @description
#' Converts the values stored in the SingleCellExperiment reducedDim slot "CA" to an cacomp object.
#' If recompute = TRUE additional parameters are recomputed from the saved values without rerunning SVD.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "SingleCellExperiment" with a LinearEmbeddingMatrix object in the reducedDim(sce, "CA") slot.
#' @param assay The name of the assay to use, e.g. "RNA".
#' @param recompute Should all values be recomputed (std_coords_rows, U and V as well as row and column masses)? Default TRUE.
#' @export
as.cacomp.SingleCellExperiment <- function(obj, assay = NULL, recompute = TRUE){

  sce_ca <- SingleCellExperiment::reducedDim(sce, "CA")
  stopifnot("Attribute singval of dimensional reduction slot CA is empty.\nThis can happen after subsetting the sce obj." = !is.null(attr(sce_ca, "singval")))
  stopifnot("Attribute prin_coords_rows of dimensional reduction slot CA is empty.\nThis can happen after subsetting the sce obj." = !is.null(attr(sce_ca, "prin_coords_rows")))

  ca_obj <- list("std_coords_cols" = sce_ca,
                 "D" = attr(sce_ca, "singval"),
                 "prin_coords_rows" = attr(sce_ca, "prin_coords_rows"))


  attr(ca_obj$std_coords_cols, "prin_coords_rows") <- NULL
  attr(ca_obj$std_coords_cols, "singval") <- NULL
  attr(ca_obj$std_coords_cols, "percInertia") <- NULL

  ca_obj$top_rows <- nrow(ca_obj$prin_coords_rows)
  ca_obj$dims <- length(ca_obj$D)

  class(ca_obj) <- "cacomp"

  if (recompute == TRUE){

    stopifnot("Assay is needed to recompute cacomp." = !is.null(assay))
    scemat <- SingleCellExperiment::assay(sce, assay)
    scemat
    # res <-  comp_std_residuals(mat=scemat)

    ca_obj <- recompute(caobj = ca_obj, mat = scemat)

    # scemat <- var_rows(mat = scemat,
    #                    top =  ca_obj$top_rows)
    #
    # res <-  comp_std_residuals(mat=scemat)
    # S <- res$S
    # tot <- res$tot
    # rowm <- res$rowm
    # colm <- res$colm
    #
    #
    # ordidx <- match(rownames(ca_obj$prin_coords_rows), names(rowm))
    # ca_obj$row_masses <- rowm[ordidx]
    #
    # ordidx <- match(rownames(ca_obj$std_coords_cols), names(colm))
    # ca_obj$col_masses <- colm[ordidx]
    #
    # ca_obj$std_coords_rows <- sweep(ca_obj$prin_coords_rows, 2, ca_obj$D, "/")
    # ca_obj$U <- sweep(ca_obj$std_coords_rows, 1, sqrt(ca_obj$row_masses), "*")
    # ca_obj$V <- sweep(ca_obj$std_coords_cols, 1, sqrt(ca_obj$col_masses), "*")
  }

  return(ca_obj)
}



