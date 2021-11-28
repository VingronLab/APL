#' @include constructor.R
NULL


#' Recompute missing values of cacomp object.
#'
#' @description
#' The caobj needs to have the std_coords_cols, the prin_coords_rows and D calculated. From this the remainder will be calculated.
#' Future updates might extend this functionality.
#'
#' @return
#' A cacomp object with additional calculated row_masses, col_masses, std_coords_rows, U and V.
#'
#' @param calist A list with std_coords_cols, the prin_coords_rows and D.
#' @param mat A matrix from which the cacomp object is derived from.
#' @param rm_zeros Removes rows & columns containing only zeros.
#' @param ... Further arguments forwarded to cacomp.
recompute <- function(calist, mat, rm_zeros = TRUE, ...){
  stopifnot(is(calist, "list"))
  stopifnot(is(mat, "matrix"))

  if(isTRUE(rm_zeros)){
    mat <- rm_zeros(mat)
  }

  # make stock of what we have

  std_rows <- is.null(calist$std_coords_rows)
  std_cols <- is.null(calist$std_coords_cols)
  prin_rows <- is.null(calist$prin_coords_rows)
  prin_cols <- is.null(calist$prin_coords_cols)

  sp_rows <- std_rows & prin_rows
  sp_cols <- std_cols & prin_cols

  d <- is.null(calist$D)
  v <- is.null(calist$V)
  u <- is.null(calist$U)

  # mat <- var_rows(mat = mat,
                  # top = nrow(mat))
  res <- comp_std_residuals(mat=mat)

  S <- res$S
  tot <- res$tot
  rowm <- res$rowm
  colm <- res$colm

  if(std_rows & !u) {
    calist$std_coords_rows <- sweep(x = calist$U, MARGIN = 1, STATS = sqrt(rowm), FUN = "/")
    std_rows <- FALSE
  }
  if(std_cols & !v){
    calist$std_coords_cols <- sweep(x = calist$V, MARGIN = 1, STATS = sqrt(colm), FUN = "/")
    std_cols <- FALSE
  }

  call_svd <- FALSE
  done <- FALSE

  while (isFALSE(done)){
    if (std_cols){
      if (d){
        if(prin_cols){
          call_svd <- TRUE
          done <- TRUE
        } else {
          # check if we can get D with row coords, otherwise call cacomp
          if(std_rows | prin_rows){
            call_svd <- TRUE
            done <- TRUE

          } else {
            calist$D <- calist$prin_coords_rows[1,]/calist$std_coords_rows[1,]
            d <- FALSE
          }
        }
      } else if (prin_cols){
        call_svd <- TRUE
        done <- TRUE

      } else {
        # calculate std_coords
        calist$std_coords_cols <- sweep(calist$prin_coords_cols, 2, calist$D, "/")
        std_cols <- FALSE
      }
    } else if (d) {
      if (prin_cols) {
        # check if we can get d through rows, otherweise cacomp
        if(std_rows | prin_rows){
          call_svd <- TRUE
          done <- TRUE

        } else {
          calist$D <- calist$prin_coords_rows[1,]/calist$std_coords_rows[1,]
          d <- FALSE
        }
      } else {
        # calculate d from col coordinates

        # calist$D <- colMeans(sweep(calist$prin_coords_cols, 1, calist$std_coords_cols, "/"))
        calist$D <- calist$prin_coords_cols[1,]/calist$std_coords_cols[1,]
        d <- FALSE

      }
    } else if (prin_cols){
      # calculate prin_cols with D and std
      calist$prin_coords_cols <- sweep(calist$std_coords_cols, 2, calist$D, "*")
      prin_cols <- FALSE

    } else {
      # all calculated
      done <- TRUE
    }
  }


  done <- FALSE
  while (isFALSE(done)){
    if (std_rows){
      if (d){
        if(prin_rows){
          call_svd <- TRUE
          done <- TRUE

        } else {
          # check if we can get D with row coords, otherwise call cacomp
          if(std_cols | prin_cols){
            call_svd <- TRUE
            done <- TRUE

          } else {
            calist$D <- calist$prin_coords_cols[1,]/calist$std_coords_cols[1,]
            d <- FALSE
          }
        }
      } else if (prin_rows){
        call_svd <- TRUE
        done <- TRUE
      } else {
        # calculate std_coords
        calist$std_coords_rows <- sweep(calist$prin_coords_rows, 2, calist$D, "/")
        std_rows <- FALSE
      }
    } else if (d) {
      if (prin_rows) {
        # check if we can get d through rows, otherweise cacomp
        if(std_cols | prin_cols){
          call_svd <- TRUE
          done <- TRUE

        } else {
          calist$D <- calist$prin_coords_cols[1,]/calist$std_coords_cols[1,]
          d <- FALSE
        }
      } else {
        # calculate d from col coordinates

        # calist$D <- colMeans(sweep(calist$prin_coords_rows, 1, calist$std_coords_rows, "/"))
        calist$D <- calist$prin_coords_rows[1,]/calist$std_coords_rows[1,]
        d <- FALSE

      }
    } else if (prin_rows){
      # calculate prin_rows with D and std
      calist$prin_coords_rows <- sweep(calist$std_coords_rows, 2, calist$D, "*")
      prin_rows <- FALSE

    } else {
      # all calculated
      done <- TRUE
    }
  }

  if(isTRUE(call_svd)){
    message("Calling cacomp to recompute from matrix.")
    ca <- cacomp(mat, princ_coords = 3, ...)
    return(ca)
  } else {

    if (nrow(mat) != nrow(calist$std_coords_rows)){
      stop("mat does not have have the correct number of rows.")
    }

    if (ncol(mat) != nrow(calist$std_coords_cols)){
      stop("mat does not have have the correct number of columns.")
    }

    calist$std_coords_rows[is.na(calist$std_coords_rows)] <- 0
    calist$std_coords_cols[is.na(calist$std_coords_cols)] <- 0
    calist$std_coords_rows[is.infinite(calist$std_coords_rows)] <- 0
    calist$std_coords_cols[is.infinite(calist$std_coords_cols)] <- 0

    ordidx <- match(rownames(calist$prin_coords_rows), names(rowm))
    calist$row_masses <- rowm[ordidx]

    ordidx <- match(rownames(calist$std_coords_cols), names(colm))
    calist$col_masses <- colm[ordidx]

    if (u) calist$U <- sweep(calist$std_coords_rows, 1, sqrt(calist$row_masses), "*")
    if (v) calist$V <- sweep(calist$std_coords_cols, 1, sqrt(calist$col_masses), "*")

    calist$tot_inertia <- sum(calist$D^2)
    calist$row_inertia <- rowSums(S^2)
    calist$col_inertia <- colSums(S^2)

    calist$top_rows <- nrow(mat)
    calist$dims <- length(calist$D)
  }

  ca <- do.call(new_cacomp, calist)
  return(ca)
}


#' Create cacomp object from Seurat/SingleCellExperiment container
#'
#' @description
#' Converts the values stored in the Seurat/SingleCellExperiment dimensional reduction slot "CA" to a cacomp object.
#' If recompute = TRUE additional parameters are recomputed from the saved values without rerunning SVD (need to specify assay to work).
#'
#' @details
#' By default extracts std_coords_cols, D, prin_coords_rows, top_rows and dims from obj and outputs a cacomp object.
#' If recompute = TRUE the following are additionally recalculated (doesn't run SVD):
#' U, V, std_coords_rows, row_masses, col_masses.
#'
#' @return
#' A cacomp object.
#'
#' @param obj An object of class "Seurat" or "SingleCellExperiment"
#' with a dim. reduction named "CA" saved. For obj "cacomp" input is returned.
#' @param assay Character. The assay from which extract the count matrix,
#' e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#' @param ... Further arguments.
#' @export
setGeneric("as.cacomp", function(obj, ...) {
  standardGeneric("as.cacomp")
})

#' @description as.cacomp.cacomp returns input without any calculations.
#' @rdname as.cacomp
#' @export
setMethod(f = "as.cacomp", signature=(obj="cacomp"), function(obj, ...) {
  stopifnot(is(obj, "cacomp"))
  return(obj)
})


#' @description Recomputes missing values and returns cacomp object from a list.
#' If you have a *complete* cacomp object in list form,
#' use do.call(new_cacomp, obj).
#' @param mat Original input matrix.
#' @rdname as.cacomp
#' @export
#' @examples
#' #########
#' # lists #
#' #########
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#' ca_list <- as.list(ca)
#'
#' # Only keep subset of elements for demonstration
#' ca_list <- ca_list[c("U", "std_coords_rows", "std_coords_cols")]
#'
#' # convert (incomplete) list to cacomp object.
#' ca <- as.cacomp(ca_list, mat = cnts)
setMethod(f = "as.cacomp", signature=(obj="list"), function(obj, ..., mat = NULL) {

  try_obj <- try(do.call(new_cacomp, ca_list), silent = TRUE)
  if (is(try_obj, "try-error")){
    obj <- recompute(calist = obj, mat = mat)
    return(obj)
  } else if (is(try_obj, "cacomp")){
    return(try_obj)
  } else {
    stop("Unexpected output from try().")
  }
})

#' @description
#' as.cacomp.Seurat: Converts the values stored in the Seurat DimReduc slot "CA" to an cacomp object.
#' @param slot character. Slot of the Seurat assay to use. Default "counts".
#' @rdname as.cacomp
#' @export
#' @examples
#'
#' ##########
#' # Seurat #
#' ##########
#' library(Seurat)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' seu <- CreateSeuratObject(counts = cnts)
#' seu <- cacomp(seu, return_input = TRUE)
#'
#' ca <- as.cacomp(seu, assay = "RNA", slot = "counts")
setMethod(f = "as.cacomp", signature=(obj="Seurat"), function(obj, ..., assay="RNA", slot = "counts") {

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))
  stopifnot("obj doesn't contain a DimReduc object named 'CA'. Try running cacomp()." = "CA" %in% names(obj@reductions))

  if (is.null(assay)) assay <- DefaultAssay(obj)

  ca_list <- list("std_coords_cols" = Seurat::Embeddings(obj, reduction = "CA"),
                  "D" = Seurat::Stdev(obj, reduction = "CA"),
                  "prin_coords_rows" = Seurat::Loadings(obj, reduction = "CA"))
  ca_list$top_rows <- nrow(ca_list$prin_coords_rows)
  ca_list$dims <- length(ca_list$D)

  colnames(ca_list$std_coords_cols) <-  paste0("Dim", seq_len(ncol(ca_list$std_coords_cols)))
  colnames(ca_list$prin_coords_rows) <-  paste0("Dim", seq_len(ncol(ca_list$prin_coords_rows)))
  names(ca_list$D) <-  paste0("Dim", seq_len(length(ca_list$D)))

  stopifnot("Assay is needed to recompute cacomp." = !is.null(assay))

  seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = slot)
  seu <- as.matrix(seu)
  seu <- seu[rownames(ca_list$prin_coords_rows),]

  ca_obj <- recompute(calist = ca_list, mat = seu)

  # ca_obj <- do.call(new_cacomp, ca_list)

  stopifnot(validObject(ca_obj))
  return(ca_obj)
})


#' @description
#' as.cacomp.SingleCellExperiment: Converts the values stored in the SingleCellExperiment reducedDim slot "CA" to a cacomp object.
#'
#' @rdname as.cacomp
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
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' sce <- SingleCellExperiment(assays=list(counts=cnts))
#' sce <- cacomp(sce, return_input = TRUE)
#'
#' ca <- as.cacomp(sce, assay = "counts")
setMethod(f = "as.cacomp",
          signature=(obj="SingleCellExperiment"),
          function(obj, ..., assay="counts") {

  sce_ca <- SingleCellExperiment::reducedDim(obj, "CA")
  stopifnot("Attribute singval of dimension reduction slot CA is empty.\nThis can happen after subsetting the sce obj." = !is.null(attr(sce_ca, "singval")))
  stopifnot("Attribute prin_coords_rows of dimension reduction slot CA is empty.\nThis can happen after subsetting the sce obj." = !is.null(attr(sce_ca, "prin_coords_rows")))

  ca_list <- list("std_coords_cols" = sce_ca,
                  "D" = attr(sce_ca, "singval"),
                  "prin_coords_rows" = attr(sce_ca, "prin_coords_rows"))

  if(is.null(assay)) assay <- "counts"

  attr(ca_list$std_coords_cols, "prin_coords_rows") <- NULL
  attr(ca_list$std_coords_cols, "singval") <- NULL
  attr(ca_list$std_coords_cols, "percInertia") <- NULL

  ca_list$top_rows <- nrow(ca_list$prin_coords_rows)
  ca_list$dims <- length(ca_list$D)


  stopifnot("Assay is needed to recompute cacomp." = !is.null(assay))
  scemat <- SummarizedExperiment::assay(obj, assay)
  scemat <- scemat[rownames(ca_list$prin_coords_rows),]

  ca_obj <- recompute(calist = ca_list, mat = scemat)


  # ca_obj <- do.call(new_cacomp, ca_list)

  stopifnot(validObject(ca_obj))
  return(ca_obj)
})

#' Convert SingleCellObject object to cacomp object.
# setMethod("as.cacomp", "SingleCellExperiment", as.cacomp.SingleCellExperiment)
