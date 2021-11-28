
#' Helper function to check if object is empty.
#' @param x object
#' @return TRUE if x has length 0 and is not NULL. FALSE otherwise
is.empty <- function(x) return(isTRUE(length(x) == 0 & !is.null(x)))


#' Check if cacomp object was correctly created.
#'
#' @description Checks if the slots in a cacomp object are of the correct size
#' and whether they are coherent.
#' @param object A cacomp object.
#' @return TRUE if it is a valid cacomp object. FALSE otherwise.
#' @export
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
#'
#' check_cacomp(ca)
check_cacomp <- function(object) {
  errors <- character()

  dim_rows <- object@top_rows
  dims <- object@dims

  # SVD results
  if (isTRUE(!is.empty(object@U) & nrow(object@U) != dim_rows)) {
    msg <- paste0("Nr. of rows in U is ", nrow(object@U), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@U) & ncol(object@U) != dims)) {
    msg <- paste0("Nr. of columns in U is ", ncol(object@U), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@V) & ncol(object@V) != dims)) {
    msg <- paste0("Nr. of columns in V is ", ncol(object@V), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@D) & length(object@D) != dims)) {
    msg <- paste0("Length of D is ", ncol(object@D), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  # CA results

  if (isTRUE(!is.empty(object@row_masses) & length(object@row_masses) != dim_rows)) {
    msg <- paste0("Length of row_masses is ", length(object@row_masses), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@col_masses) & length(object@col_masses) != nrow(object@V))) {
    msg <- paste0("Length of col_masses is ", length(object@col_masses), ".  Should be ", nrow(object@V), ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@row_inertia) & length(object@row_inertia) != dim_rows)){
    msg <- paste0("Length of row_inertia is ", length(object@row_inertia), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@col_inertia) & length(object@col_inertia) != nrow(object@V))) {
    msg <- paste0("Length of col_inertia is ", length(object@col_inertia), ".  Should be ", nrow(object@V), ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@tot_inertia) & length(object@tot_inertia) != 1)) {
    msg <- paste0("Length of tot_inertia is ", length(object@tot_inertia), ".  Should be 1.")
    errors <- c(errors, msg)
  }

  # standardized coordinates

  if (isTRUE(!is.empty(object@std_coords_rows) & nrow(object@std_coords_rows) != dim_rows)) {
    msg <- paste0("Nr. of rows in std_coords_rows is ", nrow(object@std_coords_rows), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@std_coords_rows) & ncol(object@std_coords_rows) != dims)) {
    msg <- paste0("Nr. of columns in std_coords_rows is ", ncol(object@std_coords_rows), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@std_coords_cols) & nrow(object@std_coords_cols) != nrow(object@V))) {
    msg <- paste0("Nr. of rows in std_coords_cols is ", nrow(object@std_coords_cols), ".  Should be ", nrow(object@V), ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@std_coords_cols) & ncol(object@std_coords_cols) != dims)) {
    msg <- paste0("Nr. of columns in std_coords_cols is ", ncol(object@std_coords_cols), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }


  # principal coordinates

  if (isTRUE(!is.empty(object@prin_coords_rows) & nrow(object@prin_coords_rows) != dim_rows)) {
    msg <- paste0("Nr. of rows in prin_coords_rows is ", nrow(object@prin_coords_rows), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@prin_coords_rows) & ncol(object@prin_coords_rows) != dims)) {
    msg <- paste0("Nr. of columns in prin_coords_rows is ", ncol(object@prin_coords_rows), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@prin_coords_cols) & nrow(object@prin_coords_cols) != nrow(object@V))) {
    msg <- paste0("Nr. of rows in prin_coords_cols is ", nrow(object@prin_coords_cols), ".  Should be ", nrow(object@V), ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@prin_coords_cols) & ncol(object@prin_coords_cols) != dims)) {
    msg <- paste0("Nr. of columns in prin_coords_cols is ", ncol(object@prin_coords_cols), ".  Should be ", dims, ".")
    errors <- c(errors, msg)
  }

  # AP coordinates

  if (isTRUE(!is.empty(object@apl_rows) & nrow(object@apl_rows) != dim_rows)) {
    msg <- paste0("Nr. of rows in apl_rows is ", ncol(object@apl_rows), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@apl_rows) & ncol(object@apl_rows) != 2)) {
    msg <- paste0("Nr. of columns in apl_rows is ", ncol(object@apl_rows), ".  Should be 2.")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@apl_cols) & nrow(object@apl_cols) != nrow(object@V))) {
    msg <- paste0("Nr. of rows in apl_cols is ", ncol(object@apl_cols), ".  Should be ", nrow(object@V), ".")
    errors <- c(errors, msg)
  }

  if (isTRUE(!is.empty(object@apl_cols) & ncol(object@apl_cols) != 2)) {
    msg <- paste0("Nr. of columns in apl_cols is ", ncol(object@apl_cols), ".  Should be 2.")
    errors <- c(errors, msg)
  }

  # Salpha score
  if (isTRUE(!is.empty(object@APL_score) & ncol(object@APL_score) != 4)) {
    msg <- paste0("Nr. of columns in APL_score is ", ncol(object@APL_score), ".  Should be 4.")
    errors <- c(errors, msg)
  }
  if (isTRUE(!is.empty(object@APL_score) & nrow(object@APL_score) != dim_rows)) {
    msg <- paste0("Nr. of rows in APL_score is ", nrow(object@APL_score), ".  Should be ", dim_rows, ".")
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) TRUE else errors
}

#' An S4 class that contains all elements needed for CA.
#' @name cacomp-class
#' @rdname cacomp-class
#' @description
#' This class contains elements necessary to computer CA coordinates or Association Plot coordinates,
#' as well as other informative data such as row/column inertia, gene-wise APL-scores, etc. ...
#'
#' @slot U class "matrix". Left singular vectors of the original input matrix.
#' @slot V class "matrix". Right singular vectors of the original input matrix.
#' @slot D class "numeric". Singular values of the original inpt matrix.
#' @slot std_coords_rows class "matrix". Standardized CA coordinates of the rows.
#' @slot std_coords_cols class "matrix". Standardized CA coordinates of the columns.
#' @slot prin_coords_rows class "matrix". Principal CA coordinates of the rows.
#' @slot prin_coords_cols class "matrix". Principal CA coordinates of the columns.
#' @slot apl_rows class "matrix". Association Plot coordinates of the rows for the direction defined in slot "group"
#' @slot apl_cols class "matrix". Association Plot coordinates of the columns for the direction defined in slot "group"
#' @slot APL_score class "data.frame". Contains rows sorted by the APL score.
#' Columns: Rowname (gene name in the case of gene expression data),
#' APL score calculated for the direction defined in slot "group",
#' the original row number and the rank of the row as determined by the score.
#' @slot dims class "numeric". Number of dimensions in CA space.
#' @slot group class "numeric". Indices of the chosen columns for APL calculations.
#' @slot row_masses class "numeric". Row masses of the frequency table.
#' @slot col_masses class "numeric". Column masses of the frequency table.
#' @slot top_rows class "numeric". Number of most variable rows chosen.
#' @slot tot_inertia class "numeric". Total inertia in CA space.
#' @slot row_inertia class "numeric". Row-wise inertia in CA space.
#' @slot col_inertia class "numeric". Column-wise inertia in CA space.
#' @slot permuted_data class "list". Storage slot for permuted data.
#' @export
setClass("cacomp",
         representation(
           U = "matrix",
           V = "matrix",
           D = "numeric",
           std_coords_rows = "matrix",
           std_coords_cols = "matrix",
           prin_coords_rows = "matrix",
           prin_coords_cols = "matrix",
           apl_rows = "matrix",
           apl_cols = "matrix",
           APL_score = "data.frame",
           dims = "numeric",
           group = "numeric",
           row_masses = "numeric",
           col_masses = "numeric",
           top_rows = "numeric",
           tot_inertia = "numeric",
           row_inertia = "numeric",
           col_inertia = "numeric",
           permuted_data = "list"
         ),
         prototype(
           U = matrix(0, 0, 0),
           V = matrix(0, 0, 0),
           D = numeric(),
           std_coords_rows = matrix(0, 0, 0),
           std_coords_cols = matrix(0, 0, 0),
           prin_coords_rows = matrix(0, 0, 0),
           prin_coords_cols = matrix(0, 0, 0),
           apl_rows = matrix(0, 0, 0),
           apl_cols = matrix(0, 0, 0),
           APL_score = data.frame(),
           dims = numeric(),
           group = numeric(),
           row_masses = numeric(),
           col_masses = numeric(),
           top_rows = numeric(),
           tot_inertia = numeric(),
           row_inertia = numeric(),
           col_inertia = numeric(),
           permuted_data = list()),
         validity = check_cacomp
)

#' Create new "cacomp" object.
#' @description Creates new cacomp object.
#'
#' @param ... slot names and objects for new cacomp object.
#' @return cacomp object
#' @rdname cacomp-class
#' @export
#' @examples
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)}, x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' res <-  APL:::comp_std_residuals(mat=cnts)
#' SVD <- svd(res$S)
#' names(SVD) <- c("D", "U", "V")
#' SVD <- SVD[c(2, 1, 3)]
#'
#' ca <- new_cacomp(U = SVD$U,
#'                  V = SVD$V,
#'                  D = SVD$D,
#'                  row_masses = res$rowm,
#'                  col_masses = res$colm)
new_cacomp <- function(...) new("cacomp",...)

