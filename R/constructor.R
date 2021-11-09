



is.empty <- function(x) return(isTRUE(length(x) == 0 & !is.null(x)))

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

  # standard coordinates

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

  # APL coordinates

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

  # APL score
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
           group = "integer",
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
           group = integer(),
           row_masses = numeric(),
           col_masses = numeric(),
           top_rows = numeric(),
           tot_inertia = numeric(),
           row_inertia = numeric(),
           col_inertia = numeric(),
           permuted_data = list()),
         validity = check_cacomp
)



new_cacomp <- function(...) new("cacomp",...)


#'
#'
#'
#' #' Create new object of class "cacomp"
#' #'
#' #' @param x object to create new cacomp from.
#' #' @export
#' new_cacomp <- function(x){
#'   stopifnot(is.list(x))
#'
#'   cacomp_nms <- c("U",
#'                   "V",
#'                   "D",
#'                   "std_coords_rows",
#'                   "std_coords_cols",
#'                   "prin_coords_rows",
#'                   "prin_coords_cols",
#'                   "apl_rows",
#'                   "apl_cols",
#'                   "dims",
#'                   "group",
#'                   "row_masses",
#'                   "col_masses",
#'                   "top_rows",
#'                   "tot_inertia",
#'                   "row_inertia",
#'                   "col_inertia",
#'                   "permuted_data")
#'
#'   check_names(x, canames = cacomp_nms)
#'   x <- fill_names(x, canames = cacomp_nms)
#'   x <- structure(x, class = "cacomp")
#'
#'   return(x)
#' }
#'
#' #' Check names of potential cacomp object.
#' #'
#' #' @param x An object to check.
#' #' @param canames Allowed names.
#' check_names <- function(x, canames){
#'   nms <- names(x)
#'
#'   if (sum(!nms %in% canames) != 0){
#'     stop("x contains unsuitable names.")
#'   }
#'
#' }
#'
#' #' Create missing names
#' #' @description Adds missing names to object and sets them NULL.
#' #' @param x Object to fill names in for.
#' #' @param canames Allowed names.
#' fill_names <- function(x, canames){
#'
#'   idx <- which(!canames %in% names(x))
#'   x <- c(x, setNames(rep(list(NULL), length(idx)), canames[idx]))
#'   x <- x[canames]
#'   return(x)
#' }
