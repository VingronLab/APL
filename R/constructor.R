

#' Create new object of class "cacomp"
#'
#' @param x object to create new cacomp from.
#' @export
new_cacomp <- function(x){
  stopifnot(is.list(x))

  cacomp_nms <- c("U",
                  "V",
                  "D",
                  "std_coords_rows",
                  "std_coords_cols",
                  "prin_coords_rows",
                  "prin_coords_cols",
                  "apl_rows",
                  "apl_cols",
                  "dims",
                  "group",
                  "row_masses",
                  "col_masses",
                  "top_rows",
                  "tot_inertia",
                  "row_inertia",
                  "col_inertia",
                  "permuted_data")

  check_names(x, canames = cacomp_nms)
  x <- fill_names(x, canames = cacomp_nms)
  x <- structure(x, class = "cacomp")

  return(x)
}

#' Check names of potential cacomp object.
#'
#' @param x An object to check.
#' @param canames Allowed names.
check_names <- function(x, canames){
  nms <- names(x)

  if (sum(!nms %in% canames) != 0){
    stop("x contains unsuitable names.")
  }

}

#' Create missing names
#' @description Adds missing names to object and sets them NULL.
#' @param x Object to fill names in for.
#' @param canames Allowed names.
fill_names <- function(x, canames){

  idx <- which(!canames %in% names(x))
  x <- c(x, setNames(rep(list(NULL), length(idx)), canames[idx]))
  x <- x[canames]
  return(x)
}
