

#' Create new object of class "cacomp"
#'
#' @param x object to create new cacomp from.
new_cacomp <- function(x){
  stopifnot(is.list(x))
  structure(x, class = "cacomp")
}
