

# Only here for testing of S3 methods and documentation of them.

#' testfun
#'
#' @export
#' @param x whatev
testfun <- function(x) {
  UseMethod("testfun")
}

#' testfun
#'
#' @export
#' @describeIn testf
#' @param x whatev
#' @param objm matrix
testfun.matrix <- function(x,objm){
  print("matrix")
}

#' testfun
#'
#' @export
#' @rdname testfun
#' @param x whatev
#' @param objd df
testfun.data.frame <- function(x, objd){
  print("data.frame")
}



findMethod <- function(generic, ...) {
  ch <- deparse(substitute(generic))
  f <- X <- function(x, ...) UseMethod("X")
  for(m in methods(ch)) assign(sub(ch, "X", m, fixed = TRUE), "body<-"(f, value = m))
  X(...)
}
