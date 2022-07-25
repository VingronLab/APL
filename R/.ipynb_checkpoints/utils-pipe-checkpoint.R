#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return \code{magrittr::\link[magrittr:pipe]{\%>\%}}
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @examples
#' x <- 1:100
#' x %>% head()
NULL
