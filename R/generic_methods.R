#' @include constructor.R
NULL

#' Prints cacomp object
#'
#' @description Provides more user friendly printing of cacomp objects.
#'
#' @param object cacomp object to print
#' @returns prints summary information about cacomp object.
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
#' ca
show.cacomp <- function(object){

  if (!is.empty(object@V) && !is.empty(object@U) && !is.empty(object@D)){
    cat("cacomp object with",
        nrow(object@V),
        "columns,",
        nrow(object@U),
        "rows and",
        length(object@D),
        "dimensions.")
  } else {
    cat("Uncomplete cacomp object.",
        "Consider running as.cacomp(object, recompute=TRUE).")
  }

  cat("\nCalc. standard coord.: ",
      paste0("std_coords_rows"[!is.empty(object@std_coords_rows)],
             ifelse(!is.empty(object@std_coords_rows) &&
                    !is.empty(object@std_coords_cols),
                    ", ",
                    ""),
            "std_coords_cols"[!is.empty(object@std_coords_cols)]))

  cat("\nCalc. principal coord.:",
      paste0("prin_coords_rows"[!is.empty(object@prin_coords_rows)],
             ifelse(!is.empty(object@prin_coords_rows) &&
                    !is.empty(object@prin_coords_cols),
                    ", ",
                    ""),
             "prin_coords_cols"[!is.empty(object@prin_coords_cols)]))


  cat("\nCalc. APL coord.:      ",
      paste0("apl_rows"[!is.empty(object@apl_rows)],
            ifelse(!is.empty(object@apl_rows) && !is.empty(object@apl_cols),
                   ", ",
                   ""),
            "apl_cols"[!is.empty(object@apl_cols)]))

  if (!is.empty(object@D)){
    prinInertia <- object@D^2
    percentInertia <- prinInertia / sum(prinInertia) * 100
    cat("\nExplained inertia:     ",
        paste0(round(percentInertia[1], 1),
               "% Dim1, ",
               round(percentInertia[2], 1),
               "% Dim2\n"))
  }

}

#' @rdname show.cacomp
#' @export
setMethod(f = "show", signature(object = "cacomp"), function(object) {
  show.cacomp(object)
})

#' Convert cacomp object to list.
#' @param x A cacomp object.
#' @return A cacomp object.
#' @export
#' @examples
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
setMethod("as.list",signature(x="cacomp"),function(x) {
  mapply(function(y) {

    if (inherits(slot(x,y),"cacomp")) {
      as.list(slot(x,y))
    } else {
      slot(x,y)
    }
  },
  slotNames(class(x)),
  SIMPLIFY=FALSE)
})
