

#' Prints cacomp object
#'
#' @description Provides more user friendly printing of cacomp objects.
#'
#' @param object cacomp object to print
#' @export
show.cacomp <- function(object){

  if (!is.null(object@V) && !is.null(object@U) && !is.null(object@D)){
    cat("cacomp object with", nrow(object@V), "columns,", nrow(object@U), "rows and", length(object@D), "dimensions.")
  } else {
    cat("Uncomplete cacomp object. Consider running as.cacomp(object, recompute=TRUE).")
  }

  cat("\nCalc. standard coord.: ", paste0("std_coords_rows"[!is.null(object@std_coords_rows)],
                                         ifelse(!is.null(object@std_coords_rows) && !is.null(object@std_coords_cols), ", ", ""),
                                         "std_coords_cols"[!is.null(object@std_coords_cols)]))

  cat("\nCalc. principal coord.:", paste0("prin_coords_rows"[!is.null(object@prin_coords_rows)],
                                         ifelse(!is.null(object@prin_coords_rows) && !is.null(object@prin_coords_cols), ", ", ""),
                                         "prin_coords_cols"[!is.null(object@prin_coords_cols)]))


  cat("\nCalc. APL coord.:      ", paste0("apl_rows"[!is.null(object@apl_rows)],
                                          ifelse(!is.null(object@apl_rows) && !is.null(object@apl_cols), ", ", ""),
                                          "apl_cols"[!is.null(object@apl_cols)]))

  if (!is.null(object@D)){
    prinInertia <- object@D^2
    percentInertia <- prinInertia / sum(prinInertia) * 100
    cat("\nExplained inertia:     ", paste0(round(percentInertia[1], 1), "% Dim1, ", round(percentInertia[2], 1), "% Dim2"))
  }

}

setMethod(f = "show", signature(object = "cacomp"), function(object) {
  show.cacomp(object)
})


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
