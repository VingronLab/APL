

#' Prints cacomp object
#'
#' @description Provides more user friendly printing of cacomp objects.
#'
#' @param x cacomp object to print
#' @param ... ignored.
#' @export
print.cacomp <- function(x, ...){

  if (!is.null(x$V) && !is.null(x$U) && !is.null(x$D)){
    cat("cacomp object with", nrow(x$V), "columns,", nrow(x$U), "rows and", length(x$D), "dimensions.")
  } else {
    cat("Uncomplete cacomp object. Consider running as.cacomp(x, recompute=TRUE).")
  }

  cat("\nCalc. standard coord.: ", paste0("std_coords_rows"[!is.null(x$std_coords_rows)],
                                         ifelse(!is.null(x$std_coords_rows) && !is.null(x$std_coords_cols), ", ", ""),
                                         "std_coords_cols"[!is.null(x$std_coords_cols)]))

  cat("\nCalc. principal coord.:", paste0("prin_coords_rows"[!is.null(x$prin_coords_rows)],
                                         ifelse(!is.null(x$prin_coords_rows) && !is.null(x$prin_coords_cols), ", ", ""),
                                         "prin_coords_cols"[!is.null(x$prin_coords_cols)]))


  cat("\nCalc. APL coord.:      ", paste0("apl_rows"[!is.null(x$apl_rows)],
                                          ifelse(!is.null(x$apl_rows) && !is.null(x$apl_cols), ", ", ""),
                                          "apl_cols"[!is.null(x$apl_cols)]))

  if (!is.null(x$D)){
    prinInertia <- x$D^2
    percentInertia <- prinInertia / sum(prinInertia) * 100
    cat("\nExplained inertia:     ", paste0(round(percentInertia[1], 1), "% Dim1, ", round(percentInertia[2], 1), "% Dim2"))
  }

}
