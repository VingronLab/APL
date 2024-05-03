
#' @import methods
#' @import SummarizedExperiment org.Hs.eg.db org.Mm.eg.db
#' @importFrom stats as.formula na.omit quantile runif var
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom ggplot2 ggplot aes geom_point guide_colorbar
#' @importFrom topGO showSigOfNodes score
#' @importFrom viridisLite viridis
#' @importFrom rlang .data
#' @importFrom irlba irlba
#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL


# Load python packages
.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
}
