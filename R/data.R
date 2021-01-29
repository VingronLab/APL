#' Subset of GTEx tissue dataset.
#'
#' A count matrix consisting of a subset of the
#' Genotype-Tissue Expression (GTEx) data.
#' The matrix contains count data from 22 tissues:
#' `unique(gsub("[[:digit:]]", "", colnames(gtex)))`
#'
#' @format A matrix with 56202 rows and 220 samples:
#'
#' @source \url{https://gtexportal.org/home/datasets}
"gtex"
