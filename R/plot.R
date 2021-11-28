#' @include constructor.R
NULL

#' Plot of the first 3D CA projection of the data.
#'
#' @description
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#'
#' @details
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standardized coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and column names respectively.
#' @return
#' Plot of class "plotly".
#'
#' @param obj  An object of class "cacomp", or alternatively an object of class "Seurat" or "SingleCellExperiment" with a dim. reduction named "CA" saved.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which
#' a label should be added (label should be stored in colnames).
#' Default NULL (no columns).
#' @param ... Further arguments.
#'
#' @export
#' @examples
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#'
#' ca_3Dplot(ca)
setGeneric("ca_3Dplot", function(obj,
                              xdim = 1,
                              ydim = 2,
                              zdim = 3,
                              princ_coords = 1,
                              row_labels = NULL,
                              col_labels = NULL,
                              ...) {
  standardGeneric("ca_3Dplot")
})


#' @rdname ca_3Dplot
#' @export
setMethod(f = "ca_3Dplot",
          signature=(obj="cacomp"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   zdim = 3,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   ...){

  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }

  if (princ_coords == 1){

    if(sum(!is.null(obj@prin_coords_rows), !is.null(obj@std_coords_cols)) != 2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj@prin_coords_rows
    cols <- obj@std_coords_cols
  } else if (princ_coords == 2){
    if(sum(!is.null(obj@prin_coords_cols), !is.null(obj@std_coords_rows)) != 2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj@std_coords_rows
    cols <- obj@prin_coords_cols
  } else {
    cat("princ_coords must be either 1 for rows or 2 for columns.")
    stop()
  }

  p <- plotly::plot_ly() %>%
    plotly::add_trace(x = cols[,xdim],
              y = cols[,ydim],
              z = cols[,zdim],
              mode = 'markers',
              text = rownames(cols),
              textposition = "left",
              opacity = 1,
              marker = list(color = '#990000', symbol = 'x', size = 2),
              name = 'Columns',
              hoverinfo = 'text',
              type = 'scatter3d') %>%
    plotly::add_trace(x = rows[,xdim],
              y = rows[,ydim],
              z = rows[,zdim],
              mode = 'markers',
              text = rownames(rows),
              opacity = 1,
              marker = list(color ='#0066FF', symbol = 'circle-open', size = 1),
              name = 'genes',
              hoverinfo = 'text',
              type = 'scatter3d')

  if (!is.null(row_labels)){
    p <- p %>%
      plotly::add_trace(x = rows[row_labels, xdim],
                y = rows[row_labels, ydim],
                z = rows[row_labels, zdim],
                mode = 'markers+text',
                text = rownames(rows)[row_labels],
                textposition = "left",
                textfont = list(color='#FF0000'),
                marker = list(symbol = 'circle', color ='#FF0000', size = 2),
                name = 'marked row(s)',
                hoverinfo = 'text',
                type = 'scatter3d')
  }

  if (!is.null(col_labels)){
    p <- p %>%
      plotly::add_trace(x = cols[col_labels, xdim],
                y = cols[col_labels, ydim],
                z = cols[col_labels, zdim],
                mode = 'markers+text',
                text = rownames(cols)[col_labels],
                textposition = "left",
                textfont = list(color='#990000'),
                marker = list(symbol = 'circle-open', color ='#990000', size = 3),
                name = 'marked column(s)',
                hoverinfo = 'text',
                type = 'scatter3d')
  }

  axx <- list(title =  paste0('Dim',xdim))

  axy <- list(title =  paste0('Dim',ydim))

  axz <- list(title =  paste0('Dim',zdim))
  p <- p %>%
    plotly::layout(autosize = TRUE,
           title = '3D CA plot',
           showlegend = FALSE,
           scene = list(xaxis = axx, yaxis = axy, zaxis = axz))

  return(p)

})


#' @rdname ca_3Dplot
#' @param assay Assay to use to recompute cacomp obj.
#' @param slot Seurat slot from assay to get count matrix from.
#' @export
setMethod(f = "ca_3Dplot",
          signature=(obj="Seurat"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   zdim = 3,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   ...,
                   assay = Seurat::DefaultAssay(obj),
                   slot = "counts"){
  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  if ("CA" %in% Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay, slot = slot)
  } else {
    stop("No 'CA' dimension reduction object found. Please run cacomp(seurat_obj, assay) first.")
  }

  p <- ca_3Dplot(obj = caobj,
                 xdim = xdim,
                 ydim = ydim,
                 zdim = zdim,
                 princ_coords = princ_coords,
                 row_labels = row_labels,
                 col_labels = col_labels)
  return(p)
})



#' @rdname ca_3Dplot
#' @param assay SingleCellExperiment assay to obtain counts from.
#' @export
setMethod(f = "ca_3Dplot",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   zdim = 3,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   ...,
                   assay = "counts"){
  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))

  if ("CA" %in% reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  } else {
    stop("No 'CA' dimension reduction object found. Please run cacomp(sce, top, coords = FALSE, return_input=TRUE) first.")
  }

  p <- ca_3Dplot(obj = caobj,
                xdim = xdim,
                ydim = ydim,
                zdim = zdim,
                princ_coords = princ_coords,
                row_labels = row_labels,
                col_labels = col_labels)
  return(p)
})

#' Plot of 2D CA projection of the data.
#'
#' @description
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#'
#' @details
#' Choosing type "plotly" will generate an interactive html plot with the package plotly.
#' Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either
#' the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return
#' Plot of class "plotly" or "ggplot".
#'
#' @param obj An object of class "cacomp" with the relevant standardized and principal coordinates calculated,
#'  or alternatively an object of class "Seurat" or "SingleCellExperiment" with a dim. reduction named "CA" saved.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows,
#' if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added
#' (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added
#' (label should be stored in colnames).
#' Default NULL (no columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
#' @param ... Further arguments.
#' @export
#' @examples
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#'
#' ca_biplot(ca)
setGeneric("ca_biplot", function(obj,
                                 xdim = 1,
                                 ydim = 2,
                                 princ_coords = 1,
                                 row_labels = NULL,
                                 col_labels = NULL,
                                 type = "plotly",
                                 ...) {
  standardGeneric("ca_biplot")
})



#' @rdname ca_biplot
#' @export
setMethod(f = "ca_biplot",
          signature=(obj="cacomp"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   type = "plotly",
                   ...){

  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }

  if (princ_coords == 1){

    if(sum(!is.null(obj@prin_coords_rows), !is.null(obj@std_coords_cols)) != 2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj@prin_coords_rows
    cols <- obj@std_coords_cols
  } else if (princ_coords == 2){
    if(sum(!is.null(obj@prin_coords_cols), !is.null(obj@std_coords_rows)) != 2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj@std_coords_rows
    cols <- obj@prin_coords_cols
  } else {
    stop("princ_coords must be either 1 for rows or 2 for columns.")
  }

  rows <- as.data.frame(rows)
  cols <- as.data.frame(cols)
  if (type == "ggplot"){

    # rows <- as.data.frame(rows)
    # cols <- as.data.frame(cols)
    #
    rnmx <- colnames(rows)[xdim]
    rnmy <- colnames(rows)[ydim]
    cnmx <- colnames(cols)[xdim]
    cnmy <- colnames(cols)[ydim]

    p <- ggplot2::ggplot()+
      ggplot2::geom_point(data=rows, ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy)), colour = "#0066FF", alpha = 0.7, shape = 1) +
      ggplot2::geom_point(data=cols, ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy)), colour = "#990000", shape = 4) +
      ggplot2::theme_bw()

    if (!is.null(row_labels)){
      p <- p +
        ggplot2::geom_point(data=rows[row_labels,], ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy)), colour = "#FF0000", shape = 16) +
        ggrepel::geom_text_repel(data=rows[row_labels,], ggplot2::aes_(x = as.name(rnmx), y = as.name(rnmy), label=rownames(rows[row_labels,])), colour = "#FF0000", max.overlaps = Inf)
    }
    if (!is.null(col_labels)){
      p <- p +
        ggplot2::geom_point(data=cols[col_labels,], ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy)), colour = "#990000", shape = 1) +
        ggrepel::geom_text_repel(data=cols[col_labels,], ggplot2::aes_(x = as.name(cnmx), y = as.name(cnmy), label=rownames(cols[col_labels,])), colour = "#990000", max.overlaps = Inf)
    }
  } else if (type == "plotly"){
    p <- plotly::plot_ly(type='scatter', source='plot2D', mode='markers') %>%
      plotly::add_trace(x = cols[,xdim],
                y = cols[,ydim],
                mode = 'markers',
                text = rownames(cols),
                textposition = "left",
                opacity = 1,
                marker = list(color = '#990000', symbol = 'x', size = 5),
                name = 'Columns',
                hoverinfo = 'text',
                type = 'scatter') %>%
      plotly::add_trace(x = rows[,xdim],
                y = rows[,ydim],
                mode = 'markers',
                text = rownames(rows),
                opacity = 0.7,
                marker = list(color ='#0066FF', symbol = 'circle-open', size = 2.5),
                name = 'genes',
                hoverinfo = 'text',
                type = 'scatter')

    if (!is.null(row_labels)){
      p <- p %>%
        plotly::add_trace(x = rows[row_labels, xdim],
                  y = rows[row_labels, ydim],
                  mode = 'markers+text',
                  text = rownames(rows)[row_labels],
                  textposition = "left",
                  textfont = list(color='#FF0000'),
                  marker = list(symbol = 'circle', color ='#FF0000', size = 5),
                  name = 'marked row(s)',
                  hoverinfo = 'text',
                  type = 'scatter')
    }

    if (!is.null(col_labels)){
      p <- p %>%
        plotly::add_trace(x = cols[col_labels, xdim],
                  y = cols[col_labels, ydim],
                  mode = 'markers+text',
                  text = rownames(cols)[col_labels],
                  textposition = "left",
                  textfont = list(color='#990000'),
                  marker = list(symbol = 'circle-open', color ='#990000', size = 6.5),
                  name = 'marked column(s)',
                  hoverinfo = 'text',
                  type = 'scatter')
    }

    p <- p %>%
      plotly::layout(autosize = TRUE,
             title = '2D CA plot',
             showlegend = FALSE,
             xaxis = list(title = paste0('Dim', xdim)),
             yaxis = list(title = paste0('Dim', ydim)))

  }

  return(p)

})


#' @rdname ca_biplot
#' @param assay Seurat assay for recomputation.
#' @param slot Seurat assay slot from which to get matrix.
#' @export
setMethod(f = "ca_biplot",
          signature=(obj="Seurat"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   type = "plotly",
                   ...,
                   assay = Seurat::DefaultAssay(obj),
                   slot = "counts"){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  if ("CA" %in% Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay, slot = slot)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(seurat_obj, assay) first.")
  }

 p <-  ca_biplot(obj = caobj,
            xdim = xdim,
            ydim = ydim,
            princ_coords = princ_coords,
            row_labels = row_labels,
            col_labels = col_labels,
            type = type)

 return(p)
})


#' @rdname ca_biplot
#' @param assay SingleCellExperiment assay for recomputation
#' @export
setMethod(f = "ca_biplot",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   xdim = 1,
                   ydim = 2,
                   princ_coords = 1,
                   row_labels = NULL,
                   col_labels = NULL,
                   type = "plotly",
                   ...,
                   assay = "counts"){

  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))

  if ("CA" %in% reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  } else {
    stop("No 'CA' dimension reduction object found. Please run cacomp(sce, top, coords = FALSE, return_input=TRUE) first.")
  }


  p <-  ca_biplot(obj = caobj,
                  xdim = xdim,
                  ydim = ydim,
                  princ_coords = princ_coords,
                  row_labels = row_labels,
                  col_labels = col_labels,
                  type = type)

  return(p)
})

#' Generates plot for results from apl_topGO
#' @description
#' Plots the results from the data frame generated via apl_topGO.
#'
#' @param genenr data.frame. gene enrichment results table.
#' @param ntop numeric. Number of elements to plot.
#'
#' @return
#' Returns a ggplot plot.
#' @export
#' @examples
#' library(Seurat)
#' set.seed(1234)
#' cnts <- GetAssayData(pbmc_small, slot = "counts")
#' cnts <- as.matrix(cnts)
#'
#' # Run CA on example from Seurat
#'
#' ca <- cacomp(pbmc_small,
#'              princ_coords = 3,
#'              return_input = FALSE,
#'              assay = "RNA",
#'              slot = "counts")
#'
#' grp <- which(Idents(pbmc_small) == 2)
#' ca <- apl_coords(ca, group = grp)
#' ca <- apl_score(ca,
#'                 mat = cnts)
#'
#' enr <- apl_topGO(ca,
#'                  ontology = "BP",
#'                  organism = "hs")
#'
#' plot_enrichment(enr)
plot_enrichment <- function(genenr, ntop = 10){

  genenr$geneRatio <- genenr$Significant/genenr$Annotated
  genenr$raw.p.value[genenr$raw.p.value == "< 1e-30"] <- 0
  genenr$raw.p.value <- as.numeric(genenr$raw.p.value)
  genenr$Term <- factor(genenr$Term , levels = rev(unique(genenr$Term)))

  ggplot2::ggplot(genenr[seq_len(ntop),],
                  aes(x = geneRatio,
                      y = Term,
                      size = Significant,
                      fill = raw.p.value)) +
    ggplot2::geom_point(shape = 21, color = "black") +
    # ggplot2::scale_color_continuous(low="red", high="blue", name = "p-value",
                           # guide=guide_colorbar(reverse=TRUE)) +
    ggplot2::scale_fill_viridis_c(name = "p-value",
                                  option = "viridis",
                                  direction = -1,
                                  guide=guide_colorbar(reverse=TRUE))+
    ggplot2::labs(x = "Significant/Annotated",
                  y = "GO Terms",
                  size = "# Significant") +
    ggplot2::theme_bw()
}
