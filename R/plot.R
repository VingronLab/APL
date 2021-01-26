
#' Plot of the first 3 CA dimensions.
#' 
#' @description 
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and columnnames respectively.
#' @return 
#' Plot of class "plotly".
#' 
#' @param obj An object of class "cacomp", "Seurat" or "SingleCellExperiment", with standard and principal coordinates calculated.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
ca_3Dplot <- function(obj, xdim=1, ydim=2, zdim = 3, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(obj$std_coords_cols))) UseMethod("ca_3Dplot")

#' Plot of the first 3 CA dimensions.
#' 
#' @description 
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and columnnames respectively.
#' @return 
#' Plot of class "plotly".
#' 
#' @param obj Any object.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
ca_3Dplot.default <- function(obj, xdim=1, ydim=2, zdim = 3, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(caobj$std_coords_cols))){
  stop(paste0("ca_3Dplot does not know how to handle objects of class ", 
              class(obj),
              ". Currently only objects of class 'cacomp', 'Seurat' or 'SingleCellExperiment' are supported."))
}
#' Plot of the first 3 CA dimensions.
#' 
#' @description 
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and columnnames respectively.
#' @return 
#' Plot of class "plotly".
#' 
#' @param obj An object of class "cacomp" with standard and principal coordinates calculated.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
ca_3Dplot.cacomp <- function(obj, xdim=1, ydim=2, zdim = 3, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(obj$std_coords_cols))){
  
  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }  
  
  if (princ_coords == 1){
    
    if(!sum(c("prin_coords_rows", "std_coords_cols") %in% names(obj))==2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj$prin_coords_rows
    cols <- obj$std_coords_cols
  } else if (princ_coords == 2){
    if(!sum(c("prin_coords_cols", "std_coords_rows") %in% names(obj))==2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj$std_coords_rows
    cols <- obj$prin_coords_cols
  } else {
    cat("princ_coords must be either 1 for rows or 2 for columns.")
    stop()
  }
  
  p <- plot_ly(type='scatter', source='plot2D', mode='markers') %>%
    add_trace(x = cols[,xdim],
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
    add_trace(x = rows[,xdim],
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
      add_trace(x = rows[row_labels, xdim],
                y = rows[row_labels, ydim],
                z = rows[row_labels, zdim],
                mode = 'markers+text',
                text = rownames(rows)[row_labels],
                textposition = "left",
                textfont = list(color='#0066FF'),
                marker = list(symbol = 'circle', color ='#0066FF', size = 2),
                name = 'marked row(s)',
                hoverinfo = 'text',
                type = 'scatter3d')
  } 
  
  if (!is.null(col_labels)){
    p <- p %>%
      add_trace(x = cols[col_labels, xdim],
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
  
  p <- p %>% 
    layout(autosize = T,
           title = '3D CA plot',
           showlegend = FALSE,
           xaxis = list(title = paste0('Dim',xdim)),
           yaxis = list(title = paste0('Dim', ydim)))
  
  return(p)
  
}


#' Plot of the first 3 CA dimensions.
#' 
#' @description 
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and columnnames respectively.
#' @return 
#' Plot of class "plotly".
#' 
#' @param obj An object of class "Seurat" with standard and principal coordinates calculated and stored in the "CA" DimReduc slot.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(caobj$std_coords_cols)) (all columns).
ca_3Dplot.Seurat <- function(obj, xdim=1, ydim=2, zdim = 3, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(caobj$std_coords_cols))){
  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))
  
  if ("CA" %in% Reductions(obj)){
    caobj <- as.cacomp.Seurat(obj, recompute = FALSE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(seurat_obj, assay) first.")
  }
  
  p <- ca_3Dplot.cacomp(obj = caobj,
                        xdim = xdim,
                        ydim = ydim,
                        zdim = zdim,
                        princ_coords = princ_coords,
                        row_labels = row_labels,
                        col_labels = col_labels)
  return(p)
}



#' Plot of the first 3 CA dimensions.
#' 
#' @description 
#' Plots the first 3 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row- and columnnames respectively.
#' @return 
#' Plot of class "plotly".
#' 
#' @param obj An object of class "SingleCellExperiment" with standard and principal coordinates calculated and stored in the "CA" reducedDim slot.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param zdim Integer. The dimension for the z-axis. Default 3.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(caobj$std_coords_cols)) (all columns).
ca_3Dplot.SingleCellExperiment <- function(obj, xdim=1, ydim=2, zdim = 3, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(caobj$std_coords_cols))){
  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))
  
  if ("CA" %in% reducedDimNames(obj)){
    caobj <- as.cacomp.SingleCellExperiment(obj,  recompute = FALSE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(sce, top, coords = FALSE, return_input=TRUE) first.")
  }
  
  p <- ca_3Dplot.cacomp(obj = caobj,
                        xdim = xdim,
                        ydim = ydim,
                        zdim = zdim,
                        princ_coords = princ_coords,
                        row_labels = row_labels,
                        col_labels = col_labels)
  return(p)
}

#' Plot of the first 2 CA dimensions.
#' 
#' @description 
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Choosing type "plotly" will generate an interactive html plot with the package plotly. Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return 
#' Plot of class "plotly" or "ggplot".
#' 
#' @param obj An object of class "cacomp", "Seurat" or "SingleCellExperiment" with standard and principal coordinates calculated.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
#' @export
ca_biplot <- function(obj, xdim=1, ydim=2, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(obj$std_coords_cols)), type = "plotly") UseMethod("ca_biplot")

#' Plot of the first 2 CA dimensions.
#' 
#' @description 
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Choosing type "plotly" will generate an interactive html plot with the package plotly. Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return 
#' Plot of class "plotly" or "ggplot".
#' 
#' @param obj obj
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
#' @export
ca_biplot.default <- function(obj, xdim=1, ydim=2, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(obj$std_coords_cols)), type = "plotly"){
  stop(paste0("ca_biplot does not know how to handle objects of class ", 
              class(obj),
              ". Currently only objects of class 'cacomp', 'Seurat' or 'SingleCellExperiment' are supported."))
}
#' Plot of the first 2 CA dimensions.
#' 
#' @description 
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Choosing type "plotly" will generate an interactive html plot with the package plotly. Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return 
#' Plot of class "plotly" or "ggplot".
#' 
#' @param obj An object of class "cacomp" with standard and principal coordinates calculated.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
ca_biplot.cacomp <- function(obj, xdim=1, ydim=2, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(obj$std_coords_cols)), type = "plotly"){
  
  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }
  
  if (princ_coords == 1){
    
    if(!sum(c("prin_coords_rows", "std_coords_cols") %in% names(obj))==2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj$prin_coords_rows
    cols <- obj$std_coords_cols
  } else if (princ_coords == 2){
    if(!sum(c("prin_coords_cols", "std_coords_rows") %in% names(obj))==2){
      stop("Principal and/or standard coordinates not found, please run ca_coords() first!")
    }
    rows <- obj$std_coords_rows
    cols <- obj$prin_coords_cols
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
    
    p <- ggplot()+
      geom_point(data=rows, aes_(x = as.name(rnmx), y = as.name(rnmy)), colour = "#0066FF", alpha = 0.7, shape = 1) +
      geom_point(data=cols, aes_(x = as.name(cnmx), y = as.name(cnmy)), colour = "#990000", shape = 4) +
      theme_bw()
    
    if (!is.null(row_labels)){ 
      p <- p + 
        geom_point(data=rows[row_labels,], aes_(x = as.name(rnmx), y = as.name(rnmy)), colour = "#0066FF", shape = 16) + 
        geom_text_repel(data=rows[row_labels,], aes_(x = as.name(rnmx), y = as.name(rnmy), label=rownames(rows[row_labels,])), colour = "#0066FF")
    }
    if (!is.null(col_labels)){
      p <- p + 
        geom_point(data=cols[col_labels,], aes_(x = as.name(cnmx), y = as.name(cnmy)), colour = "#990000", shape = 1) +
        geom_text_repel(data=cols[col_labels,], aes_(x = as.name(cnmx), y = as.name(cnmy), label=rownames(cols[col_labels,])), colour = "#990000")
    }
  } else if (type == "plotly"){
    p <- plot_ly(type='scatter', source='plot2D', mode='markers') %>%
      add_trace(x = cols[,xdim],
                y = cols[,ydim],
                mode = 'markers', 
                text = rownames(cols),
                textposition = "left",
                opacity = 1,
                marker = list(color = '#990000', symbol = 'x', size = 5),
                name = 'Columns',
                hoverinfo = 'text',
                type = 'scatter') %>%
      add_trace(x = rows[,xdim],
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
        add_trace(x = rows[row_labels, xdim],
                  y = rows[row_labels, ydim],
                  mode = 'markers+text',
                  text = rownames(rows)[row_labels],
                  textposition = "left",
                  textfont = list(color='#0066FF'),
                  marker = list(symbol = 'circle', color ='#0066FF', size = 5),
                  name = 'marked row(s)',
                  hoverinfo = 'text',
                  type = 'scatter')
    } 
    
    if (!is.null(col_labels)){
      p <- p %>%
        add_trace(x = cols[col_labels, xdim],
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
      layout(autosize = T,
             title = '2D CA plot',
             showlegend = FALSE,
             xaxis = list(title = paste0('Dim',xdim)),
             yaxis = list(title = paste0('Dim', ydim)))
    
  }
  
  return(p)
  
}


#' Plot of the first 2 CA dimensions.
#' 
#' @description 
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Choosing type "plotly" will generate an interactive html plot with the package plotly. Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return 
#' Plot of class "plotly" or "ggplot".
#' 
#' @param obj An object of class "Seurat" with standard and principal coordinates calculated and stored in the dim. Reductions slot.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
#' @export
ca_biplot.Seurat <- function(obj, xdim=1, ydim=2, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(caobj$std_coords_cols)), type = "plotly"){
  
  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))
  
  if ("CA" %in% Reductions(obj)){
    caobj <- as.cacomp.Seurat(obj, recompute = FALSE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(seurat_obj, assay) first.")
  }
  
 p <-  ca_biplot.cacomp(obj = caobj,
            xdim = xdim,
            ydim = ydim,
            princ_coords = princ_coords,
            row_labels = row_labels,
            col_labels = col_labels,
            type = type)
 
 return(p)
}


#' Plot of the first 2 CA dimensions.
#' 
#' @description 
#' Plots the first 2 dimensions of the rows and columns in the same plot.
#' 
#' @details 
#' Choosing type "plotly" will generate an interactive html plot with the package plotly. Type "ggplot" generates a static plot.
#' Depending on whether `princ_coords` is set to 1 or 2 either the principal coordinates of either the rows (1) or the columns (2)
#' are chosen. For the other the standard coordinates are plotted (assymetric biplot).
#' Labels for rows and columns should be stored in the row and column names respectively.
#' @return 
#' Plot of class "plotly" or "ggplot".
#' 
#' @param obj AAn object of class "SingleCellExperiment" with standard and principal coordinates calculated and stored in the dim. Reductions slot.
#' @param xdim Integer. The dimension for the x-axis. Default 1.
#' @param ydim Integer. The dimension for the y-axis. Default 2.
#' @param princ_coords Integer. If 1 then principal coordinates are used for the rows, if 2 for the columns. Default 1 (rows).
#' @param row_labels Numeric vector. Indices for the rows for which a label should be added (label should be stored in rownames). Default NULL.
#' @param col_labels Numeric vector. Indices for the columns for which a label should be added (label should be stored in colnames). 
#' Default seq(ncol(obj$std_coords_cols)) (all columns).
#' @param type String. Type of plot to draw. Either "ggplot" or "plotly". Default "plotly".
#' @export
ca_biplot.SingleCellExperiment <- function(obj, xdim=1, ydim=2, princ_coords = 1, row_labels=NULL, col_labels=seq(ncol(caobj$std_coords_cols)), type = "plotly"){
  
  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))
  
  if ("CA" %in% reducedDimNames(obj)){
    caobj <- as.cacomp.SingleCellExperiment(obj,  recompute = FALSE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(sce, top, coords = FALSE, return_input=TRUE) first.")
  }
  
  
  p <-  ca_biplot.cacomp(obj = caobj,
                  xdim = xdim,
                  ydim = ydim,
                  princ_coords = princ_coords,
                  row_labels = row_labels,
                  col_labels = col_labels,
                  type = type)
  
  return(p)
}
  