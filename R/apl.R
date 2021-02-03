

#' Calculate Association Plot coordinates
#'
#' @description
#' Calculates the association plot coordinates for either the rows, columns or both (default).
#'
#' @details
#' Coordinates (x,y) of row vector \eqn{\vec{r}} are defined as
#' \deqn{x(\vec{r}) := \left|\vec{r}\right|\cos(\phi(\vec{r}))}
#' \deqn{y(\vec{r}) := \left|\vec{r}\right|\sin(\phi(\vec{r}))}
#' The x-direction is determined by calculating the centroid of the columns selected with the indices in "group".
#'
#' @return
#' Returns input "cacomp" object and adds components "apl_rows" and/or "apl_cols" for row and column coordinates.
#' In "group" the indices of the columns used to calculate the centroid are saved.
#'
#' @param caobj A "cacomp" object with principal row coordinates and standard column coordinates calculated.
#' @param group Numeric/Character. Vector of indices or column names of the columns to calculate centroid/x-axis direction.
#' @param calc_rows TRUE/FALSE. Whether apl row coordinates should be calculated. Default TRUE.
#' @param calc_cols TRUE/FALSE. Whether apl column coordinates should be calculated. Default TRUE.
#' @export
apl_coords <- function(caobj, group, calc_rows = TRUE, calc_cols = TRUE){

  stopifnot(is(caobj, "cacomp"))

  rows <- t(caobj$prin_coords_rows)
  cols <- t(caobj$std_coords_cols)

  if (is(group, "numeric")){
    subgroup <- cols[,group]
  } else if (is(group, "character")){
    idx <- match(group, colnames(cols))
    idx <- na.omit(idx)
    subgroup <- cols[,idx]

    if (anyNA(idx)){
      warning("Not all names in 'group' are contained in the column names. Non-matching values were ignored.")
    }
  } else {
    stop("Parameter group hast to be either of type 'numeric' or 'character'.")
  }

  avg_group_coords <- rowMeans(subgroup) # centroid vector.
  length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
  length_vector_rows <- sqrt(colSums(rows^2))
  length_vector_cols <- sqrt(colSums(cols^2))

  if (calc_rows == TRUE){
    message("Calculating APL row coordinates ...")
    # r⋅X = |r|*|X|*cosθ
    # x(r) = (r⋅X)/|X| = |r|*cosθ
    rowx <- drop(t(rows) %*% avg_group_coords)/length_vector_group
    # pythagoras, y(r)=b²=c²-a²
    rowy <- sqrt(length_vector_rows^2 - rowx^2)

    caobj$apl_rows <- cbind("x"=rowx, "y"=rowy)
  }


  if (calc_cols == TRUE){
    message("Calculating APL column coordinates ...")

    colx <- drop(t(cols) %*% avg_group_coords)/length_vector_group
    coly <- sqrt(length_vector_cols^2 - colx^2)

    caobj$apl_cols <- cbind("x"=colx, "y"=coly)
  }

  caobj$group <- group

  return(caobj)
}

#' Find rows most highly associated with a condition
#'
#' @description
#' Ranks rows by a calculated score which balances the association of the row with the condition and how associated it is with other conditions.
#'
#' @details
#' The score is calculated by permuting the values of each row to determine the cutoff angle of the 99% quantile.
#' \deqn{S_{alpha}(x,y)=x-\frac{y}{\tan\alpha}}
#' By default the permutation is repeated 10 times, but for very large matrices this can be reduced.
#'
#' @return
#' Returns the input "cacomp" object with "APL_score" component added.
#' APL_score contains a data frame with ranked rows, their score and their original row number.
#'
#' @param caobj A "cacomp" object with principal row coordinates and standard column coordinates calculated.
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param dims Integer. Number of CA dimensions to retain. Needs to be the same as in caobj!
#' @param group Vector of indices of the columns to calculate centroid/x-axis direction.
#' @param reps Integer. Number of permutations to perform. Default = 10.
#' @param quant Numeric. Single number between 0 and 1 indicating the quantile used to calculate the cutoff. Default 0.99.
#' @param python A logical value indicating whether to use singular-value decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @export
apl_score <- function(caobj, mat, dims, group, reps=10, quant = 0.99, python = TRUE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    stop("Please run apl_coords() first!")
  }

  row_num <- nrow(caobj$prin_coords_rows)
  apl_rows_perm <- data.frame("x" = rep(0, row_num*reps), "y" = rep(0, row_num*reps)) #init. data frame

  for (k in seq(reps)){
    message("\nRunning permutation ", k, " out of ", reps, " to calculate row scores ...")

    #permute rows and rerun cacomp
    mat_perm <- t(apply(mat, 1, FUN=sample))
    colnames(mat_perm) <- colnames(mat)

    caobjp <- cacomp(obj = mat_perm,
                     python = python,
                     coords = TRUE,
                     princ_coords = 1,
                     dims = dims,
                     top = caobj$top_rows,
                     inertia = FALSE)

    caobjp <- apl_coords(caobj = caobjp, group = group, calc_cols = FALSE)
    idx <- ((1:row_num)+((k-1)*row_num))
    apl_rows_perm[idx,] <- caobjp$apl_rows
  }

  apl_rows_perm[,3] <- apl_rows_perm[,1]/apl_rows_perm[,2] # cotan between row and x axis
  apl_rows_perm[,3][is.na(apl_rows_perm[,3])] <- 0

  # cutoff from original code from Ela
  # angles_vector <- sort(apl_rows_perm[,3], decreasing = TRUE)
  # cutoff_cotan <- angles_vector[ceiling(0.01 * length(angles_vector))]

  # With 99% quantile, gives different results though!
  cutoff_cotan <- quantile(apl_rows_perm[,3], quant)

  score <- caobj$apl_rows[,1] - (caobj$apl_rows[,2] * cutoff_cotan)
  ranking <- data.frame("Rowname" = rownames(caobj$apl_rows),
                        "Score" = score,
                        "Row_num" = 1:nrow(caobj$apl_rows))

  ranking <- ranking[order(ranking$Score, decreasing = TRUE),]
  ranking$Rank <- 1:nrow(ranking)

  caobj$APL_score <- ranking

  return(caobj)

}

#' Association Plot
#'
#' @description
#' Plot an Association plot for the chosen columns.
#'
#' @details
#' For an interactive plot type="plotly" can be chosen, otherwise a static plot will returned.
#' The row and column coordinates have to be already calculated by `apl_coords()`.
#'
#' @return
#' Either a ggplot or plotly object.
#'
#' @param caobj  An object of class "cacomp" and "APL" with apl coordinates calculated.
#' @param type "ggplot"/"plotly". For a static plot a string "ggplot", for an interactive plot "plotly". Default "ggplot".
#' @param rowlabels TRUE/FALSE. Logical whether rows should be labelled. Default TRUE.
#' @param collabels TRUE/FALSE. Logical whether columns should be labelled. Default TRUE.
#' @param rows_idx numeric vector. Indices of the rows that should be labelled. Default all rows: 1:nrow(caobj$apl_rows).
#' @param cols_idx numeric vector. Indices of the columns that should be labelled. Default is only to label columns making up the centroid: caobj$group.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx shouls be labeled with text. Default TRUE.
#' @export
apl <- function(caobj, type="ggplot", rowlabels = TRUE, collabels = TRUE, rows_idx = 1:nrow(caobj$apl_rows), cols_idx = caobj$group, row_labs = TRUE, col_labs = TRUE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    stop("Please run apl_coords() first!")
  }

  group_cols <- data.frame(rownms = rownames(caobj$apl_cols)[cols_idx], x = caobj$apl_cols[cols_idx,"x"], y = caobj$apl_cols[cols_idx,"y"], row.names = NULL)
  group_rows <- data.frame(rownms = rownames(caobj$apl_rows)[rows_idx], x = caobj$apl_rows[rows_idx,"x"], y = caobj$apl_rows[rows_idx,"y"], row.names = NULL)

  if (row_labs == FALSE){
    rlabs <- 'markers'
    rowfont <- NULL
  } else {
    rlabs <- 'markers+text'
    rowfont <- list(color='#FF0000')

  }

  if (col_labs == FALSE){
    clabs <- 'markers'
    colfont <- NULL
  } else {
    clabs <- 'markers+text'
    colfont <- list(color='#000000')

  }

  apl_rows.tmp <- data.frame(rownms = rownames(caobj$apl_rows), caobj$apl_rows)
  apl_cols.tmp <- data.frame(rownms = rownames(caobj$apl_cols), caobj$apl_cols)

  if (type == "ggplot"){


    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data=apl_cols.tmp, ggplot2::aes(x=x, y=y), color = "#990000", shape = 4) +
      ggplot2::geom_point(data=apl_cols.tmp[caobj$group,], ggplot2::aes(x=x, y=y), color = "#990000", shape = 1)+
      ggplot2::geom_point(data=apl_rows.tmp, ggplot2::aes(x=x, y=y), color = "#0066FF", alpha = 0.5, shape = 1) +
      ggplot2::labs(title="Association Plot") +
      ggplot2::theme_bw()

    if(collabels == TRUE){
      p <- p +
        ggrepel::geom_text_repel(data=group_cols, ggplot2::aes(x=x, y=y, label=rownms), color = "#990000")}
    if (rowlabels == TRUE){
      p <- p +
        ggplot2::geom_point(data=group_rows, ggplot2::aes(x=x, y=y), color="#FF0000", shape = 16) +
        ggrepel::geom_text_repel(data = group_rows, ggplot2::aes(x=x, y=y, label=rownms), color = "#FF0000")}
    rm(apl_rows.tmp, apl_cols.tmp)

    return(p)

  } else if (type == "plotly"){
    p <- plotly::plot_ly() %>%
      plotly::add_trace(data=apl_cols.tmp,
                x = ~x,
                y =  ~y,
                mode = 'markers',
                text = apl_cols.tmp$rownms,
                textposition = "left",
                marker = list(color = '#124429',
                              symbol = 'x',
                              size = 5),
                name = 'samples',
                hoverinfo = 'text',
                type = 'scatter') %>%
      plotly::add_trace(data = apl_rows.tmp,
                x = ~x,
                y = ~y,
                mode = 'markers',
                text = apl_rows.tmp$rownms,
                opacity = 0.7,
                textposition = "left",
                marker = list(color ='#0066FF',
                              symbol = 'circle-open',
                              size = 5),
                name = 'genes',
                hoverinfo = 'text',
                type = 'scatter') %>%
      plotly::add_trace(data = group_rows,
                x = ~x,
                y = ~y,
                mode = rlabs,
                text = group_rows$rownms,
                textposition = "left",
                textfont=rowfont,
                marker = list(symbol = 'circle',
                              color = '#FF0000',
                              size = 5),
                name = 'marked genes',
                hoverinfo = 'text',
                type = 'scatter') %>%
      plotly::add_trace(data = group_cols,
                x = ~x,
                y = ~y,
                mode = clabs,
                text = group_cols$rownms,
                textposition = "left",
                textfont=colfont,
                marker = list(symbol = 'x',
                              color = '#990000',
                              size = 5),
                name = 'marked samples',
                hoverinfo = 'text',
                type = 'scatter') %>%
      plotly::layout(title = paste('Association Plot \n', ncol(caobj$U), ' first dimensions, ', length(caobj$group),' samples.\n'),
             xaxis = list(title = 'Distance from origin (x)', rangemode = "tozero"),
             yaxis = list(title = 'Distance from gene to sample line (y)', rangemode = "tozero"),showlegend = TRUE)

    rm(apl_rows.tmp, apl_cols.tmp)

    return(p)
  } else {
    stop("Please specify plot = \"ggplot\" or \"plotly\". Other options are not accepted.")
  }

}


#' Compute and plot Association plot
#'
#' @description
#' Computes singular value decomposition and coordinates for the Association plot.
#'
#' @details
#' The function is a wrapper that calls `cacomp()`, `apl_coords()`, `apl_score()` and finally `apl()` for ease of use.
#' The chosen defaults are most useful for genomics experiments, but for more fine grained control the functions
#' can be also run individually for the same results.
#' If score = FALSE, nrow and reps are ignored. If mark_rows is not NULL score is treated as if FALSE.
#' @return
#' Association plot (plotly object).
#'
#' @param obj A numeric matrix, Seurat or SingleCellExperiment object. For sequencing a count matrix, gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param caobj A "cacomp" object as outputted from `cacomp()`. If not supplied will be calculated. Default NULL.
#' @param dims Integer. Number of dimensions to keep. Default NULL (keeps all dimensions).
#' @param group Numeric/Character. Vector of indices or column names of the columns to calculate centroid/x-axis direction.
#' @param nrow Integer. The top nrow scored row labels will be added to the plot if score = TRUE. Default 10.
#' @param top Integer. Number of most variable rows to retain. Default 5000.
#' @param score Logical. Whether rows should be scored and ranked. Ignored when a vector is supplied to mark_rows. Default TRUE.
#' @param mark_rows Character vector. Names of rows that should be highlighted in the plot. If not NULL, score is ignored. Default NULL.
#' @param reps Integer. Number of permutations during scoring. Default 3.
#' @param python A logical value indicating whether to use singular-value decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx shouls be labeled with text. Default TRUE.
#' @param ... Arguments forwarded to methods.
#' @export
runAPL <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, ...){
  UseMethod("runAPL")
}



#' @rdname runAPL
#' @export
runAPL.default <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, ...){
  stop(paste0("runAPL does not know how to handle objects of class ",
              class(x),
              ". Currently only objects of class 'matrix' or objects coercible to one, 'Seurat' or 'SingleCellExperiment' are supported."))
}



#' @export
#' @rdname runAPL
runAPL.matrix <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, ...){

  if (!is(obj, "matrix")){
    obj <- as.matrix(obj)
  }
  stopifnot(is(obj, "matrix"))

  if (is.null(caobj)){
    caobj <- cacomp(obj = obj,
                    coords = TRUE,
                    top = top,
                    princ_coords = 1,
                    dims = dims,
                    python = TRUE)

  } else {
    if(!is.null(caobj$dims) && is.null(dims)){
      dims <- caobj$dims
    } else if (!is.null(caobj$dims) && !is.null(dims)) {
        # warning("The caobj was previously already subsetted to ", caobj$dims, " dimensions. Subsetting again!")
      if (dims < caobj$dims){
        # caobj <- ca_coords(caobj = caobj, dims = dims, princ_only = FALSE, princ_coords = 1)
        caobj <- subset_dims(caobj, dims = dims)
      } else if(dims > length(caobj$D)){
        warning("dims is larger than the number of available dimensions. Argument ignored")
      }
    }

    if (is.null(caobj$prin_coords_rows) && !is.null(caobj$std_coords_rows)){
      caobj <- ca_coords(caobj = caobj, dims = dims, princ_only = TRUE, princ_coords = 1)
    } else if (is.null(caobj$prin_coords_rows) || is.null(caobj$std_coords_cols)){
      caobj <- ca_coords(caobj = caobj, dims = dims, princ_only = FALSE, princ_coords = 1)
    }
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    caobj <- apl_coords(caobj = caobj, group = group)
  }


  if (is.null(mark_rows)){
    if (score == TRUE){
      caobj <- apl_score(caobj = caobj,
                         mat = obj,
                         dims = caobj$dims,
                         group = caobj$group,
                         reps= reps)

      mark_rows <- head(caobj$APL_score$Row_num, nrow)
    } else {
      mark_rows <- NULL
    }
  } else{
    if (is(mark_rows, "character")){
      mark_rows <- match(mark_rows, rownames(caobj$apl_rows))
      mark_rows <- na.omit(mark_rows)
      if (anyNA(mark_rows)){
        warning("Not all names in 'mark_rows' are contained in the row names. Maybe they were filtered out?\n Non-matching values were ignored.")
      }
    } else  {
      stop("Parameter mark_rows hast to be of type 'character'.")
    }
  }

  p <- apl(caobj = caobj,
           type = "plotly",
           rowlabels = TRUE,
           collabels = TRUE,
           rows_idx = mark_rows,
           cols_idx = caobj$group,
           row_labs = row_labs,
           col_labs = col_labs)

  return(p)
}



#' @description
#' runAPL.SingleCellExperiment: Computes singular value decomposition and coordinates for the Association plot from SingleCellExperiment objects with reducedDim(obj, "CA") slot (optional).
#'
#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#'
#' @rdname runAPL
#' @export
runAPL.SingleCellExperiment <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, ..., assay = "counts"){

  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))

  mat <- SummarizedExperiment::assay(obj, assay)

  if ("CA" %in% SingleCellExperiment::reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay, recompute = TRUE)
  }

  runAPL.matrix(obj = mat,
                caobj = caobj,
                dims = dims,
                group = group,
                mark_rows = mark_rows,
                nrow = nrow,
                top = top,
                score = score,
                reps = reps,
                python = python,
                row_labs = row_labs,
                col_labs = col_labs)

}


#' @description
#' runAPL.Seurat: Computes singular value decomposition and coordinates for the Association plot from Seurat objects, optionally with a DimReduc Object in the "CA" slot.
#'
#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#'
#' @rdname runAPL
#' @export
runAPL.Seurat <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, ..., assay = DefaultAssay(obj)){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = "data")

  if ("CA" %in% Seurat::Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay, recompute = TRUE)
  }

  runAPL.matrix(obj = seu,
                caobj = caobj,
                top = top,
                dims= dims,
                group = group,
                score = score,
                reps = reps,
                python = python,
                mark_rows = mark_rows,
                nrow = nrow,
                row_labs = row_labs,
                col_labs = col_labs)
}

