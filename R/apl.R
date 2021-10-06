

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
  cent <- cols


  if (is(group, "numeric")){
    subgroup <- cent[,group]
  } else if (is(group, "character")){
    idx <- match(group, colnames(cent))
    idx <- na.omit(idx)
    subgroup <- cent[,idx]

    if (anyNA(idx)){
      warning("Not all names in 'group' are contained in the column names. Non-matching values were ignored.")
    }
  } else {
    stop("Parameter group hast to be either of type 'numeric' or 'character'.")
  }

  if (length(group) == 1){
    avg_group_coords <- subgroup # single sample
  } else {
    avg_group_coords <- rowMeans(subgroup) # centroid vector.
  }
  length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
  length_vector_rows <- sqrt(colSums(rows^2))
  length_vector_cols <- sqrt(colSums(cols^2))

  if (calc_rows == TRUE){
    # message("Calculating APL row coordinates ...")
    # r⋅X = |r|*|X|*cosθ
    # x(r) = (r⋅X)/|X| = |r|*cosθ
    rowx <- drop(t(rows) %*% avg_group_coords)/length_vector_group
    # pythagoras, y(r)=b²=c²-a²
    rowy <- sqrt(length_vector_rows^2 - rowx^2)

    rowx[is.na(rowx)] <- 0
    rowy[is.na(rowy)] <- 0
    # rowx[is.infinite(rowx)] <- 0
    # rowy[is.infinite(rowy)] <- 0

    caobj$apl_rows <- cbind("x"=rowx, "y"=rowy)
  }


  if (calc_cols == TRUE){
    # message("Calculating APL column coordinates ...")

    colx <- drop(t(cols) %*% avg_group_coords)/length_vector_group
    coly <- sqrt(length_vector_cols^2 - colx^2)

    colx[is.na(colx)] <- 0
    coly[is.na(coly)] <- 0
    # colx[is.infinite(colx)] <- 0
    # coly[is.infinite(coly)] <- 0

    caobj$apl_cols <- cbind("x"=colx, "y"=coly)
  }

  if (is(group, "numeric")){
    caobj$group <- group
  } else if (is(group, "character")){
    idx <- match(group, colnames(cols))
    idx <- na.omit(idx)
    caobj$group <- idx
  }

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
#' If store_perm is TRUE the permuted data is stored in the cacomp object and can be used for future scoring.
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
#' @param store_perm Logical. Whether permuted data should be stored in the CA object.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @export
apl_score <- function(caobj, mat, dims, group, reps=10, quant = 0.99, python = TRUE, store_perm = TRUE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    stop("Please run apl_coords() first!")
  }

  names <- colnames(mat)
  row_num <- nrow(caobj$apl_rows)
  margin <- 1
  pc <- 1
  cc <- FALSE
  cr <- TRUE

  if (is(group, "character")){
    idx <- match(group, names)
    idx <- na.omit(idx)
    group <- idx
  }

  if (caobj$dims == 1 && !is.null(caobj$dims)){
    row_num <- 1
  }

  apl_perm <- data.frame("x" = rep(0, row_num*reps), "y" = rep(0, row_num*reps)) #init. data frame
  saved_ca <- list()
  pb <- txtProgressBar(min = 0, max = reps, style = 3)

  for (k in seq(reps)){
    # message("\nRunning permutation ", k, " out of ", reps, " to calculate row/column scores ...")

    #permute rows and rerun cacomp

    if(isTRUE(store_perm) & identical(reps, attr(caobj$permuted_data,'reps'))){

      caobjp <- caobj$permuted_data[[k]]

    } else {
      mat_perm <- t(apply(mat, margin, FUN=sample))
      colnames(mat_perm) <- colnames(mat)


      suppressWarnings(caobjp <- cacomp(obj = mat_perm,
                                         python = python,
                                         coords = TRUE,
                                         princ_coords = pc,
                                         dims = dims,
                                         top = caobj$top_rows,
                                         inertia = FALSE))

      if(isTRUE(store_perm)){
        x <- list("std_coords_cols" = caobjp$std_coords_cols,
                  "prin_coords_rows" = caobjp$prin_coords_rows)

        saved_ca[[k]] <- new_cacomp(x)
      }
    }




    caobjp <- apl_coords(caobj = caobjp, group = group, calc_cols = cc, calc_rows = cr)
    idx <- ((1:row_num)+((k-1)*row_num))

    apl_perm[idx,] <- caobjp$apl_rows

    setTxtProgressBar(pb, k)

  }

  close(pb)

  apl_perm[,3] <- apl_perm[,1]/apl_perm[,2] # cotan between row and x axis
  apl_perm[,3][is.na(apl_perm[,3])] <- 0

  # cutoff from original code from Ela
  # angles_vector <- sort(apl_perm[,3], decreasing = TRUE)
  # cutoff_cotan <- angles_vector[ceiling(0.01 * length(angles_vector))]

  # With 99% quantile, gives different results though!
  cutoff_cotan <- quantile(apl_perm[,3], quant)

  score <- caobj$apl_rows[,1] - (caobj$apl_rows[,2] * cutoff_cotan)
  ranking <- data.frame("Rowname" = rownames(caobj$apl_rows),
                        "Score" = score,
                        "Row_num" = 1:nrow(caobj$apl_rows))

  ranking <- ranking[order(ranking$Score, decreasing = TRUE),]
  ranking$Rank <- 1:nrow(ranking)

  caobj$APL_score <- ranking

  if(isTRUE(store_perm) & !identical(reps, attr(caobj$permuted_data,'reps'))){
    caobj$permuted_data <- saved_ca
    attr(caobj$permuted_data,'cutoff') <- cutoff_cotan
    attr(caobj$permuted_data,'reps') <- reps
  }

  return(caobj)

}


#' Find rows most highly associated with a condition through random APLs.
#'
#' @description
#' Ranks rows by a calculated score which balances the association of the row with the condition and how associated it is with other conditions.
#'
#' @details
#' The score is calculated by choosing random directions in space to calculate APLs for the rows.
#' \deqn{S_{alpha}(x,y)=x-\frac{y}{\tan\alpha}}
#' By default the permutation is repeated 300 times, and is independent of group size, so the same cutoff is applicable to all groups.
#' If store_perm is TRUE the calculated cutoff is stored as an attribute to ca$permuted_data to prevent recalculation when running with identical parameters.
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
#' @param store_perm Logical. Whether calculated cutoff should be stored. Default TRUE.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @export
apl_score_rand <- function(caobj, dims, reps=300, quant = 0.99, python = TRUE, store_perm = TRUE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    stop("Please run apl_coords() first!")
  }

  row_num <- nrow(caobj$apl_rows)
  margin <- 1
  pc <- 1
  cc <- FALSE
  cr <- TRUE

  if (caobj$dims == 1 && !is.null(caobj$dims)){
    row_num <- 1
  }

  apl_perm <- data.frame("x" = rep(0, row_num*reps), "y" = rep(0, row_num*reps)) #init. data frame

  if(isTRUE(store_perm) & identical(reps, attr(caobj$permuted_data,'reps'))){

    cutoff_cotan <- attr(caobj$permuted_data,'cutoff')


  } else {
    pb <- txtProgressBar(min = 0, max = reps, style = 3)

    rows <- t(caobj$prin_coords_rows)
    cols <- t(caobj$std_coords_cols)

    for (k in seq(reps)){

      # avg_group_coords <- rowMeans(subgroup) # centroid vector.
      avg_group_coords <- runif(n=dims, min=0, max = quantile(cols, 0.99))
      length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
      length_vector_rows <- sqrt(colSums(rows^2))

      rowx <- drop(t(rows) %*% avg_group_coords)/length_vector_group
      # pythagoras, y(r)=b²=c²-a²
      rowy <- sqrt(length_vector_rows^2 - rowx^2)

      rowx[is.na(rowx)] <- 0
      rowy[is.na(rowy)] <- 0

      idx <- ((1:row_num)+((k-1)*row_num))
      apl_perm[idx,] <- cbind("x"=rowx, "y"=rowy)

      setTxtProgressBar(pb, k)

    }

    close(pb)

    apl_perm[,3] <- apl_perm[,1]/apl_perm[,2] # cotan between row and x axis
    apl_perm[,3][is.na(apl_perm[,3])] <- 0

    # cutoff from original code from Ela
    # angles_vector <- sort(apl_perm[,3], decreasing = TRUE)
    # cutoff_cotan <- angles_vector[ceiling(0.01 * length(angles_vector))]

    # With 99% quantile, gives different results though!
    cutoff_cotan <- quantile(apl_perm[,3], quant)

  }

  score <- caobj$apl_rows[,1] - (caobj$apl_rows[,2] * cutoff_cotan)
  ranking <- data.frame("Rowname" = rownames(caobj$apl_rows),
                        "Score" = score,
                        "Row_num" = 1:nrow(caobj$apl_rows))

  ranking <- ranking[order(ranking$Score, decreasing = TRUE),]
  ranking$Rank <- 1:nrow(ranking)

  caobj$APL_score <- ranking

  if(isTRUE(store_perm) & !identical(reps, attr(caobj$permuted_data,'reps'))){

    if(is.null(caobj$permuted_data))  caobj$permuted_data <- list()

    attr(caobj$permuted_data,'cutoff') <- cutoff_cotan
    attr(caobj$permuted_data,'reps') <- reps
  }

  return(caobj)

}

apl_topGO <- function(caobj, ontology, organism="hs", ngenes = 1000, top_res = 10, use_coords = FALSE, return_plot = TRUE){

  if (ngenes > nrow(caobj$apl_rows)){
    stop("ngenes is larger than the total number of genes.")
  } else if (ngenes == nrow(caobj$apl_rows)){
    warning("You have selected all available genes. Gene enrichment results might not be meaningful.")
  }

  if(!is(gene_sets, "list")){
    stop("gene_sets should be a list of gene sets, each gene set containing the gene symbols.")
  }


  if(!is.null(caobj$APL_score) & !isTRUE(use_coords)){


    Score_ord <- caobj$APL_score[order(caobj$APL_score$Score),]
    ranked_genes <- Score_ord$Score
    names(ranked_genes) <- Score_ord$Rowname

  } else if (isTRUE(use_coords) | is.null(caobj$APL_score)) {

    if(is.null(caobj$apl_rows)){
      stop("No APL coordinates found for rows. Please first run apl_coords.")
    }

    ranked_genes <- 1:nrow(apl_rows)
    names(ranked_genes) <- rownames(caobj$apl_rows)[order(caobj$apl_rows[,2]),]

  } else {
    stop("APL scores not present but use_coords set to FALSE")
  }

  if (organism == "hs"){
    organism <- "org.Hs.eg.db"
  } else if (organism == "mm"){
    organism <- "org.Mm.eg.db"
  } else {
    Warning("Custom organism chosen.")
    organism <- organism
  }

  if (!ontology %in% c("BP", "CP", "MF")){
    stop("Please choose one of the following ontologies: 'BP', 'CP' or 'MF'.")
  }

  rankedGenes <- ranked_genes[1:ngenes]

  gene_sets <- topGO::annFUN.org(whichOnto=ontology, feasibleGenes=NULL, mapping= organism, ID="symbol")

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = rankedGenes,
                annot = annFUN.GO2genes,
                GO2genes = gene_sets,
                geneSelectionFun = \(x) TRUE,
                nodeSize=5)


  results.ks <- topGO::runTest(GOdata, algorithm="classic", statistic="ks")

  if(top_res > length(results.ks@score)){
    warning("More top nodes selected via top_res than available. Returning max. number of nodes instead.")
  }

  top_res <- min(top_res, length(results.ks@score))
  goEnrichment <- topGO::GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes = top_res)
  # goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]

  if (isTRUE(return_plot)){
    showSigOfNodes(GOdata, score(results.ks), firstSigNodes = top_res, useInfo = 'def')
  }

  return(goEnrichment)

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
#' @param rows_idx numeric/character vector. Indices or names of the rows that should be labelled. Default NULL.
#' @param cols_idx numeric/character vector. Indices or names of the columns that should be labelled. Default is only to label columns making up the centroid: caobj$group.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx shouls be labeled with text. Default TRUE.
#' @param show_score Logical. Wheter the S-alpha score should be shown in the plot.
#' @export
apl <- function(caobj, type="ggplot", rows_idx = NULL, cols_idx = caobj$group, row_labs = FALSE, col_labs = TRUE, show_score = FALSE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.null(caobj$apl_rows) || is.null(caobj$apl_cols)){
    stop("Please run apl_coords() first!")
  }

  if (is(rows_idx, "character")){
    names <- rownames(caobj$apl_rows)
    idx <- match(rows_idx, names)
    idx <- na.omit(idx)
    rows_idx <- idx
  }

  if (is(cols_idx, "character")){
    names <- rownames(caobj$apl_cols)
    idx <- match(cols_idx, names)
    idx <- na.omit(idx)
    cols_idx <- idx
  }

  if (is.numeric(cols_idx)){
  group_cols <- data.frame(rownms = rownames(caobj$apl_cols)[cols_idx], x = caobj$apl_cols[cols_idx,"x"], y = caobj$apl_cols[cols_idx,"y"], row.names = NULL)
  }
  if(is.numeric(rows_idx)){
  group_rows <- data.frame(rownms = rownames(caobj$apl_rows)[rows_idx], x = caobj$apl_rows[rows_idx,"x"], y = caobj$apl_rows[rows_idx,"y"], row.names = NULL)
  }

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
  apl_scores <- caobj$APL_score$Score[order(caobj$APL_score$Row_num)]
  apl_rows.tmp <- data.frame(rownms = rownames(caobj$apl_rows), caobj$apl_rows, Score = apl_scores)
  apl_cols.tmp <- data.frame(rownms = rownames(caobj$apl_cols), caobj$apl_cols)

  if (type == "ggplot"){


    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data=apl_cols.tmp, ggplot2::aes(x=x, y=y), color = "#006400", shape = 4)


      if (isTRUE(show_score)){
        p <- p + ggplot2::geom_point(data=apl_rows.tmp, ggplot2::aes(x=x, y=y, color = Score), alpha = 0.7, shape = 16) +
                # scico::  scale_fill_scico(palette = 'batlow')
                ggplot2::scale_color_viridis_c(option = "D")
      } else {
        p <- p + ggplot2::geom_point(data=apl_rows.tmp, ggplot2::aes(x=x, y=y), color = "#0066FF", alpha = 0.7, shape = 16)
      }
      p <- p +  ggplot2::geom_point(data=apl_cols.tmp[caobj$group,], ggplot2::aes(x=x, y=y), color = "#990000", shape = 4) +
                ggplot2::labs(title="Association Plot") +
                ggplot2::theme_bw()

    if(col_labs == TRUE){
      p <- p +
        ggrepel::geom_text_repel(data=group_cols, ggplot2::aes(x=x, y=y, label=rownms), color = "#990000")
      }
    if (is.numeric(rows_idx)){
      p <- p +
        ggplot2::geom_point(data=group_rows, ggplot2::aes(x=x, y=y), color="#FF0000", shape = 16)

      if(isTRUE(row_labs)){
        p <- p +
          ggrepel::geom_text_repel(data = group_rows, ggplot2::aes(x=x, y=y, label=rownms), color = "#FF0000", max.overlaps = Inf)
      }
    }
    rm(apl_rows.tmp, apl_cols.tmp)

    return(p)

  } else if (type == "plotly"){

    if (isTRUE(show_score)){
      colors_fun <- 'Viridis' #"YlGnBu"
      color_fix <- as.formula("~Score")
      sym <- "cyrcle"
    } else {
      color_fix <- '#0066FF'
      colors_fun <- NULL
      sym <- "cyrcle-open"

    }

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
                marker = list(color = color_fix, # '#0066FF'
                              colorscale = colors_fun,
                              symbol = sym,
                              colorbar=list(title = "Score", len=0.5),
                              size = 5),
                name = 'genes',
                hoverinfo = 'text',
                type = 'scatter')

    if (is.numeric(rows_idx)){

      p <- p %>% plotly::add_trace(data = group_rows,
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
                                   type = 'scatter')
    }
   if (is.numeric(cols_idx)){
     p <- p %>% plotly::add_trace(data = group_cols,
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
                                  type = 'scatter')
    }

    p <- p %>% plotly::layout(title = paste('Association Plot \n', ncol(caobj$U), ' first dimensions, ', length(caobj$group),' samples.\n'),
             xaxis = list(title = 'Distance from origin (x)', rangemode = "tozero"),
             yaxis = list(title = 'Distance from gene to sample line (y)', rangemode = "tozero"), showlegend = TRUE)

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
#' @param mark_cols Character vector. Names of cols that should be highlighted in the plot.
#' @param reps Integer. Number of permutations during scoring. Default 3.
#' @param python A logical value indicating whether to use singular-value decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx shouls be labeled with text. Default TRUE.
#' @param type "ggplot"/"plotly". For a static plot a string "ggplot", for an interactive plot "plotly". Default "plotly".
#' @param ... Arguments forwarded to methods.
#' @export
runAPL <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, mark_cols = caobj$group, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, type = "plotly", ...){
  UseMethod("runAPL")
}



#' @rdname runAPL
#' @export
runAPL.default <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, mark_cols = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, type = "plotly", ...){
  stop(paste0("runAPL does not know how to handle objects of class ",
              class(x),
              ". Currently only objects of class 'matrix' or objects coercible to one, 'Seurat' or 'SingleCellExperiment' are supported."))
}



#' @export
#' @rdname runAPL
runAPL.matrix <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, mark_cols = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, type = "plotly", ...){

  if (!is(obj, "matrix")){
    obj <- as.matrix(obj)
  }
  stopifnot(is(obj, "matrix"))

  if (is.null(caobj)){
    caobj <- cacomp(obj = obj,
                    coords = TRUE,
                    top = top,
                    princ_coords = 3,
                    dims = dims,
                    python = python)

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
                         reps= reps,
                         python = python)

      mark_rows <- head(caobj$APL_score$Row_num, nrow)

      if (is(mark_cols, "character")){
        mark_cols <- match(mark_cols, rownames(caobj$apl_cols))
        mark_cols <- na.omit(mark_cols)
        if (anyNA(mark_cols)){
          warning("Not all names in 'mark_cols' are contained in the row names. Maybe they were filtered out?\n Non-matching values were ignored.")
        }
      }
    }
  } else{
    if (is(mark_rows, "character")){
      mark_rows <- match(mark_rows, rownames(caobj$apl_rows))
      mark_rows <- na.omit(mark_rows)
      if (anyNA(mark_rows)){
        warning("Not all names in 'mark_rows' are contained in the row names. Maybe they were filtered out?\n Non-matching values were ignored.")
      }
    }

    if (is(mark_cols, "character")){
      mark_cols <- match(mark_cols, rownames(caobj$apl_cols))
      mark_cols <- na.omit(mark_cols)
      if (anyNA(mark_cols)){
        warning("Not all names in 'mark_cols' are contained in the row names. Maybe they were filtered out?\n Non-matching values were ignored.")
      }
    }
  }

  # if(!is(mark_rows, "character")) stop("Parameter mark_rows hast to be of type 'character'.")
  # if(!is(mark_cols, "character")) stop("Parameter mark_cols hast to be of type 'character'.")

  p <- apl(caobj = caobj,
           type = type,
           rows_idx = mark_rows,
           cols_idx = mark_cols,
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
runAPL.SingleCellExperiment <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, mark_cols = NULL,  reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, type = "plotly", ..., assay = "counts"){

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
                mark_cols = mark_cols,
                nrow = nrow,
                top = top,
                score = score,
                reps = reps,
                python = python,
                row_labs = row_labs,
                col_labs = col_labs,
                type = type)

}


#' @description
#' runAPL.Seurat: Computes singular value decomposition and coordinates for the Association plot from Seurat objects, optionally with a DimReduc Object in the "CA" slot.
#'
#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#'
#' @rdname runAPL
#' @export
runAPL.Seurat <- function(obj, group, caobj = NULL, dims = NULL, nrow = 10, top = 5000, score = TRUE, mark_rows = NULL, mark_cols = NULL, reps = 3, python = TRUE, row_labs = TRUE, col_labs = TRUE, type = "plotly", ..., assay = DefaultAssay(obj)){

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
                mark_cols = mark_cols,
                nrow = nrow,
                row_labs = row_labs,
                col_labs = col_labs,
                type = type)
}

