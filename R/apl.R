#' @include constructor.R
NULL

#' Calculate Association Plot coordinates
#'
#' @description
#' Calculates the Association Plot coordinates for either the rows, 
#' columns or both (default).
#'
#' @details
#' Coordinates (x,y) of row vector \eqn{\vec{r}} are defined as
#' \deqn{x(\vec{r}) := \left|\vec{r}\right|\cos(\phi(\vec{r}))}
#' \deqn{y(\vec{r}) := \left|\vec{r}\right|\sin(\phi(\vec{r}))}
#' The x-direction is determined by calculating the centroid of the columns 
#' selected with the indices in "group".
#'
#' @return
#' Returns input "cacomp" object and adds components "apl_rows" and/or 
#' "apl_cols" for row and column coordinates.
#' In "group" the indices of the columns used to calculate the 
#' centroid are saved.
#' 
#' @references
#' Association Plots: Visualizing associations in high-dimensional 
#' correspondence analysis biplots
#' Elzbieta Gralinska, Martin Vingron
#' bioRxiv 2020.10.23.352096; doi: https://doi.org/10.1101/2020.10.23.352096
#'
#' @param caobj A "cacomp" object with principal row coordinates and 
#' standardized column coordinates calculated.
#' @param group Numeric/Character. Vector of indices or column names of 
#' the columns to calculate centroid/x-axis direction.
#' @param calc_rows TRUE/FALSE. Whether apl row coordinates should 
#' be calculated. Default TRUE.
#' @param calc_cols TRUE/FALSE. Whether apl column coordinates should 
#' be calculated. Default TRUE.
#' @export
#'
#' @examples
#' set.seed(1234)
#' # Simulate scRNAseq data
#' cnts <- data.frame(cell_1 = rpois(10, 5),
#'                    cell_2 = rpois(10, 10),
#'                    cell_3 = rpois(10, 20),
#'                    cell_4 = rpois(10, 20))
#' rownames(cnts) <- paste0("gene_", 1:10)
#' cnts <- as.matrix(cnts)
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#' # Calculate APL coordinates
#' ca <- apl_coords(ca, group = 3:4)
apl_coords <- function(caobj, group, calc_rows = TRUE, calc_cols = TRUE){

  stopifnot(is(caobj, "cacomp"))

  rows <- caobj@prin_coords_rows
  cols <- caobj@std_coords_cols
  cent <- cols


  if (is(group, "numeric")){
    subgroup <- cent[group,]
  } else if (is(group, "character")){
    idx <- match(group, rownames(cent))
    idx <- na.omit(idx)
    group <- idx
    subgroup <- cent[idx,]

    if (anyNA(idx)){
      warning("Not all names in 'group' are contained in the column names. ",
              "Non-matching values were ignored.")
    }
  } else {
    stop("Parameter group has to be either of type 'numeric' or 'character'.")
  }

  if (length(group) == 1){
    avg_group_coords <- subgroup # single sample
  } else {
    avg_group_coords <- colMeans(subgroup) # centroid vector.
  }
  length_vector_group <- sqrt(drop(avg_group_coords %*% avg_group_coords))
  length_vector_rows <- sqrt(rowSums(rows^2))
  length_vector_cols <- sqrt(rowSums(cols^2))

  if (calc_rows == TRUE){
    # message("Calculating APL row coordinates ...")
    # r⋅X = |r|*|X|*cosθ
    # x(r) = (r⋅X)/|X| = |r|*cosθ
    rowx <- drop(rows %*% avg_group_coords)/length_vector_group
    # pythagoras, y(r)=b²=c²-a²
    rowy <- sqrt(length_vector_rows^2 - rowx^2)

    rowx[is.na(rowx)] <- 0
    rowy[is.na(rowy)] <- 0
    # rowx[is.infinite(rowx)] <- 0
    # rowy[is.infinite(rowy)] <- 0

    caobj@apl_rows <- cbind("x"=rowx, "y"=rowy)
  }


  if (calc_cols == TRUE){
    # message("Calculating APL column coordinates ...")

    colx <- drop(cols %*% avg_group_coords)/length_vector_group
    coly <- sqrt(length_vector_cols^2 - colx^2)

    colx[is.na(colx)] <- 0
    coly[is.na(coly)] <- 0
    # colx[is.infinite(colx)] <- 0
    # coly[is.infinite(coly)] <- 0

    caobj@apl_cols <- cbind("x"=colx, "y"=coly)
  }

  if (is(group, "numeric")){
    caobj@group <- group
  } else if (is(group, "character")){
    idx <- match(group, rownames(cols))
    idx <- na.omit(idx)
    caobj@group <- idx
  }

  stopifnot(validObject(caobj))
  return(caobj)
}

#' Find rows most highly associated with a condition
#'
#' @description
#' Ranks rows by a calculated score which balances the association of the row 
#' with the condition and how associated it is with other conditions.
#'
#' @details
#' The score is calculated by permuting the values of each row to determine the 
#' cutoff angle of the 99% quantile.
#' \deqn{S_{alpha}(x,y)=x-\frac{y}{\tan\alpha}}
#' By default the permutation is repeated 10 times, but for very large matrices 
#' this can be reduced.
#' If store_perm is TRUE the permuted data is stored in the cacomp object and 
#' can be used for future scoring.
#' @return
#' Returns the input "cacomp" object with "APL_score" component added.
#' APL_score contains a data frame with ranked rows, their score and their 
#' original row number.
#' 
#' @references
#' Association Plots: Visualizing associations in high-dimensional 
#' correspondence analysis biplots \cr
#' Elzbieta Gralinska, Martin Vingron \cr
#' bioRxiv 2020.10.23.352096; doi: https://doi.org/10.1101/2020.10.23.352096
#'
#' @param caobj A "cacomp" object with principal row coordinates and 
#' standardized column coordinates calculated.
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression 
#' values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param dims Integer. Number of CA dimensions to retain. Needs to be the same 
#' as in caobj!
#' @param group Vector of indices of the columns to calculate centroid/x-axis 
#' direction.
#' @param reps Integer. Number of permutations to perform. Default = 10.
#' @param quant Numeric. Single number between 0 and 1 indicating the quantile 
#' used to calculate the cutoff. Default 0.99.
#' @param python A logical value indicating whether to use singular-value 
#' decomposition from the python package torch.
#' @param store_perm Logical. Whether permuted data should be stored in the CA 
#' object.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @export
#'
#' @examples
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:20, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis.
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#'
#' # Calculate APL coordinates:
#' ca <- apl_coords(ca, group = 1:10)
#'
#' # Rank genes by S-alpha score
#' ca <- apl_score(ca, mat = cnts)
apl_score <- function(caobj,
                      mat,
                      dims = caobj@dims,
                      group = caobj@group,
                      reps=10,
                      quant = 0.99,
                      python = FALSE,
                      store_perm = TRUE){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.empty(caobj@apl_rows)){
    stop("Please run apl_coords() first!")
  }

  names <- colnames(mat)
  row_num <- nrow(caobj@apl_rows)
  margin <- 1
  pc <- 1
  cc <- FALSE
  cr <- TRUE

  if (is(group, "character")){
    idx <- match(group, names)
    idx <- na.omit(idx)
    group <- idx
  }

  if (caobj@dims == 1 && !is.empty(caobj@dims)){
    row_num <- 1
  }

  apl_perm <- data.frame("x" = rep(0, row_num*reps),
                         "y" = rep(0, row_num*reps))
  saved_ca <- list()
  pb <- txtProgressBar(min = 0, max = reps, style = 3)

  for (k in seq(reps)){

    #permute rows and rerun cacomp

    if(isTRUE(store_perm) & identical(reps, attr(caobj@permuted_data,'reps'))){
      calist <- caobj@permuted_data[[k]][seq_len(3)]
      mat <- caobj@permuted_data[[k]]$mat
      mat <- mat[rownames(mat) %in% rownames(calist$std_coords_rows),]
      
      caobjp <- recompute(calist, mat)

    } else {
      mat_perm <- t(apply(mat, margin, FUN=sample))
      colnames(mat_perm) <- colnames(mat)


      suppressWarnings(caobjp <- cacomp(obj = mat_perm,
                                         python = python,
                                         coords = TRUE,
                                         princ_coords = pc,
                                         dims = dims,
                                         top = caobj@top_rows,
                                         inertia = FALSE))

      if(isTRUE(store_perm)){
        x <- list("std_coords_cols" = caobjp@std_coords_cols,
                  "std_coords_rows" = caobjp@std_coords_rows,
                  "D" = caobjp@D,
                  "mat" = mat_perm)

        saved_ca[[k]] <- x

      }
    }

    caobjp <- apl_coords(caobj = caobjp,
                         group = group,
                         calc_cols = cc,
                         calc_rows = cr)
    idx <- ((seq_len(row_num)+((k-1)*row_num)))

    apl_perm[idx,] <- caobjp@apl_rows

    setTxtProgressBar(pb, k)

  }

  close(pb)

  apl_perm[,3] <- apl_perm[,1]/apl_perm[,2] # cotan between row and x axis
  apl_perm[,3][is.na(apl_perm[,3])] <- 0

  cutoff_cotan <- quantile(apl_perm[,3], quant)

  score <- caobj@apl_rows[,1] - (caobj@apl_rows[,2] * cutoff_cotan)
  ranking <- data.frame("Rowname" = rownames(caobj@apl_rows),
                        "Score" = score,
                        "Row_num" = seq_len(nrow(caobj@apl_rows)))

  ranking <- ranking[order(ranking$Score, decreasing = TRUE),]
  ranking$Rank <- seq_len(nrow(ranking))

  caobj@APL_score <- ranking

  if(isTRUE(store_perm) & !identical(reps, attr(caobj@permuted_data,'reps'))){
    caobj@permuted_data <- saved_ca
    attr(caobj@permuted_data,'cutoff') <- cutoff_cotan
    attr(caobj@permuted_data,'reps') <- reps
  }

  stopifnot(validObject(caobj))
  return(caobj)

}




#' Run Gene overrepresentation analysis with topGO
#'
#' @description
#' This function uses the Kolmogorov-Smirnov test as implemented by the package 
#' topGO to test for overrepresentation in Gene Ontology gene sets.
#'
#' @details
#' For a chosen group of cells/samples,
#' the top 'ngenes' group specific genes are used for gene overrepresentation 
#' analysis.
#' The genes are ranked either by the precomputed APL score, or, if
#' not available by their APL x-coordinates.
#'
#' @return
#' A data.frame containing the gene sets with the highest overrepresentation.
#'
#' @references 
#' Adrian Alexa and Jorg Rahnenfuhrer \cr
#' topGO: Enrichment Analysis for Gene Ontology. \cr
#' R package version 2.42.0.
#' 
#' @param caobj A "cacomp" object with principal row coordinates and 
#' standardized column coordinates calculated.
#' @param ontology Character string. Chooses GO sets for 'BP' 
#' (biological processes), 'CC' (cell compartment) or 'MF' (molecular function).
#' @param organism Character string. Either 'hs' (homo sapiens), 'mm' 
#' (mus musculus) or the name of the organism package such as 'org.*.eg.db'.
#' @param ngenes Numeric. Number of top ranked genes to test for 
#' overrepresentation.
#' @param score_cutoff numeric. S-alpha score cutoff. Only genes with a score 
#' larger will be tested.
#' @param use_coords Logical. Whether the x-coordinates of the row APL 
#' coordinates should be used for ranking.
#' Only recommended when no S-alpha score (see apl_score()) can be calculated.
#' @param return_plot Logical. Whether a plot of significant gene sets should 
#' be additionally returned.
#' @param top_res Numeric. Number of top scoring genes to plot.
#' 
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
apl_topGO <- function(caobj,
                      ontology,
                      organism="hs",
                      ngenes = 1000,
                      score_cutoff = 0,
                      use_coords = FALSE,
                      return_plot = FALSE,
                      top_res = 15){

  topGO::groupGOTerms()

  if(!is.empty(caobj@APL_score) & !isTRUE(use_coords)){

    Score_ord <- caobj@APL_score[order(caobj@APL_score$Score, decreasing=TRUE),]
    sel <- which(!Score_ord$Score >= score_cutoff)[1]
    sel <- sel-1
    ranked_genes <- seq_len(nrow(Score_ord))
    names(ranked_genes) <- Score_ord$Rowname

  } else if (isTRUE(use_coords) & is.empty(caobj@APL_score)) {

    if(is.empty(caobj@apl_rows)){
      stop("No APL coordinates found for rows. Please first run apl_coords.\n")
    }
    if (ngenes > nrow(caobj@apl_rows)){
      stop("ngenes is larger than the total number of genes.\n")
    } else if (ngenes == nrow(caobj@apl_rows)){
      warning("You have selected all available genes.\n")
    }
    ranked_genes <- seq_len(nrow(caobj@apl_rows))
    names(ranked_genes) <- rownames(caobj@apl_rows)[order(caobj@apl_rows[,1],
                                                          decreasing = TRUE)]
    sel <- ngenes
  } else {
    stop("APL scores not present but use_coords set to FALSE.\n")
  }

  if (organism == "hs"){
    organism <- "org.Hs.eg.db"
  } else if (organism == "mm"){
    organism <- "org.Mm.eg.db"
  } else {
    warning("Custom organism chosen.\n")
    organism <- organism
  }

  if (!ontology %in% c("BP", "CC", "MF")){
    stop("Please choose one of the following ontologies: 'BP', 'CC' or 'MF'.\n")
  }


  gene_sets <- topGO::annFUN.org(whichOnto=ontology,
                                 feasibleGenes=NULL,
                                 mapping=organism,
                                 ID="symbol")

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = ranked_genes,
                annot = topGO::annFUN.GO2genes,
                GO2genes = gene_sets,
                geneSelectionFun = function(x) x<sel,
                nodeSize=5)

  results_test <- topGO::runTest(GOdata,
                                 algorithm = "elim",
                                 statistic = "fisher")

  goEnrichment <- topGO::GenTable(GOdata,
                                  raw.p.value = results_test,
                                  topNodes = length(results_test@score),
                                  numChar = 1200)
  if (isTRUE(return_plot)){

    if(top_res > length(results_test@score)){
      warning("More top nodes selected via top_res than available. ",
              "Returning max. number of nodes instead.\n")
      top_res <- min(top_res, length(results_test@score))
    }
    topGO::showSigOfNodes(GOdata,
                          score(results_test),
                          firstSigNodes = top_res,
                          useInfo = 'def')
  }

  return(goEnrichment)

}

#' 
#' Plot Association Plot with ggplot
#' @description 
#' Uses ggplot to plot an Association Plot
#' 
#' @param rows Row APL-coordinates
#' @param rows_group Row AP-coordinates to highlight
#' @param rows_scored Row AP-coordinates of rows above a score cutoff.
#' @param cols Column AP-coordinates
#' @param cols_group Column AP-coordinates for the group to be highlighted.
#' @param rows_color Color for rows
#' @param rows_high_color Color for rows to be highlighted.
#' @param cols_color Column points color.
#' @param cols_high_color Color for column points to be highlighted..
#' @param score_color Color scheme for row points with a score.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx 
#' should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx 
#' shouls be labeled with text. Default FALSE.
#' @param show_score Logical. Whether the S-alpha score should be shown in 
#' the plot.
#' @param show_cols Logical. Whether column points should be plotted.
#' @param show_rows Logical. Whether row points should be plotted.
#' 
#' @returns
#' ggplot Association Plot
apl_ggplot <- function(rows,
                       rows_group = NULL,
                       cols,
                       cols_group = NULL,
                       rows_scored = NULL,
                       rows_color = "#0066FF",
                       rows_high_color = "#FF0000",
                       cols_color = "#601A4A",
                       cols_high_color = "#EE442F",
                       score_color = "rainbow",
                       row_labs = FALSE,
                       col_labs = FALSE,
                       show_score = FALSE,
                       show_cols = FALSE,
                       show_rows = TRUE){
  p <- ggplot2::ggplot()
  
  
  if (isTRUE(show_rows)){
    p <- p +
      ggplot2::geom_point(data=rows,
                          ggplot2::aes(x=.data$x, y=.data$y),
                          color = rows_color,
                          alpha = 0.7,
                          shape = 1) #16 point, 1 circle.
    
    if (isTRUE(show_score)){
      p <- p +
        ggplot2::geom_point(data=rows_scored,
                            ggplot2::aes(x=.data$x,
                                         y=.data$y,
                                         fill = .data$Score),
                            alpha = 0.7,
                            shape = 21,
                            stroke = 0.5)
      if (score_color == "rainbow"){
        hex <- c("#00007F",
                 "blue",
                 "#007FFF",
                 "cyan",
                 "#7FFF7F",
                 "yellow",
                 "#FF7F00",
                 "red",
                 "#7F0000")
        p <- p +
          ggplot2::scale_fill_gradientn(colours = hex,
                                        name = expression(S*alpha))
      } else if (score_color == "viridis"){
        p <- p +
          ggplot2::scale_color_viridis_c(option = "D")
      }
      
    }
    
    if (!is.null(rows_group)){
      p <- p +
        ggplot2::geom_point(data=rows_group,
                            ggplot2::aes(x=.data$x, y=.data$y),
                            color = rows_high_color,
                            shape = 19)
      
      if(isTRUE(row_labs)){
        p <- p +
          ggrepel::geom_text_repel(data = rows_group,
                                   ggplot2::aes(x=.data$x,
                                                y=.data$y,
                                                label=.data$rownms),
                                   color = rows_high_color,
                                   max.overlaps = Inf)
      }
    }
  }
  
  if (isTRUE(show_cols)){
    p <- p +
      ggplot2::geom_point(data=cols,
                          ggplot2::aes(x=.data$x, y=.data$y),
                          color = cols_color,
                          shape = 4)
    
      if (!is.null(cols_group)){
        p <- p +
        ggplot2::geom_point(data=cols_group,
                            ggplot2::aes(x=.data$x, y=.data$y),
                            color = cols_high_color,
                            shape = 4)
        
        if(col_labs == TRUE){
          p <- p +
            ggrepel::geom_text_repel(data=cols_group,
                                     ggplot2::aes(x=.data$x,
                                                  y=.data$y,
                                                  label=.data$rownms),
                                     color = cols_high_color)
        }
      }

  }
  
  p <- p +
    ggplot2::labs(title="Association Plot") +
    ggplot2::theme_bw()
  
  return(p)
}


#' 
#' Plot Association Plot with plotly
#' @description 
#' Uses plotly to generate an interactive Association Plot
#' 
#' @param rows Row APL-coordinates
#' @param rows_group Row AP-coordinates to highlight
#' @param rows_scored Row AP-coordinates of rows above a score cutoff.
#' @param cols Column AP-coordinates
#' @param cols_group Column AP-coordinates for the group to be highlighted.
#' @param rows_color Color for rows
#' @param rows_high_color Color for rows to be highlighted.
#' @param cols_color Column points color.
#' @param cols_high_color Color for column points to be highlighted.
#' @param score_color Color scheme for row points with a score.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx 
#' should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx 
#' shouls be labeled with text. Default FALSE.
#' @param show_score Logical. Whether the S-alpha score should be shown in 
#' the plot.
#' @param show_cols Logical. Whether column points should be plotted.
#' @param show_rows Logical. Whether row points should be plotted.
#' 
#' @returns
#' Interactive plotly Association Plot
apl_plotly <- function(rows,
                       rows_group =NULL,
                       cols,
                       cols_group,
                       rows_scored = NULL,
                       rows_color = "#0066FF",
                       rows_high_color = "#FF0000",
                       cols_color = "#601A4A",
                       cols_high_color = "#EE442F",
                       score_color = "rainbow",
                       row_labs = FALSE,
                       col_labs = FALSE,
                       show_score = FALSE,
                       show_cols = FALSE,
                       show_rows = TRUE){
  
  if (row_labs == FALSE){
    rlabs <- 'markers'
    rowfont <- NULL
  } else {
    rlabs <- 'markers+text'
    rowfont <- list(color=rows_high_color)
    
  }
  
  if (col_labs == FALSE){
    clabs <- 'markers'
    colfont <- NULL
  } else {
    clabs <- 'markers+text'
    colfont <- list(color = cols_high_color)
  }
  
  p <- plotly::plot_ly()
  
  if(isTRUE(show_rows)){
    
    color_fix <- rows_color
    colors_fun <- NULL
    color_bar <- NULL
    sym <- "circle-open"
    
    p <- p %>%
      plotly::add_trace(data = rows,
                        x = ~x,
                        y = ~y,
                        mode = 'markers',
                        text = rows$rownms,
                        opacity = 0.7,
                        textposition = "left",
                        marker = list(color = color_fix, 
                                      colorscale = colors_fun,
                                      symbol = sym,
                                      colorbar = color_bar,
                                      size = 5),
                        name = 'genes',
                        hoverinfo = 'text',
                        type = 'scatter')
    
    if (isTRUE(show_score)){
      
      if (score_color == "rainbow"){
        colors_fun <- "Jet"
      } else if (score_color == "viridis"){
        colors_fun <- 'Viridis'
      }
      color_fix <- as.formula("~Score")
      color_bar <- list(title = "S\u03B1", len=0.5)
      sym <- "circle"
      
      p <- p %>%
        plotly::add_trace(data = rows_scored,
                          x = ~x,
                          y = ~y,
                          mode = 'markers',
                          text = rows_scored$rownms,
                          opacity = 0.7,
                          textposition = "left",
                          marker = list(color = color_fix,
                                        colorscale = colors_fun,
                                        symbol = sym,
                                        colorbar = color_bar,
                                        size = 5,
                                        line = list(color = '#000000', #black
                                                    width = 0.5)),
                          name = 'genes (scored)',
                          hoverinfo = 'text',
                          type = 'scatter')
      
    }
    
    if (!is.null(rows_group)){
      
      p <- p %>%
        plotly::add_trace(data = rows_group,
                          x = ~x,
                          y = ~y,
                          mode = rlabs,
                          text = rows_group$rownms,
                          textposition = "left",
                          textfont=rowfont,
                          marker = list(symbol = 'circle',
                                        color = rows_high_color,
                                        size = 5),
                          name = 'marked genes',
                          hoverinfo = 'text',
                          type = 'scatter')
    }
  }
  
  if(isTRUE(show_cols)){
    p <- p %>%
      plotly::add_trace(data=cols,
                        x = ~x,
                        y =  ~y,
                        mode = 'markers',
                        text = cols$rownms,
                        textposition = "left",
                        marker = list(color = cols_color,
                                      symbol = 'x',
                                      size = 5),
                        name = 'samples',
                        hoverinfo = 'text',
                        type = 'scatter')
    
    if (!is.null(cols_group)){
      p <- p %>% plotly::add_trace(data = cols_group,
                                   x = ~x,
                                   y = ~y,
                                   mode = clabs,
                                   text = cols_group$rownms,
                                   textposition = "left",
                                   textfont=colfont,
                                   marker = list(symbol = 'x',
                                                 color = cols_high_color,
                                                 size = 5),
                                   name = 'marked samples',
                                   hoverinfo = 'text',
                                   type = 'scatter')
    }
  }
  
  
  p <- p %>% 
    plotly::layout(title = paste('Association Plot'),
                   xaxis = list(title = 'Distance from origin (x)',
                                rangemode = "tozero"),
                   yaxis = list(title = 'Distance from gene to sample line (y)',
                                rangemode = "tozero"),
                   showlegend = TRUE)
  
  return(p)
}

#' Association Plot
#'
#' @description
#' Plot an Association Plot for the chosen columns.
#'
#' @details
#' For an interactive plot type="plotly" can be chosen, otherwise a static plot 
#' will be returned.
#' The row and column coordinates have to be already calculated by 
#' `apl_coords()`.
#'
#' @return
#' Either a ggplot or plotly object.
#' 
#' @references
#' Association Plots: Visualizing associations in high-dimensional 
#' correspondence analysis biplots \cr
#' Elzbieta Gralinska, Martin Vingron \cr
#' bioRxiv 2020.10.23.352096; doi: https://doi.org/10.1101/2020.10.23.352096 \cr
#'
#' @param caobj  An object of class "cacomp" and "APL" with apl 
#' coordinates calculated.
#' @param type "ggplot"/"plotly". For a static plot a string "ggplot", 
#' for an interactive plot "plotly". Default "ggplot".
#' @param rows_idx numeric/character vector. 
#' Indices or names of the rows that should be labelled. Default NULL.
#' @param cols_idx numeric/character vector. 
#' Indices or names of the columns that should be labelled. 
#' Default is only to label columns making up the centroid: caobj@group.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx 
#' should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx 
#' shouls be labeled with text. Default FALSE.
#' @param show_score Logical. Whether the S-alpha score should be shown in 
#' the plot.
#' @param show_cols Logical. Whether column points should be plotted.
#' @param show_rows Logical. Whether row points should be plotted.
#' @param score_cutoff Numeric. Rows (genes) with a score >= score_cutoff will 
#' be colored according to their score if show_score = TRUE.
#' @param score_color Either "rainbow" or "viridis".
#' @export
#' @examples
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # Run correspondence analysis
#' ca <- cacomp(obj = cnts, princ_coords = 3)
#'
#' # Calculate APL coordinates for arbitrary group
#' ca <- apl_coords(ca, group = 1:10)
#'
#' # plot results
#' # Note:
#' # Due to random gene expression & group, no highly
#' # associated genes are visible.
#' apl(ca, type = "ggplot")
apl <- function(caobj,
                type="ggplot",
                rows_idx = NULL,
                cols_idx = caobj@group,
                row_labs = FALSE,
                col_labs = FALSE,
                show_score = FALSE,
                show_cols = FALSE,
                show_rows = TRUE,
                score_cutoff = 0,
                score_color = "rainbow"){

  if (!is(caobj,"cacomp")){
    stop("Not a CA object. Please run cacomp() and apl_coords() first!")
  }

  if (is.empty(caobj@apl_rows) || is.empty(caobj@apl_cols)){
    stop("Please run apl_coords() first!")
  }

  if (is(rows_idx, "character")){
    names <- rownames(caobj@apl_rows)
    idx <- match(rows_idx, names)
    idx <- na.omit(idx)
    rows_idx <- idx
  }

  if (is(cols_idx, "character")){
    names <- rownames(caobj@apl_cols)
    idx <- match(cols_idx, names)
    idx <- na.omit(idx)
    cols_idx <- idx
  }

  if (is.numeric(cols_idx)){
  group_cols <- data.frame(rownms = rownames(caobj@apl_cols)[cols_idx],
                           x = caobj@apl_cols[cols_idx,"x"],
                           y = caobj@apl_cols[cols_idx,"y"],
                           row.names = NULL)
  } else {
    group_cols <- NULL
  }
  if(is.numeric(rows_idx)){
  group_rows <- data.frame(rownms = rownames(caobj@apl_rows)[rows_idx],
                           x = caobj@apl_rows[rows_idx,"x"],
                           y = caobj@apl_rows[rows_idx,"y"],
                           row.names = NULL)
  } else {
    group_rows <- NULL
  }
  
  rows_color <- "#0066FF"
  rows_high_color <- "#FF0000"
  cols_color <- "#601A4A"
  cols_high_color <- "#EE442F"



  apl_rows.tmp <- data.frame(rownms = rownames(caobj@apl_rows),
                             caobj@apl_rows)
  apl_cols.tmp <- data.frame(rownms = rownames(caobj@apl_cols),
                             caobj@apl_cols)

  if (isTRUE(show_score)) {
    apl_scores <- caobj@APL_score$Score[order(caobj@APL_score$Row_num)]
    apl_rows.tmp$Score <- apl_scores
    idx <- which(apl_scores >= score_cutoff)

    if (is.empty(idx)) {
      show_score <- FALSE
      warning("No rows with a Score >= score_cutoff of ",
              score_cutoff, 
              ". Choose lower cutoff.")
    } else {
      apl_scored.tmp <- apl_rows.tmp[idx,]
      apl_rows.tmp <- apl_rows.tmp[-idx,]
    }

  }


  if (type == "ggplot"){
    
    p <- apl_ggplot(rows = apl_rows.tmp,
                    rows_group = group_rows,
                    cols = apl_cols.tmp,
                    cols_group = group_cols,
                    rows_scored = apl_scored.tmp,
                    rows_color = rows_color,
                    rows_high_color = rows_high_color,
                    cols_color = cols_color,
                    cols_high_color = cols_high_color,
                    score_color = score_color,
                    row_labs = row_labs,
                    col_labs = col_labs,
                    show_score = show_score,
                    show_cols = show_cols,
                    show_rows = show_rows)
    return(p)

  } else if (type == "plotly"){
    p <- apl_plotly(rows = apl_rows.tmp,
                    rows_group = group_rows,
                    cols = apl_cols.tmp,
                    cols_group = group_cols,
                    rows_scored = apl_scored.tmp,
                    rows_color = rows_color,
                    rows_high_color = rows_high_color,
                    cols_color = cols_color,
                    cols_high_color = cols_high_color,
                    score_color = score_color,
                    row_labs = row_labs,
                    col_labs = col_labs,
                    show_score = show_score,
                    show_cols = show_cols,
                    show_rows = show_rows)
    return(p)
    
  } else {
    stop("Please specify plot = \"ggplot\" or \"plotly\".",
         " Other options are not accepted.")
  }

}

#' Compute and plot Association Plot
#'
#' @description
#' Computes singular value decomposition and coordinates for 
#' the Association Plot.
#'
#' @details
#' The function is a wrapper that calls `cacomp()`, `apl_coords()`, 
#' `apl_score()` and finally `apl()` for ease of use.
#' The chosen defaults are most useful for genomics experiments, but for more 
#' fine grained control the functions
#' can be also run individually for the same results.
#' If score = FALSE, nrow and reps are ignored. If mark_rows is not NULL score 
#' is treated as if FALSE.
#' @return
#' Association Plot (plotly object).
#' 
#' @references
#' Association Plots: Visualizing associations in high-dimensional 
#' correspondence analysis biplots \cr
#' Elzbieta Gralinska, Martin Vingron \cr
#' bioRxiv 2020.10.23.352096; doi: https://doi.org/10.1101/2020.10.23.352096 \cr
#'
#' @param obj A numeric matrix, Seurat or SingleCellExperiment object. For 
#' sequencing a count matrix, gene expression values with genes in rows and 
#' samples/cells in columns.
#' Should contain row and column names.
#' @param caobj A "cacomp" object as outputted from `cacomp()`. If not supplied 
#' will be calculated. Default NULL.
#' @param dims Integer. Number of dimensions to keep. Default NULL (keeps all 
#' dimensions).
#' @param group Numeric/Character. Vector of indices or column names of the 
#' columns to calculate centroid/x-axis direction.
#' @param nrow Integer. The top nrow scored row labels will be added to the 
#' plot if score = TRUE. Default 10.
#' @param top Integer. Number of most variable rows to retain. Default 5000 
#' rows (set NULL to keep all).
#' @param score Logical. Whether rows should be scored and ranked. Ignored when 
#' a vector is supplied to mark_rows. Default TRUE.
#' @param mark_rows Character vector. Names of rows that should be highlighted 
#' in the plot. If not NULL, score is ignored. Default NULL.
#' @param mark_cols Character vector. Names of cols that should be highlighted 
#' in the plot.
#' @param reps Integer. Number of permutations during scoring. Default 3.
#' @param python A logical value indicating whether to use singular value 
#' decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` 
#' in R.
#' @param row_labs Logical. Whether labels for rows indicated by rows_idx 
#' should be labeled with text. Default TRUE.
#' @param col_labs Logical. Whether labels for columns indicated by cols_idx 
#' should be labeled with text. Default TRUE.
#' @param type "ggplot"/"plotly". For a static plot a string "ggplot", for an 
#' interactive plot "plotly". Default "plotly".
#' @param show_cols Logical. Whether column points should be plotted.
#' @param show_rows Logical. Whether row points should be plotted.
#' @param score_cutoff Numeric. Rows (genes) with a score >= score_cutoff
#' will be colored according to their score if show_score = TRUE.
#' @param score_color Either "rainbow" or "viridis".
#' @param pd_method Which method to use for pick_dims (\link[APL]{pick_dims}).
#' @param pd_reps Number of repetitions performed when using "elbow_rule" in 
#' `pick_dims`.
#' (\link[APL]{pick_dims})
#' @param pd_use Whether to use `pick_dims` (\link[APL]{pick_dims}) to determine
#' the number of dimensions. Ignored when `dims` is set by the user.
#' @param ... Arguments forwarded to methods.
#' @export
setGeneric("runAPL", function(obj,
                              group,
                              caobj = NULL,
                              dims = NULL,
                              nrow = 10,
                              top = 5000,
                              score = TRUE,
                              mark_rows = NULL,
                              mark_cols = caobj@group,
                              reps = 3,
                              python = FALSE,
                              row_labs = TRUE,
                              col_labs = TRUE,
                              type = "plotly",
                              show_cols = FALSE,
                              show_rows = TRUE,
                              score_cutoff = 0,
                              score_color = "rainbow",
                              pd_method = "elbow_rule",
                              pd_reps = 1,
                              pd_use = TRUE,
                              ...) {
  standardGeneric("runAPL")
})

#' @export
#' @rdname runAPL
#' @examples
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' # (nonsensical) APL
#' runAPL(obj = cnts,
#'        group = 1:10,
#'        dims = 10,
#'        top = 500,
#'        score = TRUE,
#'        show_cols = TRUE,
#'        type = "ggplot")
setMethod(f = "runAPL",
          signature=(obj="matrix"),
          function(obj,
                   group,
                   caobj = NULL,
                   dims = NULL,
                   nrow = 10,
                   top = 5000,
                   score = TRUE,
                   mark_rows = NULL,
                   mark_cols = NULL,
                   reps = 3,
                   python = FALSE,
                   row_labs = TRUE,
                   col_labs = TRUE,
                   type = "plotly",
                   show_cols = FALSE,
                   show_rows = TRUE,
                   score_cutoff = 0,
                   score_color = "rainbow",
                   pd_method = "elbow_rule",
                   pd_reps = 1,
                   pd_use = TRUE,
                   ...){

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
    
    if(is.null(dims) & isTRUE(pd_use)){
      
      stopifnot("Cannot use scree plot to pick dimensions" = pd_method != "scree_plot")
      pd <- pick_dims(obj = caobj,
                      mat = obj,
                      method = pd_method,
                      reps = pd_reps,
                      python = python,
                      return_plot = FALSE)
      
      message("\n Using ",
              pd,
              " dimensions. Subsetting.")
      
      dims <- pd
      caobj <- subset_dims(caobj, dims)
      
    } else if (is.null(dims) & !isTRUE(pd_use)) {
      dims <- caobj@dims
      message("\n Using ",
              dims,
              " dimensions. No subsetting performed.")
      
    } else if (!is.null(dims) & !isTRUE(pd_use)) {
      message("\n Using ", dims, " dimensions.")
    } else if(!is.null(dims) & isTRUE(pd_use)){
      message("\nBoth dims and pd_use set. Using dimensions specified by dims.")
      message("\nUsing ", dims, " dimensions.")
    }
    

  } else {
    if(!is.empty(caobj@dims) && is.null(dims)){
      if(isTRUE(pd_use)){
        stopifnot("Cannot use scree plot to pick dimensions" = 
                    pd_method != "scree_plot")
        pd <- pick_dims(obj = caobj,
                        mat = obj,
                        method = pd_method,
                        reps = pd_reps,
                        python = python,
                        return_plot = FALSE)
        
        message("\n Using ",
                pd,
                " dimensions")
        
        dims <- pd
        caobj <- subset_dims(caobj, dims)
      } else {
        dims <- caobj@dims
        message("\n Using ",
                dims,
                " dimensions. No subsetting performed.")
      }

      
    } else if (!is.empty(caobj@dims) && !is.null(dims)) {
      if (dims < caobj@dims){

        if (!isTRUE(pd_use)) {
          message("\n Using ", dims, " dimensions.")
        } else if(isTRUE(pd_use)){
          message("\nBoth dims and pd_use set. Using dimensions specified by dims.")
          message("\nUsing ",
                  dims,
                  " dimensions.")
        }
        
        caobj <- subset_dims(caobj, dims = dims)
      } else if(dims > length(caobj@D)){
        warning("dims is larger than the number of available dimensions. ",
                "Argument ignored")
      }
    }

    if (is.empty(caobj@prin_coords_rows) && !is.empty(caobj@std_coords_rows)){
      caobj <- ca_coords(caobj = caobj,
                         dims = dims,
                         princ_only = TRUE,
                         princ_coords = 1)
    } else if (is.empty(caobj@prin_coords_rows) || 
               is.empty(caobj@std_coords_cols)){
      caobj <- ca_coords(caobj = caobj,
                         dims = dims,
                         princ_only = FALSE,
                         princ_coords = 1)
    }
  }

  if (is.empty(caobj@apl_rows) || is.empty(caobj@apl_cols)){
    caobj <- apl_coords(caobj = caobj, group = group)
  }


  if (is.null(mark_rows)){
    if (score == TRUE){
      caobj <- apl_score(caobj = caobj,
                         mat = obj,
                         dims = caobj@dims,
                         group = caobj@group,
                         reps= reps,
                         python = python)

      mark_rows <- head(caobj@APL_score$Row_num, nrow)

      if (is(mark_cols, "character")){
        mark_cols <- match(mark_cols, rownames(caobj@apl_cols))
        mark_cols <- na.omit(mark_cols)
        if (anyNA(mark_cols)){
          warning(
            "Not all names in 'mark_cols' are contained in the row names.",
            " Maybe they were filtered out?",
            "\nNon-matching values were ignored.")
        }
      }
    }
  } else{
    if (is(mark_rows, "character")){
      mark_rows <- match(mark_rows, rownames(caobj@apl_rows))
      mark_rows <- na.omit(mark_rows)
      if (anyNA(mark_rows)){
        warning(
          "Not all names in 'mark_rows' are contained in the row names.",
          " Maybe they were filtered out?",
          "\n Non-matching values were ignored.")
      }
    }

    if (is(mark_cols, "character")){
      mark_cols <- match(mark_cols, rownames(caobj@apl_cols))
      mark_cols <- na.omit(mark_cols)
      if (anyNA(mark_cols)){
        warning(
          "Not all names in 'mark_cols' are contained in the row names.",
          " Maybe they were filtered out?",
          "\n Non-matching values were ignored.")
      }
    }
  }

  if(isTRUE(col_labs) & is.empty(mark_cols)){
    mark_cols <- caobj@group
  }

  p <- apl(caobj = caobj,
           type = type,
           rows_idx = mark_rows,
           cols_idx = mark_cols,
           row_labs = row_labs,
           col_labs = col_labs,
           show_score = score,
           show_cols = show_cols,
           show_rows = show_rows,
           score_cutoff = score_cutoff,
           score_color = score_color)

  return(p)
})



#' @description
#' runAPL.SingleCellExperiment: Computes singular value decomposition and 
#' coordinates for the Association Plot from SingleCellExperiment objects with 
#' reducedDim(obj, "CA") slot (optional).
#'
#' @param assay Character. The assay from which extract the count matrix for 
#' SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#'
#' @rdname runAPL
#' @export
#' @examples
#'
#' ########################
#' # SingleCellExperiment #
#' ########################
#' library(SingleCellExperiment)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' sce <- SingleCellExperiment(assays=list(counts=cnts))
#'
#' # (nonsensical) APL
#' runAPL(obj = sce,
#'        group = 1:10,
#'        dims = 10,
#'        top = 500,
#'        score = TRUE,
#'        show_cols = TRUE,
#'        type = "ggplot",
#'        assay = "counts")
setMethod(f = "runAPL",
          signature=(obj="SingleCellExperiment"),
          function(obj,
                   group,
                   caobj = NULL,
                   dims = NULL,
                   nrow = 10,
                   top = 5000,
                   score = TRUE,
                   mark_rows = NULL,
                   mark_cols = NULL,
                   reps = 3,
                   python = FALSE,
                   row_labs = TRUE,
                   col_labs = TRUE,
                   type = "plotly",
                   show_cols = FALSE,
                   show_rows = TRUE,
                   score_cutoff = 0,
                   score_color = "rainbow",
                   pd_method = "elbow_rule",
                   pd_reps = 1,
                   pd_use = TRUE,
                   ...,
                   assay = "counts"){

  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" =
              is(obj, "SingleCellExperiment"))

  mat <- SummarizedExperiment::assay(obj, assay)

  if ("CA" %in% SingleCellExperiment::reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  }

  runAPL(obj = mat,
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
        type = type,
        show_cols = show_cols,
        show_rows = show_rows,
        score_cutoff = score_cutoff,
        score_color = score_color,
        pd_method = pd_method,
        pd_reps = pd_reps,
        pd_use = pd_use)

})


#' @description
#' runAPL.Seurat: Computes singular value decomposition and coordinates for 
#' the Association Plot from Seurat objects, optionally with a DimReduc Object 
#' in the "CA" slot.
#'
#' @param assay Character. The assay from which extract the count matrix for 
#' SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for 
#' SingleCellExperiments.
#' @param slot character. The Seurat assay slot from which to extract the 
#' count matrix.
#' @rdname runAPL
#' @export
#' @examples
#'
#' ###########
#' # Seurat  #
#' ###########
#' library(Seurat)
#' set.seed(1234)
#'
#' # Simulate counts
#' cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
#'                x = sample(1:100, 50, replace = TRUE))
#' rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
#' colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
#'
#' seu <- CreateSeuratObject(counts = cnts)
#'
#' # (nonsensical) APL
#' runAPL(obj = seu,
#'        group = 1:10,
#'        dims = 10,
#'        top = 500,
#'        score = TRUE,
#'        show_cols = TRUE,
#'        type = "ggplot",
#'        assay = "RNA",
#'        slot = "counts")
setMethod(f = "runAPL",
          signature=(obj="Seurat"),
          function(obj,
                   group,
                   caobj = NULL,
                   dims = NULL,
                   nrow = 10,
                   top = 5000,
                   score = TRUE,
                   mark_rows = NULL,
                   mark_cols = NULL,
                   reps = 3,
                   python = FALSE,
                   row_labs = TRUE,
                   col_labs = TRUE,
                   type = "plotly",
                   show_cols = FALSE,
                   show_rows = TRUE,
                   score_cutoff = 0,
                   score_color = "rainbow",
                   pd_method = "elbow_rule",
                   pd_reps = 1,
                   pd_use = TRUE,
                   ...,
                   assay = Seurat::DefaultAssay(obj),
                   slot = "counts"){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = slot)
  seu <- as.matrix(seu)

  if ("CA" %in% Seurat::Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay)
  } else {
    caobj <- NULL
  }

  runAPL(obj = seu,
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
        type = type,
        show_cols = show_cols,
        show_rows = show_rows,
        score_cutoff = score_cutoff,
        score_color = score_color,
        pd_method = pd_method,
        pd_reps = pd_reps,
        pd_use = pd_use)
})
