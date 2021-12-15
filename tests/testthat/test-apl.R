
# tab <- read.delim(file = "/home/kohl/PhD/gits/APL/tests/testthat/testdata/input_data.tsv")
# mat <- as.matrix(tab[,-1])
# rownames(mat) <- tab$Country.Name
# save(mat, file = "/home/kohl/PhD/gits/APL/tests/testthat/testdata/countries.rda")

load("./testdata/countries.rda")

grp <- c(6, 7, 8, 10, 12)


ca <- cacomp(mat, princ_coords = 3, dims = 19, top = 39)
ca <- apl_coords(ca, group = grp)



test_that("Example 1, 39 genes and 19 dimensions", {

  samples1 <- read.delim(file = "./testdata/AP_coordinates/example1/AP_coordinates_samples.txt")
  samples1 <- t(samples1)
  rownames(samples1) <- colnames(mat)
  colnames(samples1) <- c("x", "y")

  genes1 <- read.delim(file = "./testdata/AP_coordinates/example1/gene_ranking.txt")
  ord <- order(as.numeric(rownames(genes1)))
  rwnms <- rownames(mat)[as.numeric(rownames(genes1))[ord]]

  genes1_sort <- as.matrix(genes1[ord, c("x.coordinate","y.coordinate")])
  dimnames(genes1_sort) <- list(rwnms, c("x", "y"))

  ca <- cacomp(mat, princ_coords = 3, dims = 19, top = 39)
  ca <- apl_coords(ca, group = grp)

  expect_equal(ca@apl_cols, samples1, tolerance = 1e-8)
  expect_equal(ca@apl_rows, genes1_sort, tolerance = 1e-8)

 })

test_that("Example 2, 39 genes and 4 dimensions",{
  samples2 <- read.delim(file = "./testdata/AP_coordinates/example2/AP_coordinates_samples.txt")
  samples2 <- t(samples2)
  rownames(samples2) <- colnames(mat)
  colnames(samples2) <- c("x", "y")

  genes2 <- read.delim(file = "./testdata/AP_coordinates/example2/gene_ranking.txt")
  ord <- order(as.numeric(rownames(genes2)))
  rwnms <- rownames(mat)[as.numeric(rownames(genes2))[ord]]

  genes2_sort <- as.matrix(genes2[ord, c("x.coordinate","y.coordinate")])
  dimnames(genes2_sort) <- list(rwnms, c("x", "y"))

  ca <- cacomp(mat, princ_coords = 3, dims = 4, top = 39)
  ca <- apl_coords(ca, group = grp)

  expect_equal(ca@apl_cols, samples2, tolerance = 1e-8)
  expect_equal(ca@apl_rows, genes2_sort, tolerance = 1e-8)
})


test_that("Example 3, 20 genes and 4 dimensions",{
  samples3 <- read.delim(file = "./testdata/AP_coordinates/example3/AP_coordinates_samples.txt")
  samples3 <- t(samples3)
  rownames(samples3) <- colnames(mat)
  colnames(samples3) <- c("x", "y")

  genes3 <- read.delim(file = "./testdata/AP_coordinates/example3/gene_ranking.txt")
  ord <- order(as.numeric(rownames(genes3)))
  rwnms <- rownames(mat)[as.numeric(rownames(genes3))[ord]]

  genes3_sort <- as.matrix(genes3[ord, c("x.coordinate","y.coordinate")])
  dimnames(genes3_sort) <- list(rwnms, c("x", "y"))

  ca <- cacomp(mat, princ_coords = 3, dims = 4, top = 20)
  ca <- apl_coords(ca, group = grp)

  expect_equal(ca@apl_cols, samples3, tolerance = 1e-8)
  expect_equal(ca@apl_rows[order(rownames(ca@apl_rows)),], genes3_sort[order(rownames(genes3_sort)),], tolerance = 1e-8)
})

