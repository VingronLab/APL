# context("test conversion to and from cacomp")

load("./testdata/smoke.rda")
load("./testdata/smoke_scRNAseq.rda")
ca <- cacomp(smoke, top = nrow(smoke), princ_coords = 3)

test_that("check recompute function", {

  calist <- APL::as.list(ca)

  calist_sub <- calist[c("D", "std_coords_cols", "std_coords_rows")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_cols", "std_coords_rows", "prin_coords_rows")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("V", "U", "D")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_rows", "V")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_cols", "std_coords_rows", "prin_coords_rows")]
  expect_error(recompute(calist_sub, smoke[1:3, ]), "mat does not have have the correct number of rows.")
  expect_error(recompute(calist_sub, smoke[, 1:3]), "mat does not have have the correct number of columns.")

})

# seu <- CreateSeuratObject(smoke)
# seu <- cacomp(seu, princ_coords = 3, return_input = TRUE, assay = "RNA", slot = "counts")
#
# sce <- SingleCellExperiment(list(counts=smoke))
# sce <- cacomp(sce, princ_coords = 3, return_input = TRUE, assay = "counts")
# save(seu, sce, file = "/home/kohl/PhD/gits/APL/tests/testthat/testdata/smoke_scRNAseq.rda")

test_that("check Seurat integration", {
  expect_equal(as.cacomp(seu, assay="RNA", slot = "counts"), ca)
})

test_that("check SingleCellExperiment integration", {
  expect_equal(as.cacomp(sce, assay = "counts"), ca)
})
