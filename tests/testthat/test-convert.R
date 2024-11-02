# context("test conversion to and from cacomp")

load("./testdata/smoke.rda")
load("./testdata/smoke_scRNAseq.rda")
set.seed(2358)

d <- min(nrow(smoke), ncol(smoke)) - 1
ca <- cacomp(smoke, top = nrow(smoke), dims = d, princ_coords = 3)

test_that("check recompute function", {

  calist <- APL::as.list(ca)

  calist_sub <- calist[c("D",
                         "std_coords_cols",
                         "std_coords_rows",
                         "params")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_cols",
                         "std_coords_rows",
                         "prin_coords_rows",
                         "params")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("V",
                         "U",
                         "D",
                         "params")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_rows",
                         "V",
                         "params")]
  expect_equal(recompute(calist_sub, smoke), ca)

  calist_sub <- calist[c("std_coords_cols",
                         "std_coords_rows",
                         "prin_coords_rows",
                         "params")]
  expect_error(recompute(calist_sub, smoke[1:3, ]), "mat does not have have the correct number of rows.")
  expect_error(recompute(calist_sub, smoke[, 1:3]), "mat does not have have the correct number of columns.")

})

# d <- min(nrow(smoke), ncol(smoke)) - 1
# seu <- SeuratObject::CreateSeuratObject(smoke)
# seu <- cacomp(seu,
#               princ_coords = 3,
#               return_input = TRUE,
#               dims = d,
#               assay = "RNA",
#               slot = "counts")
#
# sce <- SingleCellExperiment::SingleCellExperiment(list(counts = smoke))
# sce <- cacomp(
#     sce,
#     dims = 3,
#     princ_coords = 3,
#     return_input = TRUE,
#     assay = "counts"
# )
# save(seu, sce, file = "./tests/testthat/testdata/smoke_scRNAseq.rda")

test_that("check Seurat integration", {
  expect_equal(as.cacomp(seu, assay = "RNA", slot = "counts"), ca)
})

test_that("check SingleCellExperiment integration", {
  expect_equal(as.cacomp(sce, assay = "counts"), ca)
})
