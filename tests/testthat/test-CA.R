# context("Correspondence Analysis")


# library(ca)
# data(smoke)
#
# smoke_ca <- ca(smoke)
#
# smoke_prin <- cacoord(smoke_ca,
#           type = c("principal"),
#           dim = NA,
#           rows = TRUE,
#           cols = TRUE)
#
# smoke <- as.matrix(smoke)
# save(smoke, smoke_ca, smoke_prin, file = "./tests/testthat/testdata/smoke.rda")

load("./testdata/smoke.rda")

ca_python <- cacomp(obj = smoke, top = nrow(smoke), princ_coords = 3, coords = TRUE, python = TRUE)
ca_svd <- cacomp(obj = smoke, top = nrow(smoke), princ_coords = 3, coords = TRUE, python = FALSE)
cac <- ca_coords(ca_svd, princ_coords = 3)


test_that("CA with torch svd results", {
    expect_equal(ca_python@dims, length(smoke_ca$sv))

    expect_equal(as.numeric(ca_python@D), smoke_ca$sv)
    expect_equal(ca_python@std_coords_cols, smoke_ca$colcoord)
    expect_equal(ca_python@std_coords_rows, smoke_ca$rowcoord)

    expect_equal(ca_python@prin_coords_cols, smoke_prin$columns)
    expect_equal(ca_python@prin_coords_rows, smoke_prin$rows)

    expect_equal(as.numeric(ca_python@row_masses), smoke_ca$rowmass)
    expect_equal(as.numeric(ca_python@row_inertia), smoke_ca$rowinertia)

    expect_equal(as.numeric(ca_python@col_masses), smoke_ca$colmass)
    expect_equal(as.numeric(ca_python@col_inertia), smoke_ca$colinertia)
})


test_that("CA with R svd results", {
    expect_equal(ca_svd@dims, length(smoke_ca$sv))

    expect_equal(as.numeric(ca_svd@D), smoke_ca$sv)
    expect_equal(ca_svd@std_coords_cols, smoke_ca$colcoord)
    expect_equal(ca_svd@std_coords_rows, smoke_ca$rowcoord)

    expect_equal(ca_svd@prin_coords_cols, smoke_prin$columns)
    expect_equal(ca_svd@prin_coords_rows, smoke_prin$rows)

    expect_equal(as.numeric(ca_svd@row_masses), smoke_ca$rowmass)
    expect_equal(as.numeric(ca_svd@row_inertia), smoke_ca$rowinertia)

    expect_equal(as.numeric(ca_svd@col_masses), smoke_ca$colmass)
    expect_equal(as.numeric(ca_svd@col_inertia), smoke_ca$colinertia)
})

test_that("CA coord function", {
    expect_equal(cac@std_coords_cols, smoke_ca$colcoord)
    expect_equal(cac@std_coords_rows, smoke_ca$rowcoord)

    expect_equal(cac@prin_coords_cols, smoke_prin$columns)
    expect_equal(cac@prin_coords_rows, smoke_prin$rows)
})

# cacomp test for 2x2 matrix (--> only 1 dim --> error). Error handling!
