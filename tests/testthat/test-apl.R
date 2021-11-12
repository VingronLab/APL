
# tab <- read.delim(file = "/home/kohl/PhD/gits/APL/tests/testthat/testdata/input_data.tsv")
# mat <- as.matrix(tab[,-1])
# rownames(mat) <- tab$Country.Name
# save(mat, file = "/home/kohl/PhD/gits/APL/tests/testthat/testdata/countries.rda")

load("./testdata/countries.rda")

ca <- cacomp(mat, princ_coords = 3)
ca <- apl_coords(ca, group = c(6, 7, 8, 10, 12))

test_that("check apl coordinates", {

})
