testthat::context("Diversity of order 2")

# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)

# Check Hurlbert equals Simpson
testthat::test_that("Hurlbert equals Simpson", {
  testthat::skip_on_cran()
  # No correction
  testthat::expect_equal(as.numeric(Hurlbert(Ps)-1), as.numeric(Simpson(Ps)), tolerance = 1e-8, scale = Simpson(Ps))
})


# Check GenSimpson
testthat::test_that("GenSimpson equals Simpson", {
  testthat::skip_on_cran()

  testthat::expect_equal(as.numeric(GenSimpson(Ps, 1)), as.numeric(Simpson(Ps)))
  testthat::expect_lt(bcGenSimpson(Ns, 1) - bcSimpson(Ns), 1/1000)
})

