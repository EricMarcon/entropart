testthat::context("Hurlbert")

data(Paracou618)

# Check Hurlbert diversity
testthat::test_that("Hurlbert diversity is correct", {
  testthat::skip_on_cran()

  # Order 2
  testthat::expect_equal(as.numeric(HurlbertD(Paracou618.MC$Ps, 2)), as.numeric(expq(Simpson(Paracou618.MC$Ps), 2)))
  testthat::expect_lt(bcHurlbertD(Paracou618.MC$Ns, 2) - Diversity(Paracou618.MC$Ns, q=2), 1)
})

