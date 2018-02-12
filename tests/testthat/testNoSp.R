testthat::context("Less than 2 species")


# 1 species
testthat::test_that("Diversity is 1 with a single species", {
  testthat::skip_on_cran()
  
  testthat::expect_equal(as.numeric(Shannon(2)), 0)
  testthat::expect_equal(as.numeric(Simpson(2)), 0)
  testthat::expect_equal(as.numeric(Tsallis(2)), 0)
  testthat::expect_equal(as.numeric(Tsallis(2, q=2)), 0)
  testthat::expect_equal(as.numeric(Diversity(2)), 1)
  testthat::expect_equal(as.numeric(Diversity(2), q=2), 1)
})

# 0 species
testthat::test_that("Diversity is NA with species", {
  testthat::skip_on_cran()
  
  testthat::expect_true(is.na(bcShannon(NULL)))
  testthat::expect_true(is.na(bcSimpson(NULL)))
  testthat::expect_true(is.na(bcTsallis(NULL)))
  testthat::expect_true(is.na(bcTsallis(NULL, q=2)))
  testthat::expect_true(is.na(bcDiversity(NULL)))
  testthat::expect_true(is.na(bcDiversity(NULL, q=2)))
})
