testthat::context("Less than 2 species")


# 1 species
testthat::test_that("Diversity is 1 with a single species", {
  testthat::skip_on_cran()
  
  testthat::expect_equal(as.numeric(Richness(2)), 1)
  testthat::expect_equal(as.numeric(Tsallis(2, q=0.5)), 0)
  testthat::expect_equal(as.numeric(Shannon(2)), 0)
  testthat::expect_equal(as.numeric(Simpson(2)), 0)
})

# 0 species
testthat::test_that("Diversity is NA with species", {
  testthat::skip_on_cran()
  
  c0 <- numeric(0)
  testthat::expect_true(is.na(Richness(c0)))
  testthat::expect_true(is.na(Tsallis(c0, q=0.5)))
  testthat::expect_true(is.na(Shannon(c0)))
  testthat::expect_true(is.na(Simpson(c0)))
  testthat::expect_true(is.na(bcTsallis(c0, q=0)))
  testthat::expect_true(is.na(bcTsallis(c0, q=0.5)))
  testthat::expect_true(is.na(bcTsallis(c0, q=1)))
  testthat::expect_true(is.na(bcTsallis(c0, q=2)))
  testthat::expect_true(is.na(bcDiversity(c0, q=0)))
  testthat::expect_true(is.na(bcDiversity(c0, q=0.5)))
  testthat::expect_true(is.na(bcDiversity(c0, q=1)))
  testthat::expect_true(is.na(bcDiversity(c0, q=2)))
  testthat::expect_true(is.na(AbdFreqCount(c0)))
})

# Singletons only
testthat::test_that("Diversity can be calculated when coverage is zero", {
  testthat::skip_on_cran()
  
  c111 <- rep(1, 3)
  testthat::expect_equal(as.numeric(Richness(c111)), 3)
  testthat::expect_silent(Richness(c111))
  testthat::expect_equal(as.numeric(Tsallis(c111, q=0.5)), 1.464102, tolerance = 10^-6)
  testthat::expect_warning(Tsallis(c111, q=0.5))
  testthat::expect_equal(as.numeric(Shannon(c111)), 1.098612, tolerance = 10^-6)
  testthat::expect_warning(Shannon(c111))
  testthat::expect_equal(as.numeric(Simpson(c111)), 1) # Lande estimator.
  testthat::expect_silent(Simpson(c111))
})
