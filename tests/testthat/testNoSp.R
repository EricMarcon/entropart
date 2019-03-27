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
  
  testthat::expect_true(is.na(Richness(NULL)))
  testthat::expect_true(is.na(Tsallis(NULL, q=0.5)))
  testthat::expect_true(is.na(Shannon(NULL)))
  testthat::expect_true(is.na(Simpson(NULL)))
  testthat::expect_true(is.na(bcTsallis(NULL, q=0)))
  testthat::expect_true(is.na(bcTsallis(NULL, q=0.5)))
  testthat::expect_true(is.na(bcTsallis(NULL, q=1)))
  testthat::expect_true(is.na(bcTsallis(NULL, q=2)))
  testthat::expect_true(is.na(bcDiversity(NULL, q=0)))
  testthat::expect_true(is.na(bcDiversity(NULL, q=0.5)))
  testthat::expect_true(is.na(bcDiversity(NULL, q=1)))
  testthat::expect_true(is.na(bcDiversity(NULL, q=2)))
  testthat::expect_true(is.na(AbdFreqCount(NULL)))
})

# Singletons only
c111 <- rep(1, 3)
testthat::test_that("Diversity is 1 with a single species", {
  testthat::skip_on_cran()
  
  testthat::expect_equal(as.numeric(Richness(c111)), 1)
  testthat::expect_equal(as.numeric(Tsallis(2, q=0.5)), 0)
  testthat::expect_equal(as.numeric(Shannon(2)), 0)
  testthat::expect_equal(as.numeric(Simpson(2)), 0)
})
