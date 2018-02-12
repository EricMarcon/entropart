testthat::context("Coverage")

# Check bad arguments
testthat::test_that("Warning is returned", {
  testthat::skip_on_cran()
  # Singletons only
  testthat::expect_warning(Coverage(rep(1,5)), 
                           "Sample coverage is 0, most bias corrections will return NaN.", 
                           ignore.case = TRUE)
  # Zhang-Huang
  testthat::expect_warning(Coverage(c(8, 4, 2, 1), Estimator="ZhangHuang"), 
                           "Zhang-Huang sample coverage cannot be estimated because one probability is over 1/2. Chao estimator is returned.", 
                           ignore.case = TRUE)
  # Chao no doubletons
  testthat::expect_warning(Coverage(c(1, 5), Estimator="Chao"), 
                           "Chao's sample coverage cannot be estimated because there are no doubletons. Turing estimator is returned.", 
                           ignore.case = TRUE)
  
})

# Check rarely used estimators
testthat::test_that("Coverage is estimated", {
  testthat::skip_on_cran()
  # Chao
  testthat::expect_lt(abs(Coverage(1:5, Estimator="Chao")-Coverage(1:5)), 1/1000)
  # Turing
  testthat::expect_lt(abs(Coverage(1:5, Estimator="Turing")-Coverage(1:5)), 1/100)
})