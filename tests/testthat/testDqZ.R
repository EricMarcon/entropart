testthat::context("DqZ")

# Check DqZ 
testthat::test_that("DqZ with Z=I equals Diversity", {
  testthat::skip_on_cran()
  
  data(Paracou618)
  # DqZ and Diversity
  testthat::expect_equal(Dqz(Paracou618.MC$Ps), Diversity(Paracou618.MC$Ps))
})
