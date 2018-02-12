context("DivProfile")

# Check DivProfile 
testthat::test_that("DivProfile with simulations is correct", {
  testthat::skip_on_cran()
  
  data(Paracou618)
  # Calculate DivProfile
  testthat::expect_true(is.DivProfile(dp <- DivProfile(q.seq = 0:2, Paracou618.MC, NumberOfSimulations = 10)))
  # Plot it
  testthat::expect_silent(plot(dp))
})
