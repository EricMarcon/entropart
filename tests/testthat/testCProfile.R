testthat::context("Community Profile")

# Improve coverage
testthat::test_that("Diversity profiles can be plotted", {
  testthat::skip_on_cran()
  
  data(Paracou618)
  Ns <- Paracou618.MC$Ns
  
  testthat::expect_type(CP <- CommunityProfile(Diversity, Ns, NumberOfSimulations=20), "list")
  # plot and autoplot
  testthat::expect_silent(plot(CP))
  testthat::expect_silent(autoplot(CP))

  # Check as.CommunityProfile
  testthat::expect_silent(asCP <- as.CommunityProfile(CP$x, CP$y, CP$low, CP$high))
  testthat::expect_equivalent(asCP, CP)
})
