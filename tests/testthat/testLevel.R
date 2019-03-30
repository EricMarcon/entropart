testthat::context("Interpolation and extrapolation")

data(Paracou618)
Ns <- Paracou618.MC$Ns
N <- sum(Ns)

# Check Continuity
testthat::test_that("Interpolation and extrapolation are continuous", {
  testthat::skip_on_cran()

  # value at N is the average of N-1 and N+1
  testthat::expect_equal(as.numeric((Diversity(Ns, q=0, Level=N-1) + Diversity(Ns, q=0, Level=N+1))/2), 
                         as.numeric(Diversity(Ns, q=0, Level = N)), tolerance=Diversity(Ns, q=0, Level = N)/100)
  testthat::expect_equal(as.numeric((Diversity(Ns, q=1, Level=N-1) + Diversity(Ns, q=1, Level=N+1))/2), 
                         as.numeric(Diversity(Ns, q=1, Level = N)), tolerance=Diversity(Ns, q=1, Level = N)/100)
  testthat::expect_equal(as.numeric((Diversity(Ns, q=2, Level=N-1) + Diversity(Ns, q=2, Level=N+1))/2), 
                         as.numeric(Diversity(Ns, q=2, Level = N)), tolerance=Diversity(Ns, q=2, Level = N)/100)
  testthat::expect_equal(as.numeric((Diversity(Ns, q=1.5, Level=N-1) + Diversity(Ns, q=1.5, Level=N+1))/2), 
                         as.numeric(Diversity(Ns, q=1.5, Level = N)), tolerance=Diversity(Ns, q=1.5, Level = N)/100)
})


# Check Coverage to size
testthat::test_that("Size can be obtained from coverage", {
  testthat::skip_on_cran()
  
  # Coverage2Size finds the actual size
  testthat::expect_equal(as.numeric(Coverage2Size(Ns, Coverage(Ns))), N)
  # Coverage2Size returns 1 for very small coverage
  testthat::expect_equal(as.numeric(Coverage2Size(Ns, 10^-5)), 1)
})


# Improve coverage
testthat::test_that("Accumulation can be plotted", {
  testthat::skip_on_cran()
  
  # DivAC: interpolation, extrapolation, bootstrap, q=0, 1, 2 and 1.5
  testthat::expect_type(DAC <- DivAC(Ns, q=0, n.seq = c(1, N, 2*N), NumberOfSimulations = 20), "list")
  testthat::expect_type(DAC <- DivAC(Ns, q=1, n.seq = c(1, N, 2*N), NumberOfSimulations = 20), "list")
  testthat::expect_type(DAC <- DivAC(Ns, q=2, n.seq = c(1, N, 2*N), NumberOfSimulations = 20), "list")
  testthat::expect_type(DAC <- DivAC(Ns, q=1.5, n.seq = c(1, N, 2*N), NumberOfSimulations = 20), "list")
  # plot and autoplot
  testthat::expect_silent(plot(DAC))
  testthat::expect_silent(autoplot(DAC))

  # Check as.AccumCurve
  testthat::expect_silent(asDAC <- as.AccumCurve(DAC$x, DAC$y, DAC$low, DAC$high))
  testthat::expect_equivalent(asDAC, DAC)
  
  # AbdFreqCount
  testthat::expect_silent(AbdFreqCount(Ns, Level = 2*N))
})
