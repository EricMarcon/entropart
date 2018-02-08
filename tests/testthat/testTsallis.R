context("Tsallis")

# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)

# Check Tsallis limit at 1 equals Shannon
testthat::test_that("Tsallis tends to Shannon", {
  testthat::skip_on_cran()
# No correction
  testthat::expect_equal(as.numeric(Tsallis(Ps, 1 + 1E-7)),
               as.numeric(Shannon(Ps)),
               tolerance = 1e-6)
  # Best Correction
  testthat::expect_equal(as.numeric(Tsallis(Ns, 1 + 1E-7)),
               as.numeric(Shannon(Ns)),
               tolerance = 1e-6)
})

# Check Shannon vs EntropyEstimation
Ns <- Ns[Ns>0]
N <- sum(Ns)
testthat::test_that("Shannon with Zhang-Grabchak correction equals EntropyEstimation::Entropy.z", {
  testthat::skip_on_cran()
  # No correction
  testthat::expect_equal(sum(Ns/N*(digamma(N)-digamma(Ns))),
               EntropyEstimation::Entropy.z(Ns),
               tolerance = 1e-6)
})
