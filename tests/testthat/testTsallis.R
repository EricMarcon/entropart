testthat::context("Tsallis")

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


# Improve coverage
testthat::test_that("Tsallis is computed at arbitrary levels", {
  testthat::skip_on_cran()
  
  # Interpolation
  testthat::expect_named(Tsallis(Ns, q=1.5, Level = .5), "Interp")
  # Extrapolation by rarefaction
  testthat::expect_named(Tsallis(Ns, q=1.5, Level = sum(Ns)+1), "Rarefy")
})


# Improve coverage
testthat::test_that("Various estimators are available", {
  testthat::skip_on_cran()
  
  # Shannon
  testthat::expect_named(Shannon(Ns, Correction = "ChaoJost"), "ChaoJost")
  testthat::expect_named(Shannon(Ns, Correction = "ZhangHz"), "ZhangHz")
  testthat::expect_named(Shannon(Ns, Correction = "Holste"), "Holste")
  testthat::expect_named(Shannon(Ns, Correction = "Bonachela"), "Bonachela")
  testthat::expect_named(Shannon(Ns, Correction = "GenCov"), "GenCov")
  testthat::expect_named(Shannon(Ns, Correction = "Marcon"), "Marcon")
  testthat::expect_named(Shannon(Ns, Correction = "Grassberger"), "Grassberger")
  testthat::expect_named(Shannon(Ns, Correction = "Grassberger2003"), "Grassberger2003")
  testthat::expect_named(Shannon(Ns, Correction = "Schurmann"), "Schurmann")
  testthat::expect_named(Shannon(Ns, Correction = "UnveilC"), "UnveilC")
  testthat::expect_named(Shannon(Ns, Correction = "UnveiliC"), "UnveiliC")

  # Simpson
  testthat::expect_named(Simpson(Ns, Correction = "ChaoJost"), "ChaoJost")
  
  # Tsallis
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "ChaoJost"), "ChaoJost")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "ZhangGrabchak"), "ZhangGrabchak")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "ChaoShen"), "ChaoShen")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "GenCov"), "GenCov")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "Marcon"), "Grassberger")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "Grassberger"), "Grassberger")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "UnveilC"), "UnveilC")
  testthat::expect_named(Tsallis(Ns, q=1.5, Correction = "UnveiliC"), "UnveiliC")
})

