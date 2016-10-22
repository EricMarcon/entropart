context("Tsallis")

# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)

# Check Tsallis limit at 1 equals Shannon
test_that("Tsallis tends to Shannon", {
  # No correction
  expect_equal(as.numeric(Tsallis(Ps, 1 + 1E-7)),
               as.numeric(Shannon(Ps)),
               tolerance = 1e-6)
  # Best Correction
  expect_equal(as.numeric(Tsallis(Ns, 1 + 1E-7)),
               as.numeric(Shannon(Ns)),
               tolerance = 1e-6)
})
