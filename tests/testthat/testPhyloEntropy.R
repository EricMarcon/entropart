context("PhyloEntropy")

# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Ns is the total number of trees per species
Ns <- as.AbdVector(Paracou618.MC$Ns)
# Species probabilities
Ps <- as.ProbaVector(Paracou618.MC$Ns)


# Check PhyloEntropy of order 2 equals Rao
test_that("PhyloEntropy of order 2 equals Rao", {
  # No correction
  expect_equal(as.numeric(PhyloEntropy(Ps, 2, Paracou618.Taxonomy, Normalize = FALSE)$Total), 
               as.numeric(Rao(Ps, Paracou618.Taxonomy)))
  # Best correction
  expect_equal(as.numeric(PhyloEntropy(Ns, 2, Paracou618.Taxonomy, Normalize = FALSE)$Total), 
               as.numeric(Rao(Ns, Paracou618.Taxonomy)))
})
