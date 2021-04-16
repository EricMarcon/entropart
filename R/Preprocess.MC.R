Preprocess.MC <-
function(Nsi, Wi)
{
  # n_i
  Ni <- colSums(Nsi)
  if (any(Ni == 0))
    stop("All communities must contain at least one individual. Remove empty communities before attempting to build a metacommunity.")
  # p_si
  Psi <- Nsi %*% diag(1/Ni) 
  dimnames(Psi) <- dimnames(Nsi)
  # p_s
  Ps <- apply(Psi %*% diag(Wi), 1, sum)
  # n is by choice the sum of ni. Avoid integer overflow.
  N <- sum(as.numeric(Ni))
  # n_s
  Ns <- Ps * N
  
  MC <- list(
    Nsi = Nsi,
    Ns = Ns,
    Ni = Ni,
    N = N,
    Psi = Psi,
    Wi = Wi,
    Ps = Ps,
    Nspecies = dim(Nsi)[1],
    Ncommunities = dim(Nsi)[2],
    SampleCoverage = Coverage(rowSums(Nsi)),
    SampleCoverage.communities = apply(Nsi, 2, Coverage, CheckArguments=FALSE) 
    )
  class(MC) <- "MetaCommunity"
  return (MC)
}
