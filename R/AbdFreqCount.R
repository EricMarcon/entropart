AbdFreqCount <- 
function (Ns, Level = NULL, Estimator = "Best", CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Convert into integers
  NsInt <- as.integer(round(Ns))
  if (any(abs(NsInt-Ns) > sum(Ns)*.Machine$double.eps))
    warning ("The abundance frequency count requires integer abundances. Abundances have been rounded.")
  
  # Eliminate 0
  Ns <- NsInt[NsInt > 0]
  # Calculate abundance distribution
  DistNs <- tapply(Ns, Ns, length)
  
  if (is.null(Level)) {
    # No extrapolation. Prepare a two-column matrix
    afc <- matrix(c(as.integer(names(DistNs)), DistNs), ncol = 2)
  } else {
    Singletons <- DistN["1"]
    SampleSize <- sum(Ns)
    # If Level is coverage, get size
    if (Level < 1) Level <- Coverage2Size(Ns, Level, CheckArguments=FALSE)
    if (Level <= SampleSize) {
      # Interpolation
      Snu <- sapply(1:Level, function(nu) sum(exp(lchoose(Ns, nu) + lchoose(SampleSize-Ns, Level-nu) - lchoose(SampleSize, Level))))
    } else {
      # Extrapolation. Unveiled estimator currently the best.
      if (Estimator == "Best") Estimator <- "UnveilJ"
      if (Estimator == "UnveilJ") {
        # Unveil the full distribution
        PsU <- as.ProbaVector(Ns, RCorrection="Jackknife", Correction="Chao2015", Unveiling="geom")
        # Extrapolate
        Snu <- sapply(1:Level, function(nu) sum(exp(lchoose(Level, nu) + nu*log(PsU) + (Level-nu)*log(1-PsU))))
      }
      # Make a matrix and eliminate zeros
      afc <- cbind(1:Level, Snu)[Snu > .Machine$double.eps, ]
    }
  }
  colnames(afc) <- c("Abundance", "NbSpecies")
  class(afc) <- c("AbdFreqCount", class(afc))
  return(afc)
}
