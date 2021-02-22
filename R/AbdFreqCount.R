AbdFreqCount <- 
function (Ns, Level = NULL, PCorrection="Chao2015", Unveiling="geom", RCorrection="Rarefy", CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Convert into integers
  if (length(Ns) == 0) {
    # Ns may contain no speices. Return NA
    return(NA)
  } else {
    NsInt <- as.integer(round(Ns))    
  }
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
    SampleSize <- sum(Ns)
    # If Level is coverage, get size
    if (Level < 1) Level <- Coverage2Size(Ns, SampleCoverage=Level, CheckArguments=FALSE)
    if (Level <= SampleSize) {
      # Interpolation
      Snu <- vapply(seq_len(Level), function(nu) sum(exp(lchoose(Ns, nu) + lchoose(SampleSize-Ns, Level-nu) - lchoose(SampleSize, Level))), FUN.VALUE=0.0)
      # Make a matrix with all possible abundances
      afc <- cbind(seq_len(Level), Snu)
      # Return the estimator as an attribute
      attr(afc, "Estimator") <- "Interp"
    } else {
      # Extrapolation. 
      if (length(Ns) == 1) {
        # Single species: general formula won't work: log(1-PsU)
        Snu <- c(rep(0, Level-1), 1)
      } else {
        # Unveil the full distribution
        PsU <- as.ProbaVector(Ns, Correction=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, CheckArguments=FALSE)
        # Extrapolate
        Snu <- vapply(seq_len(Level), function(nu) sum(exp(lchoose(Level, nu) + nu*log(PsU) + (Level-nu)*log(1-PsU))), FUN.VALUE=0.0)
      }  
      # Make a matrix with all possible abundances
      afc <- cbind(seq_len(Level), Snu)
      # Return the estimator as an attribute
      attr(afc, "Estimator") <- "Extrap"
    }
  }
  colnames(afc) <- c("Abundance", "NbSpecies")
  class(afc) <- c("AbdFreqCount", class(afc))
  return(afc)
}
