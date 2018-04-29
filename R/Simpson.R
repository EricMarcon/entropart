Simpson <-
function(NorP, ...) 
{
  UseMethod("Simpson")
}


Simpson.ProbaVector <-
function(NorP, ..., CheckArguments = TRUE, Ps = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      stop("An argument NorP or Ps must be provided.")
    }
  }
  if (CheckArguments)
    CheckentropartArguments()
  
  return (Tsallis.ProbaVector (NorP, q=2, CheckArguments=FALSE))
}


Simpson.AbdVector <-
function(NorP, Correction="Lande", Level = NULL, ..., CheckArguments = TRUE, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return (bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
  } else {
    return (Simpson.numeric(NorP, Correction=Correction, Level=Level, CheckArguments=CheckArguments))
  }
}


Simpson.integer <-
function(NorP, Correction="Lande", Level = NULL, ..., CheckArguments = TRUE, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return (bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
  } else {
    return (Simpson.numeric(NorP, Correction=Correction, Level=Level, CheckArguments=CheckArguments))
  }
}


Simpson.numeric <-
function(NorP, Correction="Lande", Level = NULL, ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ps or Ns must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return (Simpson.ProbaVector(NorP, CheckArguments=CheckArguments))
  } else {
    # Abundances
    if (is.null(Level)) {
      return (Simpson.AbdVector(NorP, Correction=Correction, CheckArguments=CheckArguments))
    } else {
      # Eliminate 0
      NorP <- NorP[NorP > 0]
      N <- sum(NorP)
      # If Level is coverage, get size
      if (Level < 1) Level <- Coverage2Size(NorP, SampleCoverage=Level, CheckArguments=CheckArguments)
      # Exit if Ns contains no or a single species
      if (length(NorP) < 2) {
        if (length(NorP) == 0) {
          return(NA)
        } else {
          return(0)
        }
      } else {
        entropy <- 1 - 1/Level - (1-1/Level)*sum(NorP*(NorP-1))/N/(N-1)
        names(entropy) <- "Lande"
        return (entropy)
      }
    } 
  }
}


bcSimpson <-
function(Ns, Correction="Lande", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()

  # Eliminate 0
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)
  # Exit if Ns contains no or a single species
  if (length(Ns) < 2) {
  	if (length(Ns) == 0) {
  		return(NA)
  	} else {
  		return(0)
  	}
  } else {
    # Probabilities instead of abundances
    if (N < 2) {
      warning("Bias correction attempted with probability data. Correction forced to 'None'")
      Correction <- "None"
    }
  }

  if (Correction == "Lande" | Correction == "Best") {
    entropy <- 1 - sum(Ns*(Ns-1)/N/(N-1))
    names(entropy) <- Correction
    return (entropy)
  } else {
    return (bcTsallis(Ns, q=2, Correction, CheckArguments=FALSE)) 
  }
}
