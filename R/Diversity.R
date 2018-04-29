Diversity <-
function (NorP, q = 1, ...) 
{
  UseMethod("Diversity")
}


Diversity.ProbaVector <-
function(NorP, q = 1, ..., CheckArguments = TRUE, Ps = NULL) 
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
  
  Entropy <- Tsallis(NorP, q=q, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}


Diversity.AbdVector <-
function(NorP, q = 1, Correction = "Best", Level = NULL, ..., CheckArguments = TRUE, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return (bcDiversity(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
  } else {
    Entropy <- Tsallis.AbdVector(NorP, q=q, Correction=Correction, Level=Level, CheckArguments=CheckArguments)
    return (expq(Entropy, q))  
  }
}


Diversity.integer <-
function(NorP, q = 1, Correction = "Best", Level = NULL, ..., CheckArguments = TRUE, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return (bcDiversity(Ns=NorP, q=q, Correction=Correction, CheckArguments=CheckArguments))
  } else {
    Entropy <- Tsallis.integer(NorP, q=q, Correction=Correction, Level=Level, CheckArguments=CheckArguments)
    return (expq(Entropy, q))  
  }
}


Diversity.numeric <-
function(NorP, q = 1, Correction = "Best", Level = NULL, ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return (Diversity.ProbaVector(NorP, q=q, CheckArguments=CheckArguments))
  } else {
    # Abundances
    if (is.null(Level)) {
      return (Diversity.AbdVector(NorP, q=q, Correction=Correction, Level=Level, CheckArguments=CheckArguments))
    } else {
      Entropy <- Tsallis.numeric(NorP, q=q, Correction=Correction, Level=Level, CheckArguments=CheckArguments)
      return (expq(Entropy, q))  
    }
  }
}


bcDiversity <-
function(Ns, q = 1, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- bcTsallis(Ns, q, Correction, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
