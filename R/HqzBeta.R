HqzBeta <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), ...) 
{
  UseMethod("HqzBeta")
}


HqzBeta.ProbaVector <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), ..., CheckArguments = TRUE, Ps = NULL, Pexp = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      stop("An argument NorP or Ps must be provided.")
    }
  }
  if (missing(NorPexp)){
    if (!missing(Pexp)) {
      NorPexp <- Pexp
    } else {
      stop("An argument NorPexp or Pexp must be provided.")
    }
  }
  
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(NorP) != length(NorPexp)) {
    stop("NorP and NorPexp should have the same length.")
  }  
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(NorP))) {
    if (ncol(as.matrix(Z)) != length(NorP))  # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
  } else { # Matrix and NorP are named
    # Reorder NorP and NorPexp
    if (!setequal(names(NorP), names(NorPexp)))
      stop("NorP and NorPexp should have the same names.")
    NorPexp <- NorPexp[names(NorP)]
    if (length(setdiff(names(NorP), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(NorP), names(NorP)]
  }
  
  # Calculate (Zp)
  Zps <- Z %*% NorP
  Zpexp <- Z %*% NorPexp
  # Eliminate data when NorP is equal to zero because 0lnq(0)=0
  Zps <- Zps[NorP != 0]
  Zpexp <- Zpexp[NorP != 0]
  NorPexp <- NorPexp[NorP != 0]
  NorP <- NorP[NorP != 0]
  
  dataBeta <- NorP * (lnq(1/Zpexp, q)-lnq(1/Zps, q))
  entropy <- sum(dataBeta)
  names(entropy) <- "None"
  return (entropy)
}


HqzBeta.AbdVector <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., CheckArguments = TRUE, Ns = NULL, Nexp = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (missing(NorPexp)){
    if (!missing(Nexp)) {
      NorPexp <- Nexp
    } else {
      stop("An argument NorPexp or Nexp must be provided.")
    }
  }
  return (bcHqzBeta(Ns=NorP, Nexp=NorPexp , q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


HqzBeta.integer <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., CheckArguments = TRUE, Ns = NULL, Nexp = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (missing(NorPexp)){
    if (!missing(Nexp)) {
      NorPexp <- Nexp
    } else {
      stop("An argument NorPexp or Nexp must be provided.")
    }
  }
  return (bcHqzBeta(Ns=NorP, Nexp=NorPexp, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
}


HqzBeta.numeric <-
function(NorP, NorPexp = NULL, q = 1, Z = diag(length(NorP)), Correction = "Best", ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL, Pexp = NULL, Nexp = NULL) 
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
  if (missing(NorPexp)){
    if (!missing(Pexp)) {
      NorPexp <- Pexp
    } else {
      if (!missing(Nexp)) {
        NorP <- Nexp
      } else {
        stop("An argument NorPexp or Pexp or Nexp must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return (HqzBeta.ProbaVector(NorP, NorPexp, q=q, Z=Z, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return (HqzBeta.AbdVector(NorP, NorPexp, q=q, Z=Z, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcHqzBeta <-
function(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ns) != length(Nexp)) {
    stop("Ns and Nexp should have the same length.")
  }  
  
  # No correction available yet
  if (Correction == "None" | Correction == "Best") {
    return (HqzBeta.ProbaVector(Ns/sum(Ns), Nexp/sum(Nexp), q, Z, CheckArguments=FALSE))
  }
  
  warning("Correction was not recognized")
  return (NA)
}
