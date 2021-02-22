PhyloEntropy <-
function(NorP, q = 1, Tree, Normalize = TRUE, ...) 
{
  UseMethod("PhyloEntropy")
}


PhyloEntropy.ProbaVector <-
function(NorP, q = 1, Tree, Normalize = TRUE, ..., CheckArguments = TRUE, Ps = NULL) 
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
  
  # Calculate the PhyloValue
  Entropy <- PhyloApply(Tree, Tsallis, NorP, Normalize, q=q, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- ArgumentOriginalName(NorP)
  Entropy$Tree <- ArgumentOriginalName(Tree)
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  Entropy$Correction <- names(Entropy$Total) <- "None"
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)
}


PhyloEntropy.AbdVector <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", ..., CheckArguments = TRUE, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloEntropy(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloEntropy.integer <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", ..., CheckArguments = TRUE, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcPhyloEntropy(Ns=NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
}


PhyloEntropy.numeric <-
function(NorP, q = 1, Tree, Normalize = TRUE, Correction = "Best", ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return(PhyloEntropy.ProbaVector(NorP, q=q, Tree=Tree, Normalize=Normalize, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(PhyloEntropy.AbdVector(NorP, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcPhyloEntropy <-
function(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", SampleCoverage = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # If SampleCoverage is a vector, prepare an argument dataframe for PhyloApply
  if(is.null(SampleCoverage)) {
    dfArgs <- NULL
  } else {
    dfArgs <- data.frame(SampleCoverage=SampleCoverage)
  }
  # If Correction is a vector, idem
  if(length(Correction) > 1) {
    if (is.null(dfArgs)) {
      # Create a new dataframe...
      dfArgs <- data.frame(Correction=Correction)
    } else {
      # ... or a a column
      dfArgs <- cbind(dfArgs, data.frame(Correction=Correction))    
    }
  }
  
  # Calculate the PhyloValue
  if(length(Correction) == 1) {
    # Call PhyloApply with an argument Correction
    Entropy <- PhyloApply(Tree, FUN=bcTsallis, NorP=Ns, Normalize=Normalize, dfArgs=dfArgs, q=q, Correction=Correction, CheckArguments=FALSE)
  } else {
    # Call PhyloApply without an argument Correction since it is in dfArgs
    Entropy <- PhyloApply(Tree, FUN=bcTsallis, NorP=Ns, Normalize=Normalize, dfArgs=dfArgs, q=q, CheckArguments=FALSE)
  }
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- ArgumentOriginalName(Ns)
  Entropy$Tree <- ArgumentOriginalName(Tree)
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  # Corrections. May be a vector.
  Entropy$Correction <- Correction
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)                         
}


is.PhyloEntropy <-
function (x) 
{
  inherits(x, "PhyloEntropy")
}


summary.PhyloEntropy <-
function(object, ...) 
{
  cat(object$Type, "phylogenetic or functional entropy of order", object$Order, "of distribution", object$Distribution, fill=TRUE)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Entropy is", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  }
  cat("\nEntropy equals:", object$Total, "\n")
  return(invisible(NULL))
}