Coverage <-
function(Ns, Estimator = "Best", Level = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Round values
  Ns <- as.integer(round(Ns))
  # Eliminate zeros
  Ns <- Ns[Ns>0]
  # Calculate abundance distribution
  DistN <- tapply(Ns, Ns, length)
  # Singletons. Convert named number to number.
  Singletons <- as.numeric(DistN["1"])
  SampleSize <- sum(Ns)
  
  if (is.null(Level)) {
    # Etimate C at the observed level
    
    if (Estimator == "Best") Estimator <- "ZhangHuang"
    # More accurate
    
    # No singletons, C=1
    if (is.na(Singletons)) {
      coverage <- 1
      names(coverage) <- "No singleton"
      return(coverage)
    }
    
    # Singletons only
    if (Singletons == SampleSize) {
      warning ("Sample coverage is 0, most bias corrections will return NaN.")
      coverage <- 0
      names(coverage) <- "Singletons only"
      return(coverage)
    }
  
    if (Estimator == "ZhangHuang") {
      Ps <- Ns/SampleSize
      if (any(Ps >= .5)) {
        warning ("Zhang-Huang sample coverage cannot be estimated because one probability is over 1/2. Chao estimator is returned.")
        Estimator <- "Chao"
      } else {
        Nu <- as.integer(names(DistN))
        # Use Nu%%2*2-1 for (-1)^(Nu+1)
        coverage <- 1 - sum((Nu%%2*2-1) / choose(SampleSize, Nu) * DistN)
        names(coverage) <- Estimator
        return(coverage)
      }    
    }
    if (Estimator == "Chao") {
      coverage <- 1 - Singletons / SampleSize * (1-ChaoA(Ns))
      names(coverage) <- Estimator
      return(coverage)
    }
    if (Estimator == "Turing") {
      coverage <- 1 - Singletons / SampleSize
      names(coverage) <- Estimator
      return(coverage)
    }
    
  } else {
    # Chose level. Must be an integer. CheckEntropartArguments() may have accepted a value between 0 and 1
    if (Level <=1) stop("Level must be an integer >1.")

    if (Estimator == "Best") Estimator <- "Chao"
    # Extrapolation allowed.
    
    if (Estimator == "Good") {
      if (Level >= SampleSize) stop("The Good estimator only allows interpolation: Level must be less than the observed community size.")
      coverage <- 1 - EntropyEstimation::GenSimp.z(Ns, Level)
      names(coverage) <- Estimator
      return(coverage)
    }
    if (Estimator == "Chao") {
      if (Level < SampleSize) {
        # Interpolation
        NsRestricted <- Ns[(SampleSize - Ns) >= Level]
        coverage <- 1 - sum(NsRestricted/SampleSize 
          * exp(lgamma(SampleSize - NsRestricted + 1) - lgamma(SampleSize - NsRestricted - Level + 1) - lgamma(SampleSize) + lgamma(SampleSize - Level)))
      } else {
        # Extrapolation
        if (is.na(Singletons)) {
          # No singletons, C=1
          coverage <- 1
          names(coverage) <- "No singleton"
          return(coverage)
        } else {
          coverage <- 1 - Singletons / SampleSize * (1 - ChaoA(Ns))^(Level - SampleSize + 1)
        }
      }
      names(coverage) <- Estimator
      return(coverage)
    }
  }
  
  warning("Correction has not been recognized")
  return(NA)
}



# Helper for Chao's estimator. Not exported.
# A's formula depends on the presence of singletons and doubletons.
# Ns must be a vector of positive integers (not checked).
ChaoA <- function(Ns) {
  # Calculate abundance distribution
  DistN <- tapply(Ns, Ns, length)
  Singletons <- as.numeric(DistN["1"])
  Doubletons <- as.numeric(DistN["2"])
  SampleSize <- sum(Ns)
  
  if (is.na(Singletons)) {
    A <- 0
  } else {
    if (is.na(Doubletons)) {
      Chao1S0 <- (SampleSize - 1) * Singletons * (Singletons - 1) / 2 / SampleSize
    } else {
      Chao1S0 <- (SampleSize - 1) * Singletons^2 / 2 / SampleSize / Doubletons
    }
    A <- 1- SampleSize * Chao1S0/(SampleSize * Chao1S0 + Singletons)
  }
  
  return(A)
}



Coverage2Size <-
function(Ns, SampleCoverage, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Round values
  Ns <- as.integer(round(Ns))
  # Eliminate zeros
  Ns <- Ns[Ns>0]
  # Calculate abundance distribution
  DistN <- tapply(Ns, Ns, length)
  # Singletons. Convert named number to number.
  Singletons <- as.numeric(DistN["1"])
  SampleSize <- sum(Ns)
  
  # Singletons only
  if (Singletons == SampleSize) {
    stop("Sample coverage is 0.")
  }
  
  # Actual coverage
  C <- Coverage(Ns, CheckArguments = FALSE)
  
  if (SampleCoverage >= C) {
    # Extrapolation
    Size <- round(SampleSize + (log(SampleSize/Singletons) + log(1-SampleCoverage))/log(1-ChaoA(Ns)) - 1)
  } else {
    # Interpolation. Numeric resolution: minimize Delta.
    Delta <- function(Size) abs(Coverage(Ns, Estimator = "Chao", Level=Size, CheckArguments = FALSE) - SampleCoverage)
    Size <- round(stats::optimize(Delta, lower=1, upper=SampleSize)$minimum)
  }
  return(Size)
}