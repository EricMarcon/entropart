Richness <-
function(NorP, ...) 
{
  UseMethod("Richness")
}


Richness.ProbaVector <-
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
  
  S <- sum(NorP > 0)
  names(S) <- "None"
  return(S)
}


Richness.AbdVector <-
function(NorP, Correction = "Best", Alpha = 0.05, JackOver = FALSE, JackMax = 10, Level = NULL, PCorrection = "Chao2015", Unveiling = "geom", RCorrection = "Rarefy", ..., CheckArguments = TRUE, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return(bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, JackMax=JackMax, CheckArguments=CheckArguments))
  } else {
    return (Richness.numeric(NorP, Correction=Correction, Level=Level, PCorrection=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, CheckArguments=CheckArguments))
  }
}


Richness.integer <-
function(NorP, Correction = "Best", Alpha = 0.05, JackOver = FALSE, JackMax = 10, Level = NULL, PCorrection = "Chao2015", Unveiling = "geom", RCorrection = "Rarefy", ..., CheckArguments = TRUE, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  if (is.null(Level)) {
    return(bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, JackMax=JackMax, CheckArguments=CheckArguments))
  } else {
    return (Richness.numeric(NorP, Correction=Correction, Level=Level, PCorrection=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, CheckArguments=CheckArguments))
  }
}


Richness.numeric <-
function(NorP, Correction = "Best", Alpha = 0.05, JackOver = FALSE, JackMax = 10, Level = NULL, PCorrection = "Chao2015", Unveiling = "geom", RCorrection = "Rarefy", ..., CheckArguments = TRUE, Ps = NULL, Ns = NULL)
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
  if (CheckArguments)
    CheckentropartArguments()
  
  if (abs(sum(NorP) - 1) < length(NorP)*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return(Richness.ProbaVector(NorP, CheckArguments=FALSE))
  } 
  if (is.null(Level)) {
    # Abundances
    return(bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, JackMax=JackMax, CheckArguments=FALSE))
  }
  # Eliminate 0
  NorP <- NorP[NorP > 0]
  N <- sum(NorP)
  # Number of observed species
  Sobs <- length(NorP)
  # If Level is coverage, get size
  if (Level < 1) 
    Level <- Coverage2Size(NorP, SampleCoverage=Level, CheckArguments=FALSE)
  if (Level == sum(NorP)) {
    # No interpolation/extrapolation needed: estimate with no correction
    return(Richness.ProbaVector(NorP/N, CheckArguments=FALSE))
  }
  if (Level <= N) {
    # Interpolation
    richness <- Sobs - sum(exp(lchoose(N-NorP, Level) - lchoose(N, Level)))
    names(richness) <- "SAC"
    return (richness)
  } else {
    # Extrapolation. 
    Singletons <-  sum(NorP == 1)
    if (Singletons) {
      # Estimate the number of unobserved species
      if (PCorrection == "None") {
        # Don't unveil the asymptotic distribution, use the asymptotic estimator
        S0 <- bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, JackMax=JackMax, CheckArguments=FALSE) - Sobs
      } else {
        # Unveil so that the estimation of richness is similar to that of non-integer entropy
        PsU <- as.ProbaVector.numeric(NorP, Correction=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, q=0, CheckArguments=FALSE)
        S0 <- length(PsU) - Sobs
      }
      richness <- Sobs + S0*(1 - (1 - Singletons/(N*S0+Singletons))^(Level-N))
    } else {
      # No singleton
      richness <- Sobs
    }
    names(richness) <- Correction
    return (richness)  
  }
}


bcRichness <- 
function(Ns, Correction = "Best", Alpha = 0.05, JackOver = FALSE, JackMax = 10, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "Best") Correction <- "Jackknife"
  
  # Eliminate 0
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)
  # Exit if Ns contains no or a single species
  if (length(Ns) < 2) {
    if (length(Ns) == 0) {
      entropy <- NA
      names(entropy) <- "No Species"
      return (entropy)
    } else {
      S <- 1
      names(S) <- "Single Species"
      return (S)
    }
  } else {
    # Probabilities instead of abundances
    if (N < 2) {
      warning("Bias correction attempted with probability data. Correction forced to 'None'")
      Correction <- "None"
    }
  }
  
  # Calculate basic statistics
  S <- length(Ns)

  # No correction
  if (Correction == "None") {
    names(S) <- "None"
    return(S)
  }

  # Calculate basic statistics
  AFC <- AbdFreqCount(Ns)
  S1 <- AFC[AFC[, 1] == 1][2]
  S2 <- AFC[AFC[, 1] == 2][2]
  
  
  # Chao1
  if ((Correction == "Chao1") | (Correction == "iChao1")) {
    if (is.na(S1)) {
      S1 <- 0
    }
    if (is.na(S2)) {
      S0 <- (N-1)/N*S1*(S1-1)/2
      S2 <- 0
    } else {
      S0 <- (N-1)/N*S1*S1/2/S2
    }
  }
  if (Correction == "Chao1") {
    S <- S + S0
    names(S) <- "Chao1"
    return(S)
  }
  if (Correction == "iChao1") {
    S3 <- AFC[AFC[, 1] == 3][2]
    S4 <- AFC[AFC[, 1] == 4][2]
    if (is.na(S3)) {
      S3 <- 0
    }
    if (is.na(S4)) {
      S4 <- 1
    }
    iS0 <- S3/4/S4*max(S1-S2*S3/2/S4, 0)
    S <- S + S0 + iS0
    names(S) <- "iChao1"
    return(S)
  }
  
  if (Correction == "Jackknife") {
    # Adapted from jackknife in SPECIES
    # Max possible order
    k <- min(nrow(AFC) - 1, JackMax)
    if (k == 0) {
      # No optimisation possible. Return "jackknife of order 0".
      Smallestk <- 1
      jackest <- S
    } else {
      # Max number of individuals
      m <- max(AFC[, 1])
      # Complete the abundance frequency count for all counts between 1 and m
      ntemp <- cbind(seq_len(m), rep(0, m))
      ntemp[AFC[, 1], 2] <- AFC[, 2]
      AFC <- ntemp
      # Prepare a matrix with k+1 rows and 5 columns
      gene <- matrix(0, nrow=k + 1, ncol=5)
      gene[1, 1] <- S
      for (i in seq_len(k)) {
        gene[i + 1, 1] <- S
        gene[i + 1, 4] <- S
        for (j in seq_len(i)) {
          gene[i + 1, 1] <- gene[i + 1, 1] + (-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) * AFC[j, 2]
          gene[i + 1, 4] <- gene[i + 1, 4] + (-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) * AFC[j, 2] * prod(seq_len(j))
        }
        gene[i + 1, 2] <- -gene[i + 1, 1]
        for (j in seq_len(i)) {
          gene[i + 1, 2] <- gene[i + 1, 2] + ((-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) + 1)^2 * AFC[j, 2]
        }
        gene[i + 1, 2] <- gene[i + 1, 2] + sum(AFC[(i + 1):nrow(AFC), 2])
        gene[i + 1, 2] <- sqrt(gene[i + 1, 2])
      }
      if (k > 1) {
        for (i in 2:k) {
          gene[i, 3] <- -(gene[i + 1, 1] - gene[i, 1])^2/(S - 1)
          for (j in seq_len(i - 1)) {
            gene[i, 3] <- gene[i, 3] + ((-1)^(j + 1) * 2^(i) * stats::dbinom(j, i, 0.5) - (-1)^(j + 1) * 2^(i - 1) * stats::dbinom(j, i - 1, 0.5))^2 * AFC[j, 2] * S/(S - 1)
          }
          gene[i, 3] <- gene[i, 3] + AFC[i, 2] * S/(S - 1)
          gene[i, 3] <- sqrt(gene[i, 3])
          gene[i, 5] <- (gene[i + 1, 1] - gene[i, 1])/gene[i, 3]
        }
      }
      # Threshold for Burnham and Overton's test
      coe <- stats::qnorm(1 - Alpha/2, 0, 1)
      # Which orders pass the test?
      x <- (gene[2:(k + 1), 5] < coe)
      if (sum(x, na.rm=TRUE) == 0) {
        # If none, keep the max value of k
        Smallestk <- k + 1
        jackest <- gene[Smallestk, 1]
        # Estimated standard error
        sej <- gene[Smallestk, 2]
      }
      else {
        # Else, keep the smallest value (+1 because jack1 is in line 2)
        Smallestk <- which(x)[1] + 1
        if (JackOver & (Smallestk < k+1))
          # Add 1 if overestimation is requested
          Smallestk <- Smallestk +1
        # Estimated value
        jackest <- gene[Smallestk, 1]
        # Estimated standard error
        sej <- gene[Smallestk, 2]
      }
    }
    names(jackest) <- paste("Jackknife", Smallestk-1)
    return(jackest) 
  }
}
