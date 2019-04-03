EntAC <- function(Ns, q = 0, n.seq = 1:sum(Ns), PCorrection = "Chao2015", Unveiling = "geom", RCorrection = "Rarefy", 
                  NumberOfSimulations = 0, Alpha = 0.05, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  if (!is.IntValues(Ns)) {
    warning("Entropy accumulation requires integer values in argument Ns. Abundances have been rounded.")
    Ns <- round(Ns)
  }
  N <- sum(Ns)
  
  # Prepare the vector of results
  Entropy <- numeric(length(n.seq))
  ProgressBar <- utils::txtProgressBar(min=0, max=length(n.seq))
  # i must be initialized if the accumulation contains extrapolation only
  i <- 0
  
  # Interpolation
  n.seqInt <- n.seq[n.seq < N]
  # Calculate entropy at each level
  for(Level in n.seqInt) {
    # Calculate Entropy
    i <- which(n.seq==Level)
    Entropy[i] <- Tsallis.numeric(Ns, q=q, Level=Level, PCorrection=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, CheckArguments=FALSE)
    if(interactive()) utils::setTxtProgressBar(ProgressBar, i)
  }
  # Level == Sample Size
  if (any(n.seq==N)) {
    i <- which(n.seq==N)
    Entropy[i] <- Tsallis.ProbaVector(Ns/N, q=q, CheckArguments=FALSE)
    if(interactive()) utils::setTxtProgressBar(ProgressBar, i)
  }
  # Extrapolation. Don't use Tsallis for speed.
  n.seqExt <- n.seq[n.seq > N]
  PsU <- NULL
  if (length(n.seqExt)) {
    # Unveil the full distribution that rarefies to the observed entropy (or other options)
    PsU <- as.ProbaVector.numeric(Ns, Correction=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, q=q, CheckArguments=FALSE)
    # Richness
    if (q == 0) {
      Singletons <-  sum(Ns == 1)
      if (Singletons) {
        # Estimate the number of unobserved species
        Sobs <- sum(Ns > 0)
        S0 <- length(PsU) - Sobs
        # Extrapolate richness (the vector is n.seqExt)
        Entropy[(i+1):length(n.seq)] <- Sobs + S0*(1 - (1 - Singletons/(N*S0+Singletons))^(n.seqExt-N)) -1
      } else {
        # No singleton
        Entropy[(i+1):length(n.seq)] <- Sobs -1
      }
      if(interactive()) utils::setTxtProgressBar(ProgressBar, length(n.seq))
    } else {
      # Shannon
      if (q == 1) {
        # Estimate the asymptotic entropy
        Hinf <- Shannon.ProbaVector(PsU, CheckArguments=FALSE)
        # Estimate observed entropy
        Hn <- Shannon.ProbaVector(Ns/N, CheckArguments=FALSE)
        # Interpolation (the vector is n.seqExt)
        Entropy[(i+1):length(n.seq)] <- N/n.seqExt*Hn + (n.seqExt-N)/n.seqExt*Hinf
        if(interactive()) utils::setTxtProgressBar(ProgressBar, length(n.seq))
      } else {
        # Simpson
        if (q == 2) {
          # Exit if Ns contains no or a single species
          if (length(Ns) < 2) {
            if (length(Ns) == 0) {
              Entropy[(i+1):length(n.seq)] <- NA
            } else {
              Entropy[(i+1):length(n.seq)] <- 0
            }
          } else {
            # Valid extrapolation (the vector is n.seqExt)
            Entropy[(i+1):length(n.seq)] <- 1 - 1/n.seqExt - (1-1/n.seqExt)*sum(Ns*(Ns-1))/N/(N-1)
          }
          if(interactive()) utils::setTxtProgressBar(ProgressBar, length(n.seq))
        } else {
          # General case: q is not 0, 1 or 2 
          for(Level in n.seqExt) {
            # Abundance frequence count at Level (Chao et al., 2014, eq. 5)
            Snu <- sapply(1:Level, function(nu) sum(exp(lchoose(Level, nu) + nu*log(PsU) + (Level-nu)*log(1-PsU))))
            # Estimate entropy (Chao et al., 2014, eq. 6)
            i <- which(n.seq==Level)
            Entropy[i]  <- (sum(((1:Level)/Level)^q * Snu) - 1) / (1-q)
            if(interactive()) utils::setTxtProgressBar(ProgressBar, i)
          }
        }
      }
    }
  }

  
  # Simulations: generate distributions from the unveiled probabilities
  if (NumberOfSimulations > 0 & (PCorrection=="None" | Unveiling=="None")) {
    warning("Accumulation confidence interval can't be estimated without unveiling the asymptotic distribution. Neither PCorrection nor Unveiling can be 'None'")
    NumberOfSimulations <- 0
  }
  if (NumberOfSimulations > 0) {
    # Prepare the result matrix
    Envelope <- matrix(0.0, nrow = length(n.seq), ncol = 2)
    if (is.null(PsU)) {
      # Unveil the full distribution if not done before
      PsU <- as.ProbaVector.numeric(Ns, Correction=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, q=q, CheckArguments=FALSE)
    }
    for(Level in n.seq) {
      # Generate simulated communities at each level
      Communities <- stats::rmultinom(NumberOfSimulations, size=Level, prob=PsU)
      # Probabilities
      Communities <- Communities/Level
      # Calculate entropy
      Entropies <- apply(Communities, 2, Tsallis.ProbaVector, q=q, CheckArguments=FALSE)
      i <- which(n.seq==Level)
      # Store quantiles
      Envelope[i, ] <- stats::quantile(Entropies, probs = c(Alpha/2, 1-Alpha/2))
      if(interactive()) utils::setTxtProgressBar(ProgressBar, i)
    }
    entAC <- list(x=n.seq, y=Entropy, low=Envelope[, 1], high=Envelope[, 2])
  } else {
    entAC <- list(x=n.seq, y=Entropy)
  }
  
  close(ProgressBar)
  # Format the result
  class(entAC) <- c("EntAC", "AccumCurve", class(entAC))
  # Return actual values as attributes
  attr(entAC, "Size") <- N
  attr(entAC, "Value") <- Tsallis(Ns, q=q, Correction="None", CheckArguments=FALSE)
  return(entAC)
}


DivAC <- function(Ns, q = 0, n.seq = 1:sum(Ns), PCorrection = "Chao2015", Unveiling = "geom", RCorrection = "Rarefy", 
                  NumberOfSimulations = 0, Alpha = 0.05, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()

  # Estimate entropy
  Accumulation <- EntAC(Ns=Ns, q=q, n.seq=n.seq, PCorrection=PCorrection, Unveiling=Unveiling, RCorrection=RCorrection, NumberOfSimulations=NumberOfSimulations, Alpha=Alpha, CheckArguments=FALSE)
  # Calculate diversity
  Accumulation$y <- expq(Accumulation$y, q)
  if (!is.null(Accumulation$low))
    Accumulation$low <- expq(Accumulation$low, q)
  if (!is.null(Accumulation$high))
    Accumulation$high <- expq(Accumulation$high, q)
  
  # Change the class
  class(Accumulation)[which(class(Accumulation)=="EntAC")] <- "DivAC"
  attr(Accumulation, "Value") <- expq(attr(Accumulation, "Value"), q)
  return(Accumulation)
}


as.AccumCurve <-
function (x, y, low = NULL, high = NULL) 
{
  if (!is.numeric(x))
    stop("x must be a numeric vector")
  if (!is.numeric(y))
    stop("y must be a numeric vector")
  if (length(x) != length(y))
    stop("x and y must have the same length")
  
  AccumCurve <- list(x=x, y=y)
  if (!is.null(low)) {
    if (length(x) != length(low))
      stop("x and low must have the same length")
    AccumCurve$low <- low
  }
  if (!is.null(high)) {
    if (length(x) != length(high))
      stop("x and high must have the same length")
    AccumCurve$high <- high
  }
  class(AccumCurve) <- "AccumCurve"
  return(AccumCurve)
}


is.AccumCurve <-
function (x) 
{
  inherits(x, "AccumCurve")
}



plot.AccumCurve <- 
function(x, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, ylim = NULL,
         LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
{  
  if (is.null(ylim)) {
    # Evaluate ylim if not set by an argument
    if (is.null(x$low)) {
      ymin <- min(x$y)
    } else {
      ymin <- min(x$low)
    }
    if (is.null(x$high)) {
      ymax <- max(x$y)
    } else {
      ymax <- max(x$high)
    }
  } else {
    ymin <- ylim[1]
    ymax <- ylim[2]
  }
  
  if (is.null(ylab)) {
    if (inherits(x, "EntAC")) {
      ylab <- "Entropy"
    } else {
      if (inherits(x, "DivAC")) {
        ylab <- "Diversity"
      }
    }
  }
  
  graphics::plot(x=x$x, y=x$y, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(ymin, ymax), ...)
  CEnvelope(x, LineWidth=LineWidth, ShadeColor=ShadeColor, BorderColor=BorderColor)
  
  # Actual value
  graphics::abline(v=attr(x, "Size"), lty=2)
  graphics::abline(h=attr(x, "Value"), lty=2)
}



autoplot.AccumCurve <- 
function(object, ..., main = NULL, 
         xlab = "Sample Size", ylab = NULL, 
         ShadeColor = "grey75", alpha = 0.3, BorderColor = "red")
{  
  thePlot <- ggplot2::ggplot(as.data.frame.list(object), ggplot2::aes_(x=~x, y=~y))
  if (!(is.null(object$high) | is.null(object$low))) {
    thePlot <- thePlot +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~low, ymax=~high), fill=ShadeColor, alpha=alpha) +
      # Add red lines on borders of polygon
      ggplot2::geom_line(ggplot2::aes_(y=~low), colour=BorderColor, linetype=2) +
      ggplot2::geom_line(ggplot2::aes_(y=~high), colour=BorderColor, linetype=2)
  }
  if (is.null(ylab)) {
    if (inherits(object, "EntAC")) {
      ylab <- "Entropy"
    } else {
      if (inherits(object, "DivAC")) {
        ylab <- "Diversity"
      }
    }
  }
  thePlot <- thePlot +
    ggplot2::geom_line() +
    ggplot2::labs(main=main, x=xlab, y=ylab)
  
  # Actual value
  thePlot <- thePlot +
    ggplot2::geom_hline(yintercept = attr(object, "Value"), lty=2) +
    ggplot2::geom_vline(xintercept = attr(object, "Size"), lty=2)
    
  return(thePlot)
}
