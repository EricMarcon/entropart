EntAC <- function(Ns, q = 0, n.seq = 1:sum(Ns), Correction = "Best", 
                  NumberOfSimulations = 0, Alpha = 0.05, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  if (!is.IntValues(Ns)) {
    warning("Entropy accumulation requires integer values in argument Ns. Abundances have been rounded.")
    Ns <- round(Ns)
  }
  Size <- sum(Ns)
  
  if (Correction == "Best") Correction <- "UnveilJ"
  
  # Calculate entropy. Parallelize. Do not allow more forks.
  Entropy <- parallel::mclapply(n.seq, function(n) Tsallis(Ns, q=q, Correction=Correction, Level=n, CheckArguments=FALSE), mc.allow.recursive=FALSE)
  Entropy <- simplify2array(Entropy)
  # Simulations
  if (NumberOfSimulations > 0) {
    # Unveil the full distribution
    PsU <- as.ProbaVector(Ns, RCorrection="Jackknife", Correction="Chao2015", Unveiling="geom", CheckArguments=FALSE)
    # Create a MetaCommunity made of simulated communities
    MCSim <- rCommunity(NumberOfSimulations, size=Size, NorP=PsU, CheckArguments=FALSE)
    ProgressBar <- utils::txtProgressBar(min=0, max=NumberOfSimulations)
    Sims <- matrix(nrow=NumberOfSimulations, ncol=length(n.seq))
    # NumberOfSimulations with a progress bar
    for (i in 1:NumberOfSimulations) {
      # Parallelize. Do not allow more forks.
      ACasaList <- parallel::mclapply(n.seq, function(n) Tsallis(MCSim$Nsi[, i], q=q, Correction=Correction, Level=n, CheckArguments=FALSE), mc.allow.recursive=FALSE)
      Sims[i, ] <- simplify2array(ACasaList)
      utils::setTxtProgressBar(ProgressBar, i)
    }
    # Quantiles of simulations for each q
    EstEnvelope <- apply(Sims, 2, stats::quantile, probs = c(Alpha/2, 1-Alpha/2))
    # Simulation means to recenter the confidence envelope
    Means <- apply(Sims, 2, mean)
    # Prepare the object
    colnames(EstEnvelope) <- n.seq
    entAC <- list(x=n.seq,
                  y=Entropy,
                  low=EstEnvelope[1,] - Means + Entropy,
                  high=EstEnvelope[2,] - Means + Entropy)
  } else {
    entAC <- list(x=n.seq, y=Entropy)
  }
  
  # Format the result
  class(entAC) <- c("EntAC", "AccumCurve", class(entAC))
  # Return actual values as attributes
  attr(entAC, "Size") <- Size
  attr(entAC, "Value") <- Tsallis(Ns, q=q, Correction = "None", CheckArguments=FALSE)
  return(entAC)
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
