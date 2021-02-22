CommunityProfile <-
function(FUN, NorP, q.seq = seq(0, 2, 0.1), 
         NumberOfSimulations = 0, Alpha = 0.05, BootstrapMethod = "Chao2015", 
         size = 1, ..., ShowProgressBar = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments) {
    CheckentropartArguments()
  }

  # Estimated profile
  Values <- vapply(q.seq, function(q) FUN(NorP, q, ..., CheckArguments = FALSE), 0)
  
  if (NumberOfSimulations > 0) {
    if (!is.IntValues(NorP)) warning("Evaluation of the confidence interval of community profiles requires integer abundances in argument NorP. Abundances have been rounded.")
    NsInt <- round(NorP)
    # Create a MetaCommunity made of simulated communities
    if (size == 1) {
      rencenter <- TRUE
      size <- sum(NsInt)
    } else {
      rencenter <- FALSE
    }
    # The simulated communities may be of arbitrary size to obtain the confidence interval of the diversity of a smaller community
    MCSim <- rCommunity(NumberOfSimulations, size=size, NorP=NsInt, BootstrapMethod=BootstrapMethod, CheckArguments = FALSE)
    # May return NA if the bootstrap method is not recognized
    if (any(is.na(MCSim))) stop("Communities could not be simulated.")
    ProgressBar <- utils::txtProgressBar(min=0, max=NumberOfSimulations)
    Sims <- matrix(nrow=NumberOfSimulations, ncol=length(q.seq))
    # Loops are required for the progress bar, instead of:
    # Sims <- apply(MCSim$Nsi, 2, function(Nsi) CommunityProfile(FUN, Nsi, q.seq, ...)$y)
    for (i in seq_len(NumberOfSimulations)) {
      # Parallelize. Do not allow more forks in PhyloApply()
      ProfileAsaList <- parallel::mclapply(q.seq, function(q) FUN(MCSim$Nsi[, i], q, ..., CheckArguments=FALSE), mc.allow.recursive=FALSE)
      Sims[i, ] <- simplify2array(ProfileAsaList)
      if(ShowProgressBar & interactive()) 
        utils::setTxtProgressBar(ProgressBar, i)
    }
    close(ProgressBar)
    # Recenter simulated values if size is that of the community
    Means <- apply(Sims, 2, mean)
    if (rencenter) {
      Sims <- t(t(Sims)-Means+Values)
    }

    # Quantiles of simulations for each q
    EstEnvelope <- apply(Sims, 2, stats::quantile, probs = c(Alpha/2, 1-Alpha/2))
    colnames(EstEnvelope) <- q.seq
    Profile <- list(x=q.seq,
                    y=Values,
                    low=EstEnvelope[1, ],
                    high=EstEnvelope[2, ])
    if (!rencenter) {
      Profile$mid <- Means
    }
  } else {
    Profile <- list(x=q.seq,
                    y=Values)
  }
  
  class(Profile) <- "CommunityProfile"
  return (Profile)
}


as.CommunityProfile <-
function (x, y, low = NULL, high = NULL, mid = NULL) 
{
  if (!is.numeric(x))
    stop("x must be a numeric vector")
  if (!is.numeric(y))
    stop("y must be a numeric vector")
  if (length(x) != length(y))
    stop("x and y must have the same length")
  
  Profile <- list(x=x, y=y)
  if (!is.null(low)) {
    if (length(x) != length(low))
      stop("x and low must have the same length")
    Profile$low <- low
  }
  if (!is.null(high)) {
    if (length(x) != length(high))
      stop("x and high must have the same length")
    Profile$high <- high
  }
  if (!is.null(mid)) {
    if (length(x) != length(mid))
      stop("x and mid must have the same length")
    Profile$mid <- mid
  }
  class(Profile) <- "CommunityProfile"
  return(Profile)
}


is.CommunityProfile <-
function (x) 
{
  inherits(x, "CommunityProfile")
}


plot.CommunityProfile <- 
function(x, ..., main = NULL, 
         xlab = "Order of Diversity", ylab = "Diversity", ylim = NULL,
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
  
  graphics::plot(x=x$x, y=x$y, type="n", main=main, xlab=xlab, ylab=ylab, ylim=c(ymin, ymax), ...)
  CEnvelope(x, LineWidth=LineWidth, ShadeColor=ShadeColor, BorderColor=BorderColor)
}


autoplot.CommunityProfile <- 
  function(object, ..., main = NULL, 
           xlab = "Order of Diversity", ylab = "Diversity", 
           ShadeColor = "grey75", alpha = 0.3, BorderColor = "red",
           col = ggplot2::GeomLine$default_aes$colour,
           lty = ggplot2::GeomLine$default_aes$linetype,
           lwd = ggplot2::GeomLine$default_aes$size)
{  
  thePlot <- ggplot2::ggplot(as.data.frame.list(object), ggplot2::aes_(x=~x, y=~y))
  if (!(is.null(object$high) | is.null(object$low))) {
    thePlot <- thePlot +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~low, ymax=~high), fill=ShadeColor, alpha=alpha) +
      # Add red lines on borders of polygon
      ggplot2::geom_line(ggplot2::aes_(y=~low), colour=BorderColor, linetype=2) +
      ggplot2::geom_line(ggplot2::aes_(y=~high), colour=BorderColor, linetype=2)
  }
  if (!is.null(object$mid)) {
    thePlot <- thePlot +
      # Add dotted line for the mid value
      ggplot2::geom_line(ggplot2::aes_(y=~mid), linetype=2)
  }
  thePlot <- thePlot +
    ggplot2::geom_line(colour=col, linetype=lty, size=lwd) +
    ggplot2::labs(title=main, x=xlab, y=ylab)
  
  return(thePlot)
}
