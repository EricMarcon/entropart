as.SimTest <- 
function (RealValue, SimulatedValues) {
  st <- list(RealValue=RealValue, SimulatedValues=SimulatedValues)
  class(st) <- "SimTest"
  return(st)
}


is.SimTest <-
function (x) 
{
  inherits(x, "SimTest")
}


plot.SimTest <- 
function (x, Quantiles = c(0.025, 0.975), ..., colValue = "red", lwdValue = 2, ltyValue = 2, 
          colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2, 
          main = NULL, xlab = "Simulated Values", ylab = "Density")
{
  plot(stats::density(x$SimulatedValues), main=main, xlab=xlab, ylab=ylab, ...)
  graphics::abline(v=x$RealValue, col=colValue, lwd=lwdValue, lty=ltyValue)
  for (qt in Quantiles) {
    graphics::abline(v=stats::quantile(x$SimulatedValues, probs = qt), col=colQuantiles, lwd=lwdQuantiles, lty=ltyQuantiles)
  }
}


autoplot.SimTest <- 
function (object, Quantiles = c(0.025, 0.975), ..., colValue = "red", colQuantiles = "black", ltyQuantiles = 2, 
          main = NULL, xlab = "Simulated Values", ylab = "Density")
{
  df <- data.frame(SimulatedValues=object$SimulatedValues)
  thePlot <- ggplot2::ggplot() +
    ggplot2::geom_density(data=df, ggplot2::aes(x=.data$SimulatedValues, fill=factor("unique"), alpha=0.8)) +
    ggplot2::labs(title=main, x=xlab, y=ylab) +
    ggplot2::geom_vline(xintercept=object$RealValue, colour=colValue)
  for (qt in Quantiles) {
    thePlot <- thePlot +
      ggplot2::geom_vline(xintercept=stats::quantile(object$SimulatedValues, probs = qt), colour=colQuantiles, linetype=ltyQuantiles)
  }
  return(thePlot+ggplot2::theme(legend.position="none"))
}


summary.SimTest <-
function(object, Quantiles = c(0.025, 0.975), ...) 
{ 
  cat("Real value: ", object$RealValue, "\n")
  cat("Quantile in the simulated distribution: ", stats::ecdf(object$SimulatedValues)(object$RealValue), "\n")
  
  cat("Quantiles of simulations:\n")
  for (qt in Quantiles) {
    cat(sprintf("%1.2f%%", 100*qt), ": ", stats::quantile(object$SimulatedValues, probs = qt), "\n")
  }
  cat("Mean simulated value: ", mean(object$SimulatedValues), "\n")
  return(invisible(NULL))
}
