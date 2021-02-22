is.PhyloValue <-
function (x) 
{
  inherits(x, "PhyloValue")
}


plot.PhyloValue <- 
function (x, xlab = expression(italic("T")), ylab = NULL, main = NULL, ...) 
{
  Entity <- ""
  # Entity
  if (is.PhyloEntropy(x)) {
    Entity <- "Entropy"
  } else {
    if (is.PhyloDiversity(x)) {
      Entity <- "Diversity"
    }
  }
  
  # ylab
  if (is.null(ylab))
    ylab <- Entity
  
  # main
  if (is.null(main))
    main <- paste(Entity, "along the tree")
  
  graphics::plot(x=c(0, names(x$Cuts)), y=c(x$Cuts, 0), type="s", xlab=xlab, ylab=ylab, main=main, ...)
  graphics::abline(h=x$Total, lty=2)
}


autoplot.PhyloValue <- 
function (object, xlab = expression(italic("T")), ylab = NULL, main = NULL, 
          col = ggplot2::GeomLine$default_aes$colour,
          lty = ggplot2::GeomLine$default_aes$linetype,
          lwd = ggplot2::GeomLine$default_aes$size,
          ...) 
{
  Entity <- ""
  # Entity
  if (is.PhyloEntropy(object)) {
    Entity <- "Entropy"
  } else {
    if (is.PhyloDiversity(object)) {
      Entity <- "Diversity"
    }
  }
  
  # ylab
  if (is.null(ylab))
    ylab <- Entity
  
  # main
  if (is.null(main))
    main <- paste(Entity, "along the tree")
  
  # Prepare data. as.numeric() to avoid factors.
  df <- data.frame(Time=as.numeric(c(0, rep(names(object$Cuts), each=2))), Value=c(rep(object$Cuts, each=2), 0))
  
  # Plot
  thePlot <- ggplot2::ggplot(data=df, ggplot2::aes_(x=~Time, y=~Value)) +
    ggplot2::geom_line(colour=col, linetype=lty, size=lwd) +
    ggplot2::labs(title=main, x=xlab, y=ylab) +
    ggplot2::geom_hline(yintercept=object$Total, linetype=2)
  
  return(thePlot)
}


summary.PhyloValue <-
function(object, ...) 
{
  cat(object$Function, "applied to", object$Distribution, "along the tree:", object$Tree, fill=TRUE)
  cat("\nResults are", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  cat("\nThe average value is:", object$Total)
  cat("\n\nValues along the tree are:\n")
  print(object$Cuts)
  
  return(invisible(NULL))
}