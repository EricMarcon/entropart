is.MCentropy <-
function (x)
{
  inherits(x, "MCentropy")
}


plot.MCentropy <- 
function (x, ...) 
{
  graphics::barplot(c(x$Communities, 0, x$Total),
          beside = TRUE,
          width = c(x$Weights, .5, 1),
          names.arg = c(names(x$Communities), "", "Metacommunity"),
          ylab = "Entropy",
          ...
  ) 
}


autoplot.MCentropy <- 
function (object, col = ggplot2::GeomCol$default_aes$fill,
          border = ggplot2::GeomCol$default_aes$colour, ...) 
{
  theData <- data.frame(x=c(names(object$Communities), "", "Metacommunity"), y=c(object$Communities, 0, object$Total))
  # Factors to keep the order of bars (avoid sort by name)
  theData$x <- factor(theData$x, levels=theData$x)
  
  thePlot <- ggplot2::ggplot(theData, ggplot2::aes(x=.data$x, y=.data$y)) + 
    ggplot2::geom_col(width=c(object$Weights, .5, .1), fill=col, colour=border) +
    ggplot2::labs(y="Entropy")
  
  return(thePlot)
}



summary.MCentropy <-
function(object, ...) 
{
  cat(object$Method, object$Type, "entropy of order", object$Order, "of metaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n", fill=TRUE)
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Entropy is", ifelse(object$Normalized, "normalized", "not normalized"), fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated according to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  cat("Entropy of communities:", "\n")
  print(object$Communities)
  cat("Average entropy of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}