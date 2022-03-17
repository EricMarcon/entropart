is.MCdiversity <-
function (x) 
{
  inherits(x, "MCdiversity")
}


plot.MCdiversity <- 
function (x, ...) 
{
  graphics::barplot(c(x$Communities, 0, x$Total),
          beside = TRUE,
          width = c(x$Weights, .5, .1),
          names.arg = c(names(x$Communities), "", "Metacommunity"),
          ylab = "Diversity",
          ...
  )
  
}


autoplot.MCdiversity <- 
function (object, col = ggplot2::GeomCol$default_aes$fill,
          border = ggplot2::GeomCol$default_aes$colour, ...) 
{
  theData <- data.frame(x=c(names(object$Communities), "", "Metacommunity"), y=c(object$Communities, 0, object$Total))
  # Factors to keep the order of bars (avoid sort by name)
  theData$x <- factor(theData$x, levels=theData$x)
                      
  thePlot <- ggplot2::ggplot(theData, ggplot2::aes(x=.data$x, y=.data$y)) + 
    ggplot2::geom_col(width=c(object$Weights, .5, .1), fill=col, colour=border) +
    ggplot2::labs(y="Diversity")
  
  return(thePlot)
}


summary.MCdiversity <-
function(object, ...) 
{
  cat(object$Method, object$Type, "diversity of order", object$Order, "of metaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n", fill=TRUE)
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated according to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  if (!is.null(object$Communities)) {
    cat("Diversity of communities:", "\n")
    print(object$Communities)
  }
  cat("Average diversity of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}