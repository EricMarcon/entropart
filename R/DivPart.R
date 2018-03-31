DivPart <-
function(q = 1, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree. Height is 1 by default, including species neutral diversity.
  Height <- 1
  ppTree <- Preprocess.Tree(Tree)
  if (!is.null(ppTree$Height) & !Normalize) {
    Height <- ppTree$Height
  }
  
  # Correction methods
  if (Biased) Correction = "None"
  if (!is.IntValues(MC$Ns)) {
    # Only ChaoShen, Marcon (HCDT only, not Similarity-based) or None are acceptable if MC$Ns are not integers.
    if (Correction == "Best") {
      if (is.null(Z)) {
        Correction <- "Marcon"
      } else {
        Correction <- "ChaoShen"
      }
    }
  }

  # Alpha and beta entropy of communities. PhyloDetails=TRUE forces the return of a is.PhyloValue if phyloentropy is calculated.
  GammaE <- GammaEntropy(MC, q, Correction, ppTree, Normalize, Z, PhyloDetails=TRUE, CheckArguments=FALSE)
  # Alpha entropy is estimated with the same correction as Gamma
  if (is.PhyloValue(GammaE)) {
    # Phyloentropy. GammaEntropy returned a PhyloValue with Corrections along the tree
    AlphaCorrection <- GammaE$Corrections
    # Make GammaE a number again
    GammaE <- GammaE$Total
  } else {
    # GammaEntropy is a named number
    AlphaCorrection <- names(GammaE)
  }
  AlphaE <- AlphaEntropy(MC, q, AlphaCorrection, ppTree, Normalize, Z, CheckArguments=FALSE)
  # beta is calculated as gamma-alpha to ensure continuity. Community beta entropy is not calculated.
  BetaE  <- list(Communities = NA, Total = GammaE - AlphaE$Total)      

  # Total Diversities
  AlphaD <- expq(AlphaE$Total / Height, q) * Height
  GammaD <- expq(GammaE / Height, q) * Height
  BetaD  <- GammaD / AlphaD
    # equals: expq(BetaE$Total / Height / (1 - (q-1)*AlphaE$Total/Height), q) 
  
  DivPart <- (list(
    MetaCommunity = ArgumentOriginalName(MC),
    Order = q, 
    Biased = Biased, 
    Correction = AlphaCorrection,
    Normalized = Normalize,
    TotalAlphaDiversity = AlphaD,
    TotalBetaDiversity = BetaD,
    GammaDiversity = GammaD,
    CommunityAlphaDiversities = expq(AlphaE$Communities / Height, q) * Height,
    TotalAlphaEntropy = AlphaE$Total,
    TotalBetaEntropy = BetaE$Total,
    GammaEntropy = GammaE,
    CommunityAlphaEntropies = AlphaE$Communities,
    CommunityBetaEntropies = BetaE$Communities
    ))
  if(!is.null(Tree))
    DivPart$Tree <- ArgumentOriginalName(Tree)
  if(is.null(Z)) {
    DivPart$Method <- "HCDT"
  } else {
    DivPart$Method <- "Similarity-based"
    DivPart$Z <- ArgumentOriginalName(Z)
  }
  class(DivPart) <- "DivPart"
  
  return (DivPart)
}


is.DivPart <-
function (x) 
{
  inherits(x, "DivPart")
}


plot.DivPart <- 
function (x, ...) 
{
  graphics::plot(c(0, x$GammaDiversity), c(0, length(x$CommunityAlphaDiversities)), type = "n", xlab = expression(paste(alpha, " and ", gamma, " diversity")), ylab = expression(paste(beta, " diversity")), ...)
  graphics::rect(0, 0, x$GammaDiversity, 1, lty=2)
  graphics::rect(0, 0, x$TotalAlphaDiversity, x$TotalBetaDiversity, lty=2)
}



autoplot.DivPart <- 
function (object, ...) 
{
  thePlot <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes_(xmin=~xmin, ymin=~ymin, xmax=~xmax, ymax=~ymax), 
                       data.frame(xmin=0, ymin=0, xmax=object$GammaDiversity, ymax=1), alpha=0.3) +
    ggplot2::geom_rect(ggplot2::aes_(xmin=~xmin, ymin=~ymin, xmax=~xmax, ymax=~ymax), 
                       data.frame(xmin=0, ymin=0, xmax=object$TotalAlphaDiversity, ymax=object$TotalBetaDiversity), alpha=0.3) +
    ggplot2::expand_limits(y=c(0, length(object$CommunityAlphaDiversities))) +
    ggplot2::labs(x=expression(paste(alpha, " and ", gamma, " diversity")), y=expression(paste(beta, " diversity"))) +
    ggplot2::annotate(geom="text", y=0, x=object$TotalAlphaDiversity, parse=TRUE, label=as.character(expression(alpha))) +
    ggplot2::annotate(geom="text", y=0, x=object$GammaDiversity, parse=TRUE, label=as.character(expression(gamma))) +
    ggplot2::annotate(geom="text", x=0, y=object$TotalBetaDiversity, parse=TRUE, label=as.character(expression(beta)))
  
  return(thePlot)
}


summary.DivPart <-
function(object, ...) 
{    
  cat(object$Method, "diversity partitioning of order", object$Order, "of metaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated\naccording to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated\naccording to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  
  cat("Alpha diversity of communities:", "\n")
  print(object$CommunityAlphaDiversities)
  cat("Total alpha diversity of the communities:", "\n")
  print(object$TotalAlphaDiversity)
  cat("Beta diversity of the communities:", "\n")
  print(object$TotalBetaDiversity)
  cat("Gamma diversity of the metacommunity:", "\n")
  print(object$GammaDiversity)
  
  return(invisible(NULL))
}