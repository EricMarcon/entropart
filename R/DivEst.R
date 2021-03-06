DivEst <-
function(q = 0, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, Simulations = 100, ShowProgressBar = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }

  # Estimation from data
  RealEst <- DivPart(q, MC, Biased, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)

  # RedrawSpecies resamples a community according to species abundances.
  RedrawSpecies<- function(SpeciesAbundances){
    # Very simplified (for speed) version of rCommunity with BootstrapMethod="Marcon"
    stats::rmultinom(1, sum(SpeciesAbundances), SpeciesAbundances)
  }

  # SimulateEntropy resamples all communities and calculates entropy
  SimulateEntropy <- function(Progression) {
    SimNsi <- apply(MC$Nsi, 2, RedrawSpecies)
    rownames(SimNsi) <- rownames(MC$Nsi) 
    SimMC <- Preprocess.MC(SimNsi, MC$Wi)
    NewSim <- DivPart(q, SimMC, Biased, Correction, Tree, Normalize, Z, CheckArguments=FALSE)
    # update progress bar
    if(ShowProgressBar & interactive()) 
      utils::setTxtProgressBar(ProgressBar, Progression)
    c(NewSim$TotalAlphaEntropy, NewSim$TotalBetaEntropy, NewSim$GammaEntropy)
  }
  
  # Simulate entropy
  ProgressBar <- utils::txtProgressBar(min=0, max=Simulations)
  RawSimulatedEntropy <- vapply(seq_len(Simulations), SimulateEntropy, FUN.VALUE=rep(0.0, 3))
  close(ProgressBar)
  
  # Recenter entropy
  SimulatedEntropy <- RawSimulatedEntropy + 
    with(RealEst, c(TotalAlphaEntropy, TotalBetaEntropy, GammaEntropy)) - 
    apply(RawSimulatedEntropy, 1, mean)
  # Transform entropy to diversity
  if (q == 1) { 
    SimulatedDiversity <- exp(SimulatedEntropy)
  } else {
    SimulatedDiversity <- SimulatedEntropy
    SimulatedDiversity[1,] <- expq(SimulatedEntropy[1,] / Height, q) * Height
    SimulatedDiversity[3,] <- expq(SimulatedEntropy[3,] / Height, q) * Height
    SimulatedDiversity[2,] <- SimulatedDiversity[3,] / SimulatedDiversity[1,]
      # equals: expq(SimulatedEntropy[2,] / Height / (1 - (q-1)*SimulatedEntropy[1,]/Height), q) * Height
  }
  dimnames(SimulatedDiversity)[[1]] <- list("Alpha", "Beta", "Gamma")

  DivEst <- RealEst
  DivEst$SimulatedDiversity <- SimulatedDiversity
  class(DivEst) <- c("DivEst", class(RealEst))

  return (DivEst)
}


is.DivEst <-
function (x) 
{
  inherits(x, "DivEst")
}


plot.DivEst <- 
  function (x, ..., main = NULL, Which = "All", 
            Quantiles = c(0.025, 0.975), colValue = "red", lwdValue = 2, ltyValue = 2,
            colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2) 
  {
    # Save graphical parameters
    op <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow=c(2, 2))
    
    if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Alpha Diversity"
    if (Which == "All" | Which == "Alpha") {
      graphics::plot(as.SimTest(x$TotalAlphaDiversity, x$SimulatedDiversity["Alpha",]), main=main, 
                     Quantiles=Quantiles, ..., colValue=colValue, lwdValue=lwdValue, ltyValue=ltyValue,
                     colQuantiles=colQuantiles, lwdQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles) 
    }
    if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
    if (Which == "All" | Which == "Beta") {
      graphics::plot(as.SimTest(x$TotalBetaDiversity, x$SimulatedDiversity["Beta",]), main=main, 
                     Quantiles=Quantiles, ..., colValue=colValue, lwdValue=lwdValue, ltyValue=ltyValue,
                     colQuantiles=colQuantiles, lwdQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles)
    }
    if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
    if (Which == "All" | Which == "Gamma") {
      graphics::plot(as.SimTest(x$GammaDiversity, x$SimulatedDiversity["Gamma",]), main=main, 
                     Quantiles=Quantiles, ..., colValue=colValue, lwdValue=lwdValue, ltyValue=ltyValue,
                     colQuantiles=colQuantiles, lwdQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles)
    }
    
    # Legend and restore parameters
    if (Which == "All") {
      graphics::par(mar=c(0, 0, 0, 0))
      graphics::plot(0:10, 0:10, type="n", xlab=NULL, frame.plot=FALSE, xaxt="n", yaxt="n", col.lab="white")
      leg <- c("Null Distribution", "True Estimate", "95% confidence interval") 
      graphics::legend(2, 8, leg, col = c(1, 2, 1), lty = 1:3, merge = TRUE, cex=1)
      graphics::par(op)
    }
  }


autoplot.DivEst <- 
function (object, ..., main = NULL, Which = "All", labels = NULL, font.label = list(size=11, face="plain"),
          Quantiles = c(0.025, 0.975), colValue = "red", colQuantiles = "black", ltyQuantiles = 2) 
{
  if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Alpha Diversity"
  if (Which == "All" | Which == "Alpha") {
    AlphaPlot <- autoplot(as.SimTest(object$TotalAlphaDiversity, object$SimulatedDiversity["Alpha",]), main=main, 
                          Quantiles=Quantiles, ..., colValue=colValue, colQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles)
  }
  if (Which == "Alpha")
    return(AlphaPlot)

  if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
  if (Which == "All" | Which == "Beta") {
    BetaPlot <- autoplot(as.SimTest(object$TotalBetaDiversity, object$SimulatedDiversity["Beta",]), main=main, 
                         Quantiles=Quantiles, ..., colValue=colValue, colQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles)
  }
  if (Which == "Beta")
    return(BetaPlot)
  
  if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
  if (Which == "All" | Which == "Gamma") {
    GammaPlot <- autoplot(as.SimTest(object$GammaDiversity, object$SimulatedDiversity["Gamma",]), main=main, 
                          Quantiles=Quantiles, ..., colValue=colValue, colQuantiles=colQuantiles, ltyQuantiles=ltyQuantiles)
  }
  if (Which == "Gamma")
    return(GammaPlot)
  
  # Which == "All": return a multiple plot
  return(ggpubr::ggarrange(AlphaPlot + ggplot2::labs(x=NULL), BetaPlot +  ggplot2::labs(x=NULL), GammaPlot, 
                           ncol = 1, nrow = 3, labels=labels, font.label=font.label))
}


summary.DivEst <-
function(object, ...) 
{
  cat("Diversity partitioning of order", object$Order, "of MetaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  cat("Alpha diversity of communities:", "\n")
  print(object$CommunityAlphaDiversities)
  cat("Total alpha diversity of the communities:", "\n")
  print(object$TotalAlphaDiversity)
  cat("Beta diversity of the communities:", "\n")
  print(object$TotalBetaDiversity)
  cat("Gamma diversity of the metacommunity:", "\n")
  print(object$GammaDiversity)
  
  cat("Quantiles of simulations (alpha, beta and gamma diversity):\n")
  quant <- c(0, 0.01, 0.025, 0.05, 0.1, seq(0.25, 0.75, 0.25), 0.9, 0.95, 0.975, 0.99, 1)
  print(stats::quantile(object$SimulatedDiversity["Alpha", ], probs = quant))
  print(stats::quantile(object$SimulatedDiversity["Beta", ], probs = quant))
  print(stats::quantile(object$SimulatedDiversity["Gamma", ], probs = quant))
  
  return(invisible(NULL))
}