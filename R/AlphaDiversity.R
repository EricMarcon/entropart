AlphaDiversity <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()

  # Preprocess the tree to get its height
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  
  # Calculate normalized entropy, Height will be addressed later
  AlphaEntropy <- AlphaEntropy(MC, q, Correction, ppTree, Z, Normalize=TRUE)
  Diversity <- list(
    MetaCommunity = ArgumentOriginalName(MC),
    Method = AlphaEntropy$Method,
    Type = "alpha",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Communities = expq(AlphaEntropy$Communities, q) * Height,
    Total = expq(AlphaEntropy$Total, q)* Height
    )
  if(!is.null(Tree))
    Diversity$Tree <- ArgumentOriginalName(Tree)
  if(!is.null(Z))
    Diversity$Z <- ArgumentOriginalName(Z)
  class(Diversity) <- "MCdiversity"
  
  return(Diversity)  
}
