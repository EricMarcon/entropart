GammaEntropy <-
function(MC, q = 1, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  if (!is.null(Tree)) {
    # Method <- "HCDT"
    if (is.IntValues(MC$Ns) | Correction=="None") {
      # Integer abundances or no correction: just estimate.
      Entropy <- bcPhyloEntropy(MC$Ns, q, Tree, Normalize, Correction, CheckArguments=FALSE)$Total
    } else {
      # Use the fallback estimation
      SampleCoverage <- Coverage(rowSums(MC$Nsi), CheckArguments=FALSE)
      Entropy <- bcPhyloEntropy(MC$Ns, q, Tree, Normalize, Correction, SampleCoverage=SampleCoverage, CheckArguments=FALSE)$Total
    }
  } else {
    if (!is.null(Z)) {
      # Method <- "Similarity-based"
      if (is.IntValues(MC$Ns) | Correction=="None") {
        # Integer abundances or no correction: just estimate.
        Entropy <- bcHqz(MC$Ns, q, Z, Correction, CheckArguments=FALSE)
      } else {
        # Use the fallback estimation
        SampleCoverage <- Coverage(rowSums(MC$Nsi), CheckArguments=FALSE)
        Entropy <- bcHqz(MC$Ns, q, Z, Correction, SampleCoverage=SampleCoverage, CheckArguments=FALSE)
      }
    } else {
      # Method <- "HCDT"
      if (is.IntValues(MC$Ns) | Correction=="None") {
        # Integer abundances or no correction: just estimate.
        Entropy <- bcTsallis(MC$Ns, q, Correction, CheckArguments=FALSE)
      } else {
        # Use the fallback estimation
        SampleCoverage <- Coverage(rowSums(MC$Nsi), CheckArguments=FALSE)
        Entropy <- bcTsallis(MC$Ns, q, Correction, SampleCoverage=SampleCoverage, CheckArguments=FALSE)
      }
    }
  }
  
  return(Entropy)
}
