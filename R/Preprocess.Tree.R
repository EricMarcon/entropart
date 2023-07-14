Preprocess.Tree <-
function(Tree)
{
  # The tree may be NULL or already processed
  if (is.null(Tree) | inherits(Tree, "PPtree")) return (Tree)

  # Tree must be either a phylog, phylo or a hclust object
  if (inherits(Tree, "phylog")) {
    # Build an hclust object to use cutree later. 
    # Distances in $Wdist are actually sqrt(2*distance)
    # Caution: Distances in hclust count full branch lengths between species, i.e. twice the ultramtetric distance
    # See ?as.phylo.hclust
    hTree <- stats::hclust(Tree$Wdist^2/2, "average")
    # build a phylo object
    phyTree <- ape::as.phylo.hclust(hTree)
    # Double edge.lengths to correct as.phylo.hclust
    phyTree$edge.length <- 2*phyTree$edge.length
  } else {
    if (inherits(Tree, "phylo")) {
      phyTree <- Tree
      # Build an hclust object to use cutree later.
      # Edge lengths are multiplied by 2 during the conversion. Divide by 2 before that.
      Tree$edge.length <- Tree$edge.length/2
      # Make the tree binary if it contains polytomies or the conversion fails
      if (!ape::is.binary(Tree)) {
        Tree <- ape::multi2di(Tree)
      }
      hTree <- ape::as.hclust.phylo(Tree)
    } else {
      if (inherits(Tree, "hclust")) {
        # Caution: Distances in hclust count full branch lengths between species, i.e. twice the ultramtetric distance
        hTree <- Tree
        # build a phylo object to use $droot later
        phyTree <- ape::as.phylo.hclust(Tree)
        # Double edge.lengths to correct as.phylo.hclust
        phyTree$edge.length <- 2*phyTree$edge.length
      } else {
        stop("Tree must be an object of class phylo, phylog or hclust")
      }
    }
  }

  # Calculate distances between nodes and leaves
  DistancesFromLeaves <- ape::branching.times(phyTree)
  # Get a sorted list of cuts (eliminate leaves)
  Cuts <- sort(DistancesFromLeaves[setdiff(names(DistancesFromLeaves), names(phyTree$leaves))])
  # Calculate intervals between cuts (add 0 to Cuts to get the first interval)
  Intervals <- diff(c(0, Cuts))
  # Eliminate 0 intervals (when a node contains more than 2 tips), including rounding errors
  RoundingError <- max(DistancesFromLeaves)*10*.Machine$double.eps
  Cuts <- Cuts[Intervals > RoundingError]
  Intervals <- Intervals[Intervals > RoundingError]

  ppTree <- list(
    phyTree   = phyTree,
    hTree     = hTree,
    Height    = Cuts[length(Cuts)],
    Cuts      = Cuts,
    Intervals = Intervals
    )
  class(ppTree) <- "PPtree"
  return (ppTree)
}


# Internal function tips mimics geiger::tips (geiger is no longer maintained)
# Much slower than the geiger::tips
# Arguments are 
# Tree : a phylo object
# node : an internal node number of the phylo object
# Returns a character vector : the names of the tips descending from the node
tips <- function(Tree, node) 
{
  nTips <- length(Tree$tip.label)
  if (node > nTips) {
    # internal node numbers start after the last tip
    return(ape::extract.clade(Tree, node)$tip.label)
  } else {
    return(Tree$tip.label[node])
  }
}
