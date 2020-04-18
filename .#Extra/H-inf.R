# Estimation de la plus petite probabilité
EstimateHinf <- function(SampleSize=5000) {
  # Taille de la communauté
  cSize <- 1E6
  # Tirage d'une grande communauté log-normale
  Community <- rCommunity(1, size=cSize)
  plot(as.ProbaVector(Community))
  # Entropie d'ordre -infini
  x <- 1/min(as.ProbaVector(Community))
  # Echantillonnage de quelques hectares
  Sample <- rmultinom(1, size=SampleSize, Community)
  # Communauté échantillonnée en bleu
  lines(x=1:length(Sample), y=sort(as.ProbaVector(Sample), decreasing = TRUE), col="blue")
  # Reconstitution de la communauté selon Chao et al., 2015 (estimateur de S modifié)
  Unveiled <- as.ProbaVector(Sample, Correction = "Chao2015", Unveiling = "geom", RCorrection = "Jackknife")
  # Communauté reconstituée en rouge
  lines(x=1:length(Unveiled), y=sort(Unveiled, decreasing = TRUE), col="red")
  # Entropie d'ordre -infini estimée
  y <- 1/min(Unveiled)
  # Retour des deux valeurs pour comparaison
  return(c(x,y))
}

# n simulations, taille de l'échantillon en paramètre.
replicate(n=10, EstimateHinf(1000))
