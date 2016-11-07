#' Entropie zeta
#' =============
#' 
#' Fonction d'information
curve((1-x)^1, 0, 1, lty=1, xlab="p", ylab="I(p)")
curve((1-x)^2, 0, 1, lty=3, add=TRUE)
curve((1-x)^3, 0, 1, lty=3, add=TRUE)
curve((1-x)^4, 0, 1, lty=3, add=TRUE)
curve((1-x)^5, 0, 1, lty=3, add=TRUE)
curve((1-x)^25, 0, 1, lty=3, add=TRUE)

#' Contribution de chaque espèce à l'entropie
curve(x*(1-x)^1, 0, 1, lty=1, xlab="p", ylab="p I(p)")
curve(x*(1-x)^2, 0, 1, lty=3, add=TRUE)
curve(x*(1-x)^3, 0, 1, lty=3, add=TRUE)
curve(x*(1-x)^4, 0, 1, lty=3, add=TRUE)
curve(x*(1-x)^5, 0, 1, lty=3, add=TRUE)
curve(x*(1-x)^25, 0, 1, lty=3, add=TRUE)
