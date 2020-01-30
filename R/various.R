### Find how much of variance to explain
findWeights <- function(target.vec, var.explained, max.weights) {
  weights                           <- (target.vec)^2 / sum((target.vec)^2)
  sort.weight                       <- sort(weights, decreasing = TRUE, index.return = TRUE)
  cut.off                           <- min(which.min(cumsum(sort.weight$x) < var.explained), max.weights)
  return(list('ind'     = sort.weight$ix[1:cut.off],
              'cut.off' = cut.off))
}
