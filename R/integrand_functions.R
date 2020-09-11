dirichlet_for_integration <- function(x, alpha){
  out <- gtools::ddirichlet(as.vector(x), alpha)/sqrt(length(alpha)) # normalizing constant so that is sums to 1 on the unit simplex
  ifelse(is.nan(out) | is.infinite(out), 0, out)
}

logisticnormal_for_integration <- function(x, mu, sigma){
  out <- dlogisticnormal(as.vector(x), mu, sigma)/sqrt(length(mu)) # normalizing constant so that is sums to 1 on the unit simplex
  ifelse(is.nan(out) | is.infinite(out), 0, out)
}
