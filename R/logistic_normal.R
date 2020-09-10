#' Functions for the logistic normal distribution
#'
#' Functions to compute the density of or generate random deviates from the
#' logistic normal distribution.
#'
#' The logistic normal distribution arises when we draw from a multivariate
#' normal distribution and then apply the logistic transformation to the resulting
#' draws. It has been described as the Gaussian on the unit simplex. It allows
#' for more flexible forms of dependence among the components
#' than the better-known Dirichlet distribution.
#'
#' In this parameterization, exactly one element of the random vector must be fixed
#' at zero, i.e. mean zero with zero variance and zero covariance with other
#' elements of the vector. Other parameterizations leave this element out of the
#' multivariate draw. The reason for including it here is that we can then reshuffle
#' the parameters as we cycle through e.g. pairs of candidates.
#'
#' @param n Number of random vectors to generate
#' @param mu Vector of expectations (pre-logistic transformation).
#' One value must be zero, and it must correspond to the row/column of
#' \code{sigma} that is also
#' all zeros.
#' @param sigma Matrix of covariances (pre-logistic transformation).
#' There must be a \code{j} such that \code{sigma[j,]} and \code{sigma[,j]}
#' are all zeros, and \code{mu[j]} must also be zero.
#' @param x A vector containing one random deviate (i.e. one )
#' @param tol \code{dlogisticnormal} returns zero if the sum of
#' \code{x} is further than \code{tol} away from 1.
#'
#' @examples
#' mu <- c(.1, -.1, 0)
#' sigma <- rbind(c(.25, .05, 0),
#'                c(.05, .15, 0),
#'                c(0, 0, 0))
#' rlogisticnormal(10, mu = mu, sigma = sigma)
#' dlogisticnormal(rbind(c(.4, .35, .25), c(.3, .5, .2)), mu = mu, sigma = sigma)
#'
#'
#' @name logisticnormal
NULL

#' @rdname logisticnormal
#' @export
rlogisticnormal <- function(n, mu, sigma){
  zero_index <- get_zero_index_and_check_errors(mu, sigma)
  #draw from multivariate normal
  draws <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
  ed <- exp(draws) # logit transformation
  ed/apply(ed, 1, sum)
}

#' @rdname logisticnormal
#' @export
dlogisticnormal <- function(x, mu, sigma, tol = 1e-10){
  zero_index <- get_zero_index_and_check_errors(mu, sigma)
  if(!is.matrix(x)){
    if (is.data.frame(x)){
      x <- as.matrix(x)
    }else{x <- t(x)}
  }
  mu_mat <- matrix(mu, nrow = nrow(x), ncol = length(mu), byrow = T)
  stopifnot(ncol(x) == ncol(mu_mat))
  stopifnot(ncol(mu) == ncol(sigma))
  stopifnot(ncol(sigma) == nrow(sigma))

  # make versions that separate out the zero and the others
  mu_mat_1 <- mu_mat[,-zero_index]
  x_1 <- matrix(x[,-zero_index], nrow = nrow(x))
  x_0 <- matrix(x[,zero_index], nrow = nrow(x_1), ncol = ncol(x_1), byrow = F)
  sigma_1 <- sigma[-zero_index, -zero_index]

  # compute components of density in logs, then sum and exp()
  ld1 <- -(1/2)*log(det(2*pi*sigma_1)) # normalizing scalar
  ld2 <- as.matrix(- apply(log(x), 1, sum), ncol = 1) # n X 1
  mat_part <- log(x_1) - log(x_0) - mu_mat_1 # n x k-1
  ld3 <- diag(-(1/2)*mat_part%*%solve(sigma_1)%*%t(mat_part))
  out <- as.numeric(exp(ld1 + ld2 + ld3))

  # if not on the unit simplex, set to zero
  out[abs(apply(x, 1, sum) - 1) > tol] <- 0
  out
}

get_zero_index_and_check_errors <- function(mu, sigma){
  zero_index <- which(apply(sigma, 1, FUN = function(v){sum(v == 0)}) == ncol(sigma))
  if(length(zero_index) != 1){stop("One row/column of sigma should be all zeros.")}
  if(mu[zero_index] != 0){stop("mu should zero in the row/column where sigma is all zeros.")}
  zero_index
}


rlogisticnormal_canonical <- function(n, mu, sigma){
  #draw from multivariate normal
  draws <- cbind(mvtnorm::rmvnorm(n, mean = mu, sigma = sigma), 0)
  # logit transformation
  ed <- exp(draws)
  ed/apply(ed, 1, sum)
}

# vectorized version
dlogisticnormal_canonical <- function(x, mu, sigma){
  if(!is.matrix(x)){
    if (is.data.frame(x)){
      x <- as.matrix(x)
    }else{x <- t(x)}
  }
  mu_mat <- matrix(mu, nrow = nrow(x), ncol = length(mu), byrow = T)
  stopifnot(ncol(mu) == ncol(sigma))
  stopifnot(ncol(sigma) == nrow(sigma))

  ld1 <- -(1/2)*log(det(2*pi*sigma)) # normalizing scalar
  ld2 <- as.matrix(- apply(log(x), 1, sum), ncol = 1) # n X 1
  mat_part <- log(x[,1:(ncol(x) - 1)]/x[,ncol(x)]) - mu_mat # n x k-1
  ld3 <- diag(-(1/2)*mat_part%*%solve(sigma)%*%t(mat_part))
  as.numeric(exp(ld1 + ld2 + ld3))
}

# non-vectorized version
dlogisticnormal1 <- function(x, mu, sigma){
  stopifnot(length(x) - 1 == length(mu))
  stopifnot(length(mu) == nrow(sigma))
  stopifnot(nrow(sigma) == ncol(sigma))

  part1 <- 1/sqrt(det(2*pi*sigma))
  part2 <- 1/prod(x)
  part3a <- log(x[-length(x)]/x[length(x)]) - mu
  part3 <- exp(-(1/2)*t(part3a)%*%solve(sigma)%*%part3a)
  part1*part2*part3
}


# do they return the same thing? yes they do. so my log rules are correct.
test <- F
if(test){
  x <- c(.4, .4, .3, .2)
  mu <- c(1, .3, -1)
  sigma <- diag(3)*.02
  dlogisticnormal1(x, mu, sigma)
  dlogisticnormal(x, mu, sigma)


  mu <- c(1, .3, -1)
  sigmah <- diag(3)*.02
  samps <- rlogisticnormal(n = 10, mean = mu, sigma = sigmah)
  dlogisticnormal(samps, mu, sigmah)

}

