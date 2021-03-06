% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_normal.R
\name{logisticnormal}
\alias{logisticnormal}
\alias{rlogisticnormal}
\alias{dlogisticnormal}
\title{Functions for the logistic normal distribution}
\usage{
rlogisticnormal(n, mu, sigma)

dlogisticnormal(x, mu, sigma, tol = 1e-10)
}
\arguments{
\item{n}{Number of random vectors to generate}

\item{mu}{Vector of expectations (pre-logistic transformation).
One value must be zero, and it must correspond to the row/column of
\code{sigma} that is also
all zeros.}

\item{sigma}{Matrix of covariances (pre-logistic transformation).
There must be a \code{j} such that \code{sigma[j,]} and \code{sigma[,j]}
are all zeros, and \code{mu[j]} must also be zero.}

\item{x}{A vector containing one random deviate (i.e. one )}

\item{tol}{\code{dlogisticnormal} returns zero if the sum of
\code{x} is further than \code{tol} away from 1.}
}
\description{
Functions to compute the density of or generate random deviates from the
logistic normal distribution.
}
\details{
The logistic normal distribution arises when we draw from a multivariate
normal distribution and then apply the logistic transformation to the resulting
draws. It has been described as the Gaussian on the unit simplex. It allows
for more flexible forms of dependence among the components
than the better-known Dirichlet distribution.

In this parameterization, exactly one element of the random vector must be fixed
at zero, i.e. mean zero with zero variance and zero covariance with other
elements of the vector. Other parameterizations leave this element out of the
multivariate draw. The reason for including it here is that we can then reshuffle
the parameters as we cycle through e.g. pairs of candidates.
}
\examples{
mu <- c(.1, -.1, 0)
sigma <- rbind(c(.25, .05, 0),
               c(.05, .15, 0),
               c(0, 0, 0))
rlogisticnormal(10, mu = mu, sigma = sigma)
dlogisticnormal(rbind(c(.4, .35, .25), c(.3, .5, .2)), mu = mu, sigma = sigma)


}
