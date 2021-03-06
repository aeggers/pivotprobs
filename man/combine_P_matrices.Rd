% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_P_matrices.R
\name{combine_P_matrices}
\alias{combine_P_matrices}
\title{Make the P matrix (election probability matrix)}
\usage{
combine_P_matrices(election_event_probs)
}
\arguments{
\item{election_event_probs}{The output of \code{election_event_probs()} or
a standalone Monte Carlo function;
a list of election events, each of which is a list with a probability
(\code{integral}) and an event-specific P matrix (\code{P}).}
}
\description{

}
\examples{
plurality_election(k = 3) \%>\%
   election_event_probs(method = "sc", alpha = c(10, 8, 5)) \%>\%
   combine_P_matrices -> P
t(P) \%*\%  c(2, 5, 9)
}
