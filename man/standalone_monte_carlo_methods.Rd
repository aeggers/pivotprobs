% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standalone_monte_carlo_methods.R
\name{standalone_monte_carlo_methods}
\alias{standalone_monte_carlo_methods}
\alias{plurality_event_probs_from_sims}
\alias{positional_event_probs_from_sims}
\alias{irv_event_probs_from_sims}
\alias{condorcet_event_probs_from_sims}
\title{Faster Monte Carlo estimation of pivot probabilities}
\usage{
plurality_event_probs_from_sims(
  sims = NULL,
  n = 1000,
  window = 0.01,
  cand_names = NULL,
  sep = "_",
  method = "density",
  merge = F,
  drop = F,
  bw_divisor = 1,
  skip_non_pivot_events = T,
  s = NULL,
  raw = F
)

positional_event_probs_from_sims(
  sims,
  window = 0.01,
  n = 1000,
  s = 0.5,
  cand_names = NULL,
  sep = "_",
  method = "density",
  merge = F,
  drop = F,
  increments = 10,
  bw_divisor = 1,
  skip_non_pivot_events = T,
  raw = F
)

irv_event_probs_from_sims(
  sims,
  window = 0.01,
  n = 1000,
  s = 0,
  cand_names = NULL,
  sep = "_",
  method = "density",
  merge = F,
  bw_divisor = 1,
  skip_non_pivot_events = T
)

condorcet_event_probs_from_sims(
  sims,
  n = 1000,
  window = 0.01,
  cand_names = NULL,
  sep = "_",
  kemeny = T,
  method = "density",
  merge = F,
  bw_divisor = 1,
  s = NULL,
  skip_non_pivot_events = T
)
}
\arguments{
\item{sims}{A matrix of simulated election results, with one column per
ballot type. Must be 6 columns for the ordinal methods (IRV, Kemeny-Young, positional).}

\item{n}{Size of electorate}

\item{window}{Window within which two candidates are considered to be tied
when \code{method="rectangular"}. Wider window means lower variance but
more bias.}

\item{cand_names}{Names of the candidates.}

\item{sep}{Separation between candidate names.}

\item{merge}{Merge adjacent pivot events?}

\item{drop}{Drop a dimension?}

\item{skip_non_pivot_events}{Skip non pivot events?}

\item{s}{Score allocated to a second-place ranking in positional
and IRV elections.}

\item{raw}{Return counts rather than integrals/proportions?}

\item{kemeny}{For Condorcet method, compute Kemeny-Young pivot event probabilities? If \code{F}, only handles event in which the Condorcet winner
is decided.}
}
\description{
These methods take simulated election results and return a list of
pivot events, each with a pivot probability (\code{integral}) and
P matrix (\code{P}). They produce identical results to
\code{election_event_probs()} (or nearly so) but they are faster because
they are written in a less general way.
}
\details{
For plurality elections, any number of candidates can be specified.
The other voting systems can only handle three candidates.

Results have been validated against \code{election_event_probs()}.
}
\examples{
sims <- gtools::rdirichlet(100000, alpha = c(9,7,3,4,4,6))
plurality_pivot_probs_from_sims(sims)
positional_pivot_probs_from_sims(sims, s=.5)
irv_pivot_probs_from_sims(sims)
condorcet_pivot_probs_from_sims(sims)


}
