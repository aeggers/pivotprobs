#' Make the P matrix (election probability matrix)
#'
#' Takes the output of \code{election_event_probabilities()} and returns a
#' matrix (the P matrix) indicating the probability of each outcome (rows)
#' as a function of the additional ballot (columns). If \code{u} is a vector
#' of utilities (one per outcome), then \code{t(P)%*%u} gives the expected utility
#' of each ballot.
#'
#' @param election_event_probs The output of \code{election_event_probabilities()};
#' a list of election events, each of which is a list with a probability
#' (\code{integral}) and an event-specific P matrix (\code{P}).
#'
#' @examples
#' plurality_election(k = 3) %>%
#'    election_event_probabilities(method = "sc", alpha = c(10, 8, 5)) %>%
#'    combine_P_matrices -> P
#' t(P) %*%  c(2, 5, 9)

#' @export
combine_P_matrices <- function(election_event_probs){
  eep <- election_event_probs
  # could also do a dplyr way -- this way is loops
  P <- eep[[names(eep)[1]]]$P*eep[[names(eep)[1]]]$integral
  for(j in 2:length(names(eep))){
    the_integral <- eep[[names(eep)[j]]]$integral
    if(is.null(the_integral)){next}
    if(is.na(the_integral)){next}
    P <- P + eep[[names(eep)[j]]]$P*the_integral
  }
  P
}
