#' Fast computation of plurality pivot probabilities given Dirichlet beliefs
#'
#' Returns probability of a tie for first between two candidates in a
#' k-candidate plurality election given Dirichlet beliefs, not normalized for
#' electorate size. Uses the approximation introduced by Eggers & Vivyan (2020),
#' which (when there are more than three candidates) involves an
#' independence assumption on the vote shares of candidates
#' 3 to k.
#'
#' The estimates are not normalized for the size of the electorate (and thus could
#' be larger than 1). Dividing by
#' the electorate size gives (an approximation of) the true pivot probability.
#'
#' Dirichlet beliefs for a k-candidate plurality election are characterized by
#' a length-k vector \code{alpha}. It may be helpful to think of \code{alpha}
#' as the product of a vector of expected vote shares
#' \code{(mu_1, mu_2, ... mu_k)} and a
#' scalar precision parameter \code{s}. Eggers & Vivyan (2020) find that the
#' precision of UK election forecasts is characterized by \code{s=85}.
#' So, given an expected
#' result of \code{c(.4, .35, .25)} we might model the result by setting
#' \code{alpha} to \code{c(.4, .35, .25)*85}.
#'
#' For k=3 candidates the Eggers-Vivyan method produces exact pivot probabilities
#' (as \code{increments} goes to infinity, and after dividing by electorate size).
#'
#' For k>3 candidates the probability of unlikely ties for first is overstated due to
#' the independence assumption: the probability of candidates 1 and 2 tying for first
#' at a given vote share x is approximated by the product of
#' \itemize{
#' \item the probability of candidates 1 and 2 receiving vote share x
#' \item the probability of candidate 3 receiving below x (given 1 and 2 each get x)
#' \item the probability of candidate 4 receiving below x (given 1 and 2 each get x)
#' \item etc up to k}
#'
#' @examples
#' eggers_vivyan_probability_of_tie_for_first(c(10, 7, 5)) # non-normalized pivot prob for 3 candidates
#' eggers_vivyan_probability_of_tie_for_first(c(10, 7, 5))/100000 # normalized for electorate size
#' eggers_vivyan_probability_of_tie_for_first(c(10, 7, 5, 3)) # 4 candidates
#' eggers_vivyan_probability_of_tie_for_first(c(10, 7, 5, 3), increments = 100) # more precise
#'
#' @param alpha Length-k parameter vector for Dirichlet belief distribution.
#' See Details.
#' @param increments Number of points at which to compute the probability
#' of candidates 1 and 2 tying for first. More increments means
#' a more precise estimate.
#'
#' @export
eggers_vivyan_probability_of_tie_for_first <- function(alpha, increments = 50){

  boundary_points <- seq(from = 1/length(alpha), to = .5, length = increments+1) # the least a pair of parties can get and be tied for first is 1/k each. the most they can get is .5.

  midpoints <- apply(cbind(boundary_points[1:increments], boundary_points[2:(increments+1)]), 1, mean) # the midpoints of those grid segments

  # storage for a for-loop (could use map etc but trying to keep base R)
  tie_vec <- probs <- rep(NA, length(midpoints)) # probs holds the probability of the two parties being tied for first at each point; tie_vec holds the probability of the two parties being tied at a given point
  prob_mat <- matrix(NA, nrow = length(midpoints), ncol = length(alpha)-2) # this holds the probabilities of being below y for each of the parties other than the ones who are tied.

  for(i in 1:length(midpoints)){
    this <- probability_of_tie_for_first_at_y_naive(alpha, midpoints[i])
    probs[i] <- this$probability
    tie_vec[i] <- this$prob.of.tie.at.y  # cryptic comment here: "this seems wrong, though it's not used".
    prob_mat[i,] <- this$probs.of.being.below.y
  }
  increment <- (midpoints[2] - midpoints[1])*sqrt(6) # if d is the step length in the ys vector (i.e. ys[2] - ys[1]), then the segment length on the unit simplex is sqrt(d^2 + d^2 + (-2d)^2) = sqrt(6) d
  list(estimate = sum(probs)*increment,
       ys = midpoints, tie_vec = tie_vec, prob_mat = prob_mat) # I include these for diagnostic purposes
}

probability_of_tie_for_first_at_y_naive <- function(alpha, y){
  # probability that parties 1 and 2 tie at y can be calculated without specifying the values of the other parties, by the aggregation property.
  prob.of.tie.at.y <- gtools::ddirichlet(x = c(y,y,1-2*y), alpha = c(alpha[1], alpha[2], sum(alpha) - alpha[1] - alpha[2]))/sqrt(3) # correction necessary because integration on whole d-simplex yields sqrt(d)
  # now, what is the probability that _no other party_ would have gotten above y?
  # we do it separately for each party -- a naive approach not recognizing how the vote shares for the other parties are related.
  # obvious they are correlated because they must sum to 1-2*y.
  probs <- rep(NA, length(alpha) - 2)
  for(i in 3:length(alpha)){
    this.prob <- pbeta(y/(1-2*y), alpha[i], sum(alpha[-c(1,2,i)])) # this is the part I'm least sure of
    probs[i - 2] <- this.prob
  }
  list(prob.of.tie.at.y = prob.of.tie.at.y, probs.of.being.below.y = probs, probability = prod(c(prob.of.tie.at.y, probs)))
}
