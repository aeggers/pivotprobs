
eggers_vivyan_plurality_pivot_probs <- function(alpha = NULL, mu = NULL, precision = NULL, increments = 50, cand_names = NULL, sep = ""){
  if(is.null(alpha)){
    if(is.null(mu) | is.null(precision)){
      stop("You must pass either 'alpha' OR 'mu' and 'precision'.")
    }
    alpha <- mu*precision
  }
  if(is.null(cand_names)){cand_names <- letters[1:3]}
  out <- list()
  for(i in 1:(length(alpha) - 1)){
    for(j in (i + 1):length(alpha)){
      out[[paste0(cand_names[i], sep, cand_names[j])]] <- probability_of_tie_for_first(alpha = c(alpha[c(i, j)], alpha[-c(i,j)]), increments = increments)$estimate # I once divided by n here -- but this can be done later.
    }
  }
  out
}


napotffay <- function(alpha, y){
  # naive.analytical.probability.of.tie.for.first.at.y
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

# now we aggregate this up (integrate over values of y)
probability.of.tie.for.first <- function(alpha, increments = 50){
  ys = seq(1/length(alpha), .5, length = increments) # the least a pair of parties can get and be tied for first is 1/k each. the most they can get is .5.  We are going to calculate the probability of a tie for first at each value in that sequence and then integrate.
  probs = c() # this holds the probability of the two parties being tied for first at each y
  tie.vec = c()  # this holds the probability of the two parties being tied at a given y
  prob.mat = matrix(NA, nrow = 0, ncol = length(alpha)-2) # this holds the probabilities of being below y for each of the parties other than the ones who are tied.
  increment = (ys[2] - ys[1])*sqrt(6) # if d is the step length in the ys vector (i.e. ys[2] - ys[1]), then the segment length on the unit simplex is sqrt(d^2 + d^2 + (-2d)^2) = sqrt(6) d
  for(i in 1:(length(ys) - 1)){
    this = napotffay(alpha, (ys[i] + ys[i+1])/2)
    probs = c(probs, increment*this$probability)
    tie.vec = increment*c(tie.vec, this$prob.of.tie.at.y)  # this seems wrong, though it's not used.
    prob.mat = rbind(prob.mat, this$probs.of.being.below.y)
  }
  list(estimate = sum(probs),
       ys = ys, tie.vec = tie.vec, prob.mat = prob.mat) # I include these for diagnostic purposes
}


probability_of_tie_for_first <- function(alpha, increments = 50){
  boundary_points <- seq(from = 1/length(alpha), to = .5, length = increments+1) # the least a pair of parties can get and be tied for first is 1/k each. the most they can get is .5.
  midpoints <- apply(cbind(boundary_points[1:increments], boundary_points[2:(increments+1)]), 1, mean) # the midpoints of those grid segments
  # storage for a for-loop (could use map etc but trying to keep base R)
  tie_vec <- probs <- rep(NA, length(midpoints)) # probs holds the probability of the two parties being tied for first at each point; tie_vec holds the probability of the two parties being tied at a given point
  prob_mat <- matrix(NA, nrow = length(midpoints), ncol = length(alpha)-2) # this holds the probabilities of being below y for each of the parties other than the ones who are tied.
  for(i in 1:length(midpoints)){
    this <- napotffay(alpha, midpoints[i])
    probs[i] <- this$probability
    tie_vec[i] <- this$prob.of.tie.at.y  # cryptic comment here: "this seems wrong, though it's not used".
    prob_mat[i,] <- this$probs.of.being.below.y
  }
  increment <- (midpoints[2] - midpoints[1])*sqrt(6) # if d is the step length in the ys vector (i.e. ys[2] - ys[1]), then the segment length on the unit simplex is sqrt(d^2 + d^2 + (-2d)^2) = sqrt(6) d
  list(estimate = sum(probs)*increment,
       ys = midpoints, tie_vec = tie_vec, prob_mat = prob_mat) # I include these for diagnostic purposes
}

