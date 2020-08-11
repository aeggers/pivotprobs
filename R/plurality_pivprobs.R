# analytical plurality piv probs
# Eggers & Vivyan 2020 approximation

# TODO: clean this up and compare to the grid version


# simulation based
pivotal.probabilities = function(samp, tol = .01){
  rank.mat = t(apply(-samp, 1, rank)) # why is there a transpose?
  out <- list()
  for(i in 1:ncol(samp)){
    for(j in 1:ncol(samp)){
      if(i >= j){next} # this will be a lot faster.
      out[[paste0(colnames(samp)[i], colnames(samp)[j])]] = sum((rank.mat[,i] + rank.mat[,j]) == 3 & abs(samp[,i] - samp[,j]) < tol/2)/(tol*nrow(samp))
    }
  }
  out
}

win.shares = function(samp){
  max.of.row = apply(samp, 1, max, na.rm = T)
  out = c()
  for(j in 1:ncol(samp)){
    out = c(out, sum(samp[,j] == max.of.row))
  }
  names(out) = colnames(samp)
  out/sum(out)
}


### new analytical (naive) version

library(gtools)
# require(ExtDist)

# first note that pBeta_ab gives the same thing as pbeta when we rescale the variable.
# y = .27
# alpha.vec = c(40, 30, 20, 10, 5)
# i = 3
# pBeta_ab(y, alpha.vec[i], sum(alpha.vec[-c(1,2,i)]), a = 0, b = 1 - 2*y)
# pbeta(y/(1 - 2*y), alpha.vec[i], sum(alpha.vec[-c(1,2,i)]))

napotffay = function(alpha.vec, y){
  # naive.analytical.probability.of.tie.for.first.at.y
  # probability that parties 1 and 2 tie at y can be calculated without specifying the values of the other parties, by the aggregation property.
  prob.of.tie.at.y = ddirichlet(x = c(y,y,1-2*y), alpha = c(alpha.vec[1], alpha.vec[2], sum(alpha.vec) - alpha.vec[1] - alpha.vec[2]))
  # now, what is the probability that _no other party_ would have gotten above y?
  # we do it separately for each party -- a naive approach not recognizing how the vote shares for the other parties are related.
  # obvious they are correlated because they must sum to 1-2*y.
  probs = c()
  for(i in 3:length(alpha.vec)){
    this.prob = pbeta(y/(1-2*y), alpha.vec[i], sum(alpha.vec[-c(1,2,i)])) # this is the part I'm least sure of
    # 		this.prob = pBeta_ab(y, alpha.vec[i], sum(alpha.vec[-c(1,2,i)]), a = 0, b = 1 - 2*y) # this is the part I'm least sure of
    probs = c(probs, this.prob)
  }
  list(prob.of.tie.at.y = prob.of.tie.at.y, probs.of.being.below.y = probs, probability = prod(c(prob.of.tie.at.y, probs)))
}

# now we aggregate this up (integrate over values of y)
probability.of.tie.for.first = function(alpha.vec, increments = 50){
  ys = seq(1/length(alpha.vec), .5, length = increments) # the least a pair of parties can get and be tied for first is 1/k each. the most they can get is .5.  We are going to calculate the probability of a tie for first at each value in that sequence and then integrate.
  probs = c() # this holds the probability of the two parties being tied for first at each y
  tie.vec = c()  # this holds the probability of the two parties being tied at a given y
  prob.mat = matrix(NA, nrow = 0, ncol = length(alpha.vec)-2) # this holds the probabilities of being below y for each of the parties other than the ones who are tied.
  increment = ys[2] - ys[1]
  for(i in 1:(length(ys) - 1)){
    this = napotffay(alpha.vec, (ys[i] + ys[i+1])/2)
    probs = c(probs, increment*this$probability)
    tie.vec = increment*c(tie.vec, this$prob.of.tie.at.y)  # this seems wrong, though it's not used.
    prob.mat = rbind(prob.mat, this$probs.of.being.below.y)
  }
  list(estimate = sum(probs), # apply(cbind(tie.vec, prob.mat), 1, prod)),  # for each y, the probability of
       ys = ys, tie.vec = tie.vec, prob.mat = prob.mat) # I include these for diagnostic purposes
}

## and here's a function that yields pivotal probabilities in a vector, given an alpha input.
pivotal.probabilities.analytical = function(alpha.vec, increments = 50, n = 50000){
  out <- list()
  for(i in 1:length(alpha.vec)){
    for(j in 1:length(alpha.vec)){
      if(i >= j){next}
      out[[paste0(names(alpha.vec)[i], names(alpha.vec)[j])]] = probability.of.tie.for.first(alpha.vec = c(alpha.vec[i], alpha.vec[j], alpha.vec[-c(i,j)]), increments = increments)$estimate/n
    }
  }
  out
}



# myatt-fisher approach

incomplete.beta = function(x,a,b){
  pbeta(x,a,b)*beta(a,b)
}

posterior.pivotal.probability = function(pie.vec, s, which.ones = c(1,1,0), prior = 1){
  pie = pie.vec[which(which.ones == 0)]
  inc.beta.part = incomplete.beta(1/3, pie*s + prior, (1 - pie)*s + prior)
  2^(pie*s)*inc.beta.part/2^(s + prior)  # the 2016 wp version. to ask Myatt: is the +1 in the denominator the prior? is 2^(pie*s) actually 2^(pie*s + 1 - prior)?
  # (2^-((1 - pie)*s + 1))*inc.beta.part  # their old specification -- they are the same.
}

# check against my thing
test = F
if(test){
  v.vec = c(.41, .39, .2)
  s = 85
  # MF version
  prior = 1
  posterior.pivotal.probability(pie.vec = v.vec, s = s, prior = prior)/posterior.pivotal.probability(pie.vec = v.vec, s = s, which.ones = c(0, 1, 1), prior = prior)
  # my analytical version
  alpha.vec = v.vec*s + prior
  increments = 200
  probability.of.tie.for.first(v.vec* s, increments = increments)$estimate/probability.of.tie.for.first(alpha.vec[c(2,3,1)], increments = increments)$estimate

  # hmm that's not good.
  # issue 1: the increment matters a lot when the probability of a tie is low (i.e. for smaller parties).
  ## OK I seem to have gotten rid of that variation by having y calculated at the midpoint of the bins. there is a bit of variation but it is smaller.
  # issue 2 is the absolute size of these things. but I think this is arbitrary -- we're calculating this along a surface, and effectively assuming a width of 1 in one direction and increment in the other. So the 1 can be replaced by .00001 or whatever. And even the Myatt-Fisher thing is "propto" not exactly equal.
  # issue 3 is that the ratios are not the same.
  # in simulations the ratio bounces around, but my best answer is that the simulations are a little closer to the Myatt-Fisher version.
  # the fact that the relationship is linear but not exactly linear is interesting. I don't know what this could be.
  # one remaining thing is to check the prior again. I set the prior to 1. I do get closer when I have a prior of 1, but I don't know if this is because the MF equation assumes a prior of 1 or what.

  # one thing I could do as a test is:

  m = 500
  v.vecs = rdirichlet(m, alpha = c(10,9,6))
  s = 30
  sto = matrix(NA, nrow = 0, ncol = 3)
  colnames(sto) = c("mf", "me", "sim")
  n = 50000
  sim = F

  for(i in 1:m){
    if(i %% 10 == 0){cat(i, " ", sep = "")}
    v.vec = v.vecs[i,]
    v.vec = sort(v.vec, decreasing = T)
    # MF
    mf.ratio = posterior.pivotal.probability(pie.vec = v.vec, s = s, prior = 1)/posterior.pivotal.probability(pie.vec = v.vec, s = s, which.ones = c(0, 1, 1), prior = 1)
    # ME
    alpha.vec = v.vec*s + 1
    my.ratio = probability.of.tie.for.first(v.vec* s)$estimate/probability.of.tie.for.first(alpha.vec[c(2,3,1)], increments = 100)$estimate
    sim.ratio = NA
    if(sim){
      # SIM
      dir.sample = rdirichlet(n, alpha = alpha.vec)
      colnames(dir.sample) = c("A", "B", "C")
      pp.vec = pivotal.probabilities(dir.sample, out.as.vector = T)
      sim.ratio = pp.vec[1]/pp.vec[3]

      if(is.infinite(sim.ratio)){next}
    }
    sto = rbind(sto, c(mf.ratio, my.ratio, sim.ratio))
    if(i > 5){
      plot(sto[,1], sto[,2], pch = 19, col = rgb(.3, .3, .3, alpha = .5), cex = .5, xlab = "MF", ylab = "ME & SIM")
      abline(a = 0, b =1, lty = 2)
      abline(lm(sto[,2] ~ sto[,1]))
      if(sim){
        points(sto[,1], sto[,3], pch = 19, col = rgb(.8, .1, .1, alpha = .5), cex = .5)
        abline(lm(sto[,3] ~ sto[,1]), col = "red")
      }
    }
  }
  # sto = sto[!is.infinite(sto[,3]), ]


  # very systematic discrepancy -- nearly perfectly linear.
  # it barely depends on the increment.
  # it does depend on the prior assumed in the Myatt-Fisher thing.
  # I wonder if the error is from not normalizing. right now we are integrating rectangles. maybe this should be a big triangle. i.e. we should be weighting the stuff -- no that seems arbitrary.

  v.vecs[sto[,1] > 400,]

}
# what we can say now is that the correlation is almost 1. but note that the ratio being different means that the predictions/implications could be different.


