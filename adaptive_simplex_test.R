#### README
## 20200804-20200805 I quickly coded up the direct Monte Carlo for positional methods and Kemeny. For positional methods I tried a two-step Monte-Carlo where I separately sample second preference shares (from Beta distribution) and then first preference shares (from Dirichlet on the 3-simplex). See alternative_dirichlet_draws. I checked that this seems to give the same answers as the direct Dirichlet sampling. I had an idea to compute pivot probabilities more efficiently by drawing the second preference shares and then integrating first-preference shares along the appropriate loci. This takes quite a bit longer than the simulations and does not give the same answer. Possible explanations:
## maybe it doesn't actually make sense at a statistical level. I thought that if you could generate the same sample by first drawing second pref shares and then drawing first pref shares, then you could get the same pivot probs by first drawing second pref shares, then computing pivot probs (integrating over first pref shares), and then averaging. It still seems to make sense. But maybe it doesn't. Could spend time formalizing a bit to see if that reveals anything.
## maybe I didn't implement it correctly.
  ## simple things: are the lines I am integrating over the correct lines? would be easy to mess up the order.
  ## more complex: am I using the code correctly? could code up my own integration. there is a question of whether we are operating on the canonical simplex (a plane in 3d, where things are symmetric) or on the unit simplex (a plane in 2d, non-symmetric). does it matter?

f1 <- function(x, alpha){
  gtools::ddirichlet(as.numeric(x), alpha)
}
S1 <- cbind(c(.5, .5, 0), c(1/3, 1/3, 1/3))
S1_total <- cbind(c(1,0,0), c(0,1,0), c(0,0,1))

f2 <- function(x, alpha){
  x1 <- as.numeric(x)
  gtools::ddirichlet(c(x1, 1 - sum(x1)), alpha)
}
S2 <- cbind(c(.5, .5), c(1/3, 1/3))
S2_total <- cbind(c(1,0), c(0,1), c(0,0))

alpha <- c(10, 9, 3)
v1 <- SimplicialCubature::adaptIntegrateSimplex(f1, S1, alpha = alpha)
total_v1 <- SimplicialCubature::adaptIntegrateSimplex(f1, S1_total, alpha = alpha)
v2 <- SimplicialCubature::adaptIntegrateSimplex(f2, S2, alpha = alpha)
total_v2 <- SimplicialCubature::adaptIntegrateSimplex(f2, S2_total, alpha = alpha)

# they give different answers but when normalized they are the same:
v1$integral/total_v1$integral # total_v1$integral sums to sqrt(3)
v2$integral/total_v2$integral  # total_v2$integral sums to 1



## and here is the code for testing what I did:



# get the brexit preferences
v_vec <- votevizr::extract_vector_of_vote_shares_and_candidate_names_from_result(votevizr::brexit_prefs, split = " > ")$vector_of_vote_shares[1:6]

## now get Borda piv probs two ways
## analytical version, quasi-MC
pps <- piv_probs_positional_2_step(v_vec, s = 20, positional_s = .5, N = 1000)
apply(pps, 1, mean)/sqrt(3)

# compare to simulation version
sims <- draw_dirichlet_sims(100000, v_vec, 20)
positional_piv_probs_simulation(sims, s = .5)
# much faster, and they disagree.

# some similarity (though reversed AC and BC -- check if that's consistent
v_vec <- c(7, 25, 12, 16, 11, 15)
v_vec <- v_vec/sum(v_vec)

pps <- piv_probs_positional_2_step(v_vec, s = 20, positional_s = .5, N = 1000)
apply(pps, 1, mean)/sqrt(3)

# compare to simulation version
sims <- draw_dirichlet_sims(100000, v_vec, 20)
positional_piv_probs_simulation(sims, s = .5)
# these are similar without the reversal!

## OK checking plurality
pps <- piv_probs_positional_2_step(v_vec, s = 20, positional_s = 0, N = 1)
pps_sim <- positional_piv_probs_simulation(sims, s = 0)
pps/sqrt(3)
pps_sim
(pps/sqrt(3))/unlist(pps_sim) # not quite linear
plot(pps_sim, pps/sqrt(3))

# So I can come back to this and try to understand what's happening.
# The numerical approach is taking a long time. But if it was getting precise answers for small piv probs, that would be worth it.
# What I really wish existed was a function where you pass a function (Dirichlet distribution) and some inequalities and it integrated the resulting region.
# SimplicialCubature is getting close, but I find that I need to know a lot about the geometry to get an answer.

## SO leave this here for now.








p_vec <- c(v_vec[1]/sum(v_vec[1:2]), v_vec[3]/sum(v_vec[3:4]), v_vec[5]/sum(v_vec[5:6]))



# adaptIntegrate example

ab_piv_prob <- function(x, positional_s = .5, v_vec, s, width = .1){

  score_a <- x[1] + x[2] + positional_s*(x[3] + x[5])
  score_b <- x[3] + x[4] + positional_s*(x[1] + x[6])
  score_c <- x[5] + x[6] + positional_s*(x[2] + x[4])

  if(score_a > score_b & score_a - score_b < width & score_b > score_c){
    gtools::ddirichlet(x, v_vec*s)
  }else{
    0
  }

}

library(SimplicialCubature)
S <- CanonicalSimplex(5)
out <- adaptIntegrateSimplex(f = ab_piv_prob, S = S, fDim = 6L, positional_s = .5, v_vec = c(2,1,1,1,1,2)/7, s = 30, width = .1)


# OK simpler example. let's get the probability of a candidate winning.

a_win_prob <- function(x, v_vec, s){
  if(x[1] > x[2] & x[1] > x[3]){
    gtools::ddirichlet(x, v_vec*s)
  }else{
    0
  }
}

S <- CanonicalSimplex(2)
out <- adaptIntegrateSimplex(f = a_win_prob, S = S, fDim = 1L, v_vec = c(5,4,2)/11, s = 30)


### try to get his examples to work, then adapt

# compute the area of a unit simplex
f4 <- function(x){1}
# 2-dim integral, exact answer area of unit simplex = sqrt(3)/2 = 0.8660254...
adaptIntegrateSimplex(f4, UnitSimplexV(3))
# AE: a little confused why sqrt(3)/2 is the answer.

# Can I adapt f to integrate over a subset of this?
# how often does 1 win?
# OK minimal adaptation of this
f4a <- function(x) {if(x[1] > x[2] & x[1] > x[3]){1}else{0}}
adaptIntegrateSimplex(f4a, UnitSimplexV(3))
# yes, but it takes "too many function evaluations"

# can I adapt the simplex instead?
adaptIntegrateSimplex(f4a, cbind(c(1,0,0), c(1/3, 2/3, 0), c(1/3, 0, 2/3))) # same result, same error message.

# let's integrate dirichlet distribution over simplex?
f_dir <- function(x, v_vec = c(2, 1, 2)/5, s = 20){
  gtools::ddirichlet(as.numeric(x), v_vec*s)
}

dir_out <- adaptIntegrateSimplex(f_dir, UnitSimplexV(3))
# yields sqrt 3. why?

# can I integrate within a smaller area to get the probability of a tie, or a win?
# two options -- change the function or change the simplex.
f_dir_2 <- function(x, v_vec = c(2, 1, 2)/5, s = 20){
  if(x[1] > x[2] & x[1] > x[3]){
    gtools::ddirichlet(as.numeric(x), v_vec*s)
  }else{
    0
  }
}

adaptIntegrateSimplex(f_dir_2, UnitSimplexV(3))
# "too many function evaluations" but I get an answer -- .777.

restricted_simplex <- cbind(c(1,0,0), c(1/3, 2/3, 0), c(1/3, 0, 2/3))
adaptIntegrateSimplex(f_dir_2, restricted_simplex) # an answer -- .8288



# by the way what is the answer?
sims <- gtools::rdirichlet(1000000, alpha = 20*c(2, 1, 2)/5)
mean(sims[,1] > sims[,2] & sims[,1] > sims[,3])

# that second answer is very close if I normalize by sqrt(3).


# can I restrict more?
more_restricted <- cbind(c(1,0,0), c(.5, .5, 0), c(1/3, 1/3, 1/3), c(.5, 0, .5))
adaptIntegrateSimplex(f_dir_2, more_restricted)
# no that yields numerical zero.
maj_mat <- cbind(c(1,0,0), c(.5, .5,0), c(.5, 0, .5))
plur_mat <- cbind(c(.5, .5,0), c(.5, 0, .5), c(1/3, 1/3, 1/3))
win_array <- array(c(maj_mat, plur_mat), dim = c(3,3,2))
win_out <- adaptIntegrateSimplex(f_dir, S = win_array)
# yes that does work.
# so I am combining two triangles and I get the right answer.

# so how about integrating along a line?


# here's another example
f4 <- function(x) { 1 }
# this seems like a weird way to enter data, but it just makes a diamond around the origin.
S4 <- array( c( 1,0, 0,1, 0,1, -1,0, -1,0, 0,-1, 0,-1, 1,0) , dim=c(2,2,4) )
adaptIntegrateSimplex( f4, S4 ) # so this is the length of that -- sqrt(2)*4, and that's what we get.

# so I can do Dirichlet tie probs like this.
dir_f <- function(x, alpha){gtools::ddirichlet(c(as.numeric(x), 1 - sum(x)), alpha)}
ab_tie_for_first_2d <- cbind(c(1/3, 1/3), c(1/2, 1/2))
adaptIntegrateSimplex(dir_f, ab_tie_for_first_2d, alpha = c(20, 18, 6))

# OK that looks pretty good. how can I extend that?
# what about plurality in four dimensions?

# this gives a numerical zero
ab_tie_for_first_3d <- cbind(c(1/2, 1/2, 0), c(1/3, 1/3, 0), c(1/4, 1/4, 1/4), c(1/3, 1/3, 1/3))
adaptIntegrateSimplex(dir_f, ab_tie_for_first_3d, alpha = c(20, 18, 6, 4))

# what if I do like I did before?
mat_1 <- cbind(c(.5,.5,0,0), c(1/3,1/3,0, 1/3), c(1/3, 1/3, 1/3, 0))
mat_2 <- cbind(c(1/3,1/3,0, 1/3), c(1/3, 1/3, 1/3, 0), rep(1/4, 4))
ab_tie_array <- array(c(mat_1, mat_2), dim = c(4,3,2))
v_vec <- c(5,4,3,2)/14
tie_out_ab <- adaptIntegrateSimplex(f_dir, S = ab_tie_array, v_vec = v_vec, s = 20)
# wow that does something. is it correct?
tie_out_ac <- adaptIntegrateSimplex(f_dir, S = ab_tie_array, v_vec = v_vec[c(1,3,2,4)], s = 20)
tie_out_ad <- adaptIntegrateSimplex(f_dir, S = ab_tie_array, v_vec = v_vec[c(1,4,2,3)], s = 20)
# looks sensible.

# If this is working (and I can check more carefully), then the question is whether I can extend this to e.g. Borda.
# so let's check it
source("~/Dropbox/research/partisanship_corruption_survey_experiment/analysis/data_prep/pivotality_functions.R")
vv <- v_vec
names(vv) <- LETTERS[1:4]
pp <- pivotal.probabilities.analytical(vv*20)

plot(pp[c("p.B.A", "p.C.A", "p.D.A")], c(tie_out_ab$integral, tie_out_ac$integral, tie_out_ad$integral), pch = 19)

# so that looks good.

# now, can I make any progress on the Borda count case?
# I can't figure it out.

# what about this idea of drawing the second preference ratios first, from a beta distribution, and then the first preferences. Do we get the same thing?

alternative_dirichlet_draws <- function(N, v_vec, s){
  alpha_vec <- v_vec*s
  fp_shares <- gtools::rdirichlet(N, alpha = c(sum(alpha_vec[1:2]), sum(alpha_vec[3:4]), sum(alpha_vec[5:6])))
  ps <- data.frame(
    p_ab = rbeta(N, alpha_vec[1], alpha_vec[2]),
    p_ba = rbeta(N, alpha_vec[3], alpha_vec[4]),
    p_ca = rbeta(N, alpha_vec[5], alpha_vec[6])
  )
  sims <- data.frame(
    ab = fp_shares[,1]*ps$p_ab,
    ac = fp_shares[,1]*(1 - ps$p_ab),
    ba = fp_shares[,2]*ps$p_ba,
    bc = fp_shares[,2]*(1 - ps$p_ba),
    ca = fp_shares[,3]*ps$p_ca,
    cb = fp_shares[,3]*(1 - ps$p_ca)
    )
  list(fp_shares = fp_shares,
       ps = ps,
       sims = sims)
}

v_vec <- c(25, 12, 25, 27, 12, 30)
v_vec <- v_vec/sum(v_vec)
N <- 100000
s <- 30
sims_1 <- data.frame(gtools::rdirichlet(N, alpha = v_vec*s))
colnames(sims_1) <- c("ab", "ac", "ba", "bc", "ca", "cb")
sims_2_list <- alternative_dirichlet_draws(N, v_vec, s)
sims_2 <- sims_2_list$sims

# OK so I think this does work. I did some checks on specific intervals and they agreed.

# this implies that we can draw the p's, then numerical integration.
# this helps because once we have the p's, it's just integrating along a line. so to get this to work I would need to get the mapping from p's to endpoints of the positional tie lines.
# can I get those endpoints from votevizr?



