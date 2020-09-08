
# Note below we have functions that handle the general case with truncated ballots. is just the simpler case.


# These functions return pivot probabilities that are not normalized for population. i.e. divide by N to get the true pivot probability.

# second-round pivot events
eggers_nowacki_second_round_pivot_probability <- function(alpha, increments = 20){
  stopifnot(length(alpha) == 6)
  midpoints.y <- get.increment.midpoints(1/4, 1/2, increments)  # the y values we will be evaluating
  int.limits <- 2*midpoints.y - .5
  int.limits[midpoints.y > 1/3] <- midpoints.y[midpoints.y > 1/3]/2
  beta.parts <- pbeta(2*int.limits, alpha[6], sum(alpha[3:4])) # vectorized: calculating the probability of \pi_{CB}/(\pi_CB + \pi_B) being below
  y.mat <- cbind(midpoints.y, 1/2 - midpoints.y, 1/2)
  dir.parts <- apply(y.mat, 1, gtools::ddirichlet, alpha = c(sum(alpha[c(1,2)]), alpha[5], sum(alpha[c(3,4,6)])))/sqrt(3) # correction on gtools::ddirichlet() necessary because integrating that function over whole d-simplex yields sqrt(d)
  # NOTE I previously was dividing by 2 in the line above. I don't remember why. Was it about stretching pbeta part?
  stopifnot(length(beta.parts) == length(dir.parts))
  width_of_channel <- sqrt(6)/4  # (times N)
  distance_between_midpoints <- sqrt(2)*(midpoints.y[2] - midpoints.y[1])
  width_of_channel*distance_between_midpoints*sum(dir.parts*beta.parts) # see appendix of AJPS (?) paper for explanation of width_of_channel term. basically, the width is the distance between the points (1/4, 1/4, 1/2) and (1/4 -1/4N, 1/4 - 1/4N, 1/2 + 1/2N).
  # the distance_between_midpoints term is just a numerical integration point -- distance between (y, 1/2 - y, 1/2) and (y + d, 1/2 - y - d, 1/2) is sqrt(2)*d.
}

# first-round pivot events
eggers_nowacki_first_round_pivot_probabilities <- function(alpha, increments = 100){
  stopifnot(length(alpha) == 6)
  fp_alpha <- c(sum(alpha[1:2]), sum(alpha[3:4]), sum(alpha[5:6]))
  midpoints.y = get.increment.midpoints(1/4, 1/3, increments)  # the y values we will be evaluating
  pr.tie.for.second.at.y <- pr.1.beats.3.given.y <- pr.2.beats.3.given.y <- rep(NA, length(midpoints.y))
  for(i in 1:length(midpoints.y)){
    y <- midpoints.y[i]
    pr.tie.for.second.at.y[i] <- gtools::ddirichlet(c(y, y, 1-2*y), alpha = fp_alpha)/sqrt(3) # normalization because integral over whole d-dimensional unit simplex yields sqrt(d)
    pr.1.beats.3.given.y[i] <- pbeta(2 - 1/(2*y), alpha[4], alpha[3])
    pr.2.beats.3.given.y[i] <- pbeta(2 - 1/(2*y), alpha[2], alpha[1])
  }
  distance_between_points <- (midpoints.y[2] - midpoints.y[1])*sqrt(6) # because distance between (y, y, 1-2y) and (y+d, y+d, 1-2(y+d)) is sqrt(6)*d
  width_of_channel <- 1/sqrt(2)  # (times N)
  list(
    "i_j|ij" = distance_between_points*width_of_channel*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*pr.2.beats.3.given.y),
    "i_j|kj" = distance_between_points*width_of_channel*sum(pr.tie.for.second.at.y*(1 - pr.1.beats.3.given.y)*pr.2.beats.3.given.y),
    "i_j|ik" = distance_between_points*width_of_channel*sum(pr.tie.for.second.at.y*pr.1.beats.3.given.y*(1 - pr.2.beats.3.given.y))
  )
}

## utility used above
get.increment.midpoints = function(start, end, increments, increments.at.end = NULL){
	inc.length = (end - start)/increments
	out = seq(start + inc.length/2, end - inc.length/2, length = increments)
	if(!is.null(increments.at.end)){
		to.add.at.end = get.increment.midpoints(mean(out[c(length(out) - 1, length(out))]), end, increments.at.end)
		out = c(out[-length(out)], to.add.at.end)
	}
	out
}



# Below we have code that generalizes to the case with truncated ballots, but it's not yet used.



##################################################################
######### key command to get all pivotal probabilities ###########
##################################################################

av.pivotal.event.probs.general = function(v.vec, s.vec, increments = 20, increments.at.end = 10, increments.for.dir.mat = 20){
  c(second.round.pivotal.event.probs.general(v.vec, s.vec, increments, increments.at.end), first.round.pivotal.event.probs.general(v.vec, s.vec, increments, increments.for.dir.mat))
}


#####################################
#### second-round pivotal events ####
#####################################


pr.second.round.pivotal.cands.1.and.2.general = function(v.vec, s.vec, increments = 20, increments.at.end = 10, alpha.min = .5, noisy = F, list.please = F){
	stopifnot(length(v.vec) == 9)
	stopifnot(length(s.vec) == 4)
	if(length(unique(s.vec)) == 1 & sum(v.vec[c(7,8,9)]) == 0){ # no truncated ballots, one precision parameter
		if(noisy){cat("Using efficient version (numerical integration along only one dimension).\n")}
		midpoints.y = get.increment.midpoints(1/4, 1/2, increments)  # the y values we will be evaluating
		y.inc = midpoints.y[2] - midpoints.y[1]
		alpha.vec = v.vec*s.vec[1]
		# more efficient approach for this case
		dir.parts = c()
		int.limits = 2*midpoints.y - .5
		int.limits[midpoints.y > 1/3] = midpoints.y[midpoints.y > 1/3]/2
		beta.parts = pbeta(int.limits/.5, alpha.vec[6], sum(alpha.vec[3:4])) # vectorized: calculating the probability of \pi_{CB}/(\pi_CB + \pi_B) being below
		y.mat = cbind(midpoints.y, 1/2 - midpoints.y, 1/2)
		dir.parts = apply(y.mat, 1, ddirichlet, alpha = c(sum(alpha.vec[c(1,2)]), alpha.vec[5], sum(alpha.vec[c(3,4,6)])))/2  # NOTE dividing by 2
		stopifnot(length(beta.parts) == length(dir.parts))
		if(list.please){return(list(
			pr.ab.tie.unconditional = y.inc*sum(dir.parts), # probability of TPP ties: A and B tying
			pr.ab.tie.unconditional.2 = dbeta(.5, sum(alpha.vec[c(1,2,5)]), sum(alpha.vec[c(3,4,6)])),
			pr.ab.tie.conditional = y.inc*(sum(dir.parts*beta.parts)), # probability of TPP ties and these are the two who advance
			pr.c.last.conditional = cbind(midpoints.y, beta.parts)  # probability at each value of y that C is last (y indicates \pi_A, and 1/2 - y is \pi_{CA})
			))}
		return(y.inc*sum(dir.parts*beta.parts))
	}else{    # truncated ballots and/or more than one precision parameter
		# we use a grid approach
		do.beta = v.vec[9]*s.vec[4] < alpha.min  # this is done if there are no or few truncated ballots in the key spot.
		if(noisy){cat("Using grid approach with ", ifelse(do.beta, "beta density for second pref calcs (no or few truncated ballots).\n", "integration along dirichlet line for second pref calcs.\n"), sep = "")}
		midpoints.ab = get.increment.midpoints(1/4, 1/2, increments)  # the A and B values we will be evaluating in the FP plot
		ab.inc = midpoints.ab[2] - midpoints.ab[1]
		fp.grid = sp.grid = matrix(0, length(midpoints.ab), length(midpoints.ab))
		fp.alpha = s.vec[1]*c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
		# now we cycle through the grid
		for(i in 1:length(midpoints.ab)){
			a.val = midpoints.ab[i]  # pi_A
			for(j in 1:length(midpoints.ab)){
				y = midpoints.ab[j]  # pi_B
				x = 1 - a.val - y # pi_C
				if(x > y | x > a.val){next}
				vec.3 = c(1-y-x, y, x)
				# cat("At piA ", round(vec.3[1], 2), ", piB ", round(vec.3[2], 2), ", piC ", round(vec.3[3], 2), ".\n", sep = "")
				fp.grid[i, j] = ab.inc^2*ddirichlet(vec.3, fp.alpha) # the probability of being at this particular combination of first-preference shares.
				if(do.beta){
					# if the proportion of truncated ballots for C (cand 3) is 0 or very low, we just get the density at the key value on the beta distribution
					# so this would always be the route taken if truncated ballots are disallowed.
					sp.grid[i,j] = (1/(2*x))*dbeta((1 - 2*y)/(2*x), v.vec[6]*s.vec[4], v.vec[5]*s.vec[4])  # dividing by 2 because of the tie approximation
				}else{
					# otherwise we need to integrate along a line through C's dirichlet distribution
					lower.limit = 0
					upper.limit = 1/2 - y
					how.many.z = round((upper.limit - lower.limit)/ab.inc)
					sum = 0
					if(how.many.z > 1){
						midpoints.z = get.increment.midpoints(lower.limit, upper.limit, how.many.z, increments.at.end = increments.at.end)  # the z values (\pi_{CB}) we will be evaluating
						basic.inc = midpoints.z[2] - midpoints.z[1]
						increments.z = c(rep(basic.inc, how.many.z - 1), rep(basic.inc/increments.at.end, increments.at.end))
						picbs = midpoints.z
						picxs = 1-2*y - 2*picbs
						picas = x + picbs - 1 + 2*y
						which.to.do = which(picxs > 0 & picxs < x & picas > 0 & picas < x)
						mat = cbind(picbs, picxs, picas)
						norm.mat = mat/apply(mat, 1, sum)
						dir.parts = apply(norm.mat[which.to.do,], 1, ddirichlet, alpha = v.vec[c(6,9,5)]*s.vec[4])
						increments.vec = increments.z[which.to.do]
						sum = sum(increments.vec*dir.parts)
					}
					sp.grid[i,j] = (1/(x^2))*sum # normalization because the density is originally on 0-x, 0-x space and is stretched out 0-1, 0-1 space when I normalize -- so the area has increased by a factor of 1/x^2.  this means the density must be multiplied by the same factor. *****NOT***** dividing by 2 here because we are integrating along \pi_{CB}, and the height of the desired area is 1/n here.
				}
			}
		}
		# so at this point we have cycled over the relevant area, and we scaled as we went. now just output it.
		if(list.please){return(list(sum = sum(fp.grid*sp.grid), fp.grid = fp.grid, sp.grid = sp.grid))}
		return(sum(fp.grid*sp.grid))
	}
}

# permutations of the Dirichlet special case, when we have truncated ballots
second.round.pivotal.event.probs.general = function(v.vec, s.vec, increments = 20, increments.at.end = 10){
	ab.thing = pr.second.round.pivotal.cands.1.and.2.general(v.vec, s.vec, increments = increments, increments.at.end = increments.at.end)
	ac.thing = pr.second.round.pivotal.cands.1.and.2.general(v.vec[c(2,1,5,6,3,4,7,9,8)], s.vec, increments = increments, increments.at.end = increments.at.end)
	bc.thing = pr.second.round.pivotal.cands.1.and.2.general(v.vec[c(4,3,6,5,1,2,8,9,7)], s.vec, increments = increments, increments.at.end = increments.at.end)
	list("AB" = ab.thing,
	"AC" = ac.thing,
	"BC" = bc.thing
	)
}



#####################################
#### first-round pivotal events #####
#####################################


# these functions are used below in first.round.pivotal.events.cands.1.and.2.truncated()
# we want to compute the dirichlet once over the whole possible area, and then add up the sums as appropriate.
dirichlet.matrix.to.beat.3 = function(alpha.segment, increments = 20){
	paxy = get.increment.midpoints(0, 1, increments) # \pi_{AX}/y, thus paxy
	pacy = get.increment.midpoints(0, 1, increments) # \pi_{AC}/y, thus pacy
	dir.mat = matrix(0, ncol = length(paxy), nrow = length(pacy)) # this matrix covers the simplex of \pi_{AX}/y, \pi_{AC}/y.
	colnames(dir.mat) = round(paxy, 4); rownames(dir.mat) = round(pacy, 4)
	# now calculate the dirichlet at each point in the grid
	for(i in 1:length(pacy)){
		for(j in 1:length(paxy)){
			eval.at = c(1 - paxy[j] - pacy[i], pacy[i], paxy[j])
			if(min(eval.at) <= 0 | max(eval.at) >= 1){next} # don't need to calculate when something is below or above 1.
			dir.mat[i,j] = ddirichlet(eval.at, alpha = alpha.segment)
		}
	}
	# normalize -- this should add up to 1. previously multiplied by the tile size
	dir.mat/sum(dir.mat) # *(paxy[2] - paxy[1])*(pacy[2] - pacy[1])
}

use.mat.to.beat.3 = function(y, increments = 20){
	# a matrix of integers: "For this y, should this combination of parameters be included in the Dirichlet sum?"
	paxy = get.increment.midpoints(0, 1, increments) # \pi_{AX}/y, thus paxy
	pacy = get.increment.midpoints(0, 1, increments) # \pi_{AC}/y, thus pacy
	use.mat = matrix(NA, ncol = 0, nrow = length(pacy))
	for(j in 1:length(paxy)){
		use.mat = cbind(use.mat, as.integer(pacy <= 2 - 1/(2*y) - paxy[j]/2))
	}
	use.mat
}


pr.first.round.pivotal.events.cands.1.and.2.general = function(v.vec, s.vec, increments = 20, increments.for.dir.mat = 20){
	stopifnot(length(v.vec) == 9)
	stopifnot(length(s.vec) == 4)
	fp.alpha = s.vec[1]*c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
	midpoints.z = get.increment.midpoints(1/4, 1/3, increments)  # the y values we will be evaluating
	pr.tie.for.second.at.z = pr.1.beats.3.given.z = pr.2.beats.3.given.z = c()
	sp.alpha.vec = c(v.vec[c(1,2)]*s.vec[2], v.vec[c(3,4)]*s.vec[3], v.vec[c(5,6)]*s.vec[4], v.vec[7:9]*s.vec[2:4])  # this allows the general case with second pref precision different from first pref precision
	# if truncated ballots are included, we calculate these to use below.
	if(v.vec[8] > 0){
		dir.mat.13 = dirichlet.matrix.to.beat.3(sp.alpha.vec[c(3,4,8)], increments.for.dir.mat) # this matrix shows the probability of getting entries in a simplex with BX along the horizontal and BC on the vertical.
	}
	if(v.vec[7] > 0){
		dir.mat.23 = dirichlet.matrix.to.beat.3(sp.alpha.vec[c(1,2,7)], increments.for.dir.mat)
	}
	for(z in midpoints.z){
		pr.tie.for.second.at.z = c(pr.tie.for.second.at.z, ddirichlet(c(z, z, 1-2*z), alpha = fp.alpha))
		# 1 beats 3
		if(v.vec[8] > 0 | v.vec[7] > 0){this.use.mat = use.mat.to.beat.3(z, increments.for.dir.mat)}
		if(v.vec[8] == 0){
			pr.1.beats.3.given.z = c(pr.1.beats.3.given.z, pbeta(2 - 1/(2*z), sp.alpha.vec[4], sp.alpha.vec[3]))
		}else{ # have to do the integration -- truncated ballots
			pr.1.beats.3.given.z = c(pr.1.beats.3.given.z, sum(this.use.mat*dir.mat.13))
		}
		if(v.vec[7] == 0){
			pr.2.beats.3.given.z = c(pr.2.beats.3.given.z, pbeta(2 - 1/(2*z), sp.alpha.vec[2], sp.alpha.vec[1]))
		}else{ # have to do the integration -- truncated ballots
			pr.2.beats.3.given.z = c(pr.2.beats.3.given.z, sum(this.use.mat*dir.mat.23))
		}
	}
	inc.size = (midpoints.z[2] - midpoints.z[1]) # edit this.
	list(
		"12.12" = inc.size*sum(pr.tie.for.second.at.z*pr.1.beats.3.given.z*pr.2.beats.3.given.z),  # NOT dividing by 2 here: for this Dirichlet, start from tie and add 1/n to the leader -- can do this at every point.
		"12.32" = inc.size*sum(pr.tie.for.second.at.z*(1 - pr.1.beats.3.given.z)*pr.2.beats.3.given.z),
		"12.13" = inc.size*sum(pr.tie.for.second.at.z*pr.1.beats.3.given.z*(1 - pr.2.beats.3.given.z))
	)
}

first.round.pivotal.event.probs.general = function(v.vec, s.vec, increments = 100, increments.for.dir.mat = 20){
	frpe.12 = pr.first.round.pivotal.events.cands.1.and.2.general(v.vec, s.vec, increments = increments, increments.for.dir.mat = increments.for.dir.mat)
	frpe.13 = pr.first.round.pivotal.events.cands.1.and.2.general(v.vec = v.vec[c(2,1,5,6,3,4,7,9,8)], s.vec = s.vec, increments = increments, increments.for.dir.mat = increments.for.dir.mat)
	frpe.23 = pr.first.round.pivotal.events.cands.1.and.2.general(v.vec = v.vec[c(4,3,6,5,1,2,8,9,7)], s.vec = s.vec, increments = increments, increments.for.dir.mat = increments.for.dir.mat)
	list("AB.AB" = frpe.12[["12.12"]], "AB.AC" = frpe.12[["12.13"]], "AB.CB" = frpe.12[["12.32"]],
	"AC.AC" = frpe.13[["12.12"]], "AC.AB" = frpe.13[["12.13"]], "AC.BC" = frpe.13[["12.32"]],
	"BC.BC" = frpe.23[["12.12"]], "BC.BA" = frpe.23[["12.13"]], "BC.AC" = frpe.23[["12.32"]])
}




