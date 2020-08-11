### grid-based estimation
###

# same thing but not pruning at every opportunity
plurality_tie_grid <- function(increment = .025, k = 4){
  stopifnot(k >= 3)
  # beware using a small increment and high k
  standard_seq <- seq(increment/2, .5 - increment/2, by = increment)

  to_expand <- as.data.frame(matrix(standard_seq, ncol = k - 2, nrow = length(standard_seq), byrow = F))

  the_grid <- expand.grid(to_expand)
  the_grid$V0 <- the_grid$V1
  the_grid <- the_grid[which(apply(the_grid, 1, sum) < 1), ]
  the_grid$last <- 1 - apply(the_grid, 1, sum)
  the_grid <- the_grid[which(the_grid$V0 == apply(the_grid, 1, max)),]
  the_grid[,c("V0", "V1", names(the_grid)[which(!names(the_grid) %in% c("V0", "V1", "last"))], "last")]
#   the_grid[,c(ncol(the_grid)-1,1,2:(ncol(the_grid) -2), ncol(the_grid))]
}

# example:
# ptg5 <- plurality_tie_grid(k = 5, incr = .01)

unnormalized_pivot_prob_from_ptg_and_alpha <- function(ptg, alpha){
  sum(gtools::ddirichlet(as.matrix(ptg), alpha))
}

plurality_pivot_probs_grid_based <- function(alpha_vec, increment = .025, cand_names = NULL, sep = ""){
  if(is.null(cand_names)){
    cand_names <- names(alpha_vec)
  }
  if(is.null(cand_names)){
    stop("Must provide either cand_names argument or named alpha_vec.")
  }
  k <- length(alpha_vec)
  ptg <- plurality_tie_grid(increment, k)
  # cat(dim(ptg))
  out <- list()
  for(i in 1:(length(cand_names)-1)){
    for(j in (i+1):length(cand_names)){
      this_alpha <- c(alpha_vec[c(i,j)], alpha_vec[-c(i,j)])
      out[[paste0(cand_names[i], sep, cand_names[j])]] = unnormalized_pivot_prob_from_ptg_and_alpha(ptg, this_alpha)*increment^(k - 2)
    }
  }

  out

}

# example: plurality_pivot_probs_grid_based(alpha_vec*20, increment = .0025)

# So implementation-wise this seems to be working.
# next to do: compare with version based on approximation (used in Eggers and Vivyan), and simulation.
# and if it's working, can apply to other methods.


# But:
  ## the increment matters a lot, in some cases when the runtime also differs
  ## e.g.
  ## plurality_pivot_probs_grid_based(c(.35, .3, .2, .1, .05)*20, increment = .01, cand_names = c("A", "B", "C", "D", "E"), sep = "_")

  ## $A_B
  ## [1] 1.516363

  ## $A_C
  ## [1] 0.5390683

  ## $A_D
  ## [1] 0.08465921

  # etc

  ## plurality_pivot_probs_grid_based(c(.35, .3, .2, .1, .05)*20, increment = .025, cand_names = c("A", "B", "C", "D", "E"), sep = "_")
  ## $A_B
  ## [1] 1.348127

  ## $A_C
  ## [1] 0.4946683

  ## $A_D
  ## [1] 0.07917274

  # although the ratios are closer. and close is good enough for this.








