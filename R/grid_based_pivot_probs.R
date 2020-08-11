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

PTG_store <- list()

# NB: this could be made not specific to Dirichlet.
unnormalized_pivot_prob_from_ptg_and_alpha <- function(ptg, alpha){
  sum(gtools::ddirichlet(as.matrix(ptg), alpha))
}

# NB: this could be made not specific to Dirichlet.
#' @export
plurality_pivot_probs_grid_based <- function(alpha_vec, increment = .025, cand_names = NULL, sep = ""){
  if(is.null(cand_names)){
    cand_names <- names(alpha_vec)
  }
  if(is.null(cand_names)){
    stop("Must provide either cand_names argument or named alpha_vec.")
  }
  k <- length(alpha_vec)
  # TODO: allow for these to be stored so I don't have to regenerate.
  this_ptg_name <- paste0("k_", k, "_increment_", increment)
  if(this_ptg_name %in% names(PTG_store)){
    # cat("using stored ptg!\n")
    ptg <- PTG_store[[this_ptg_name]]
  }else{
    # cat("making ptg from scratch!\n")
    PTG_store[[this_ptg_name]] <<- ptg <- plurality_tie_grid(increment, k)
  }
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
# There are more efficient ways to do this that take advantage of Dirichlet properties: the approximation in Eggers & Vivyan; probably also better to do something that takes advantage of Dirichlet properties (aggregation).
# But this could be made not specific to Dirichlet.
# TODO:
  ## -- compare with version based on approximation (used in Eggers and Vivyan), and simulation.
  ## -- provide a way to grab grids from storage so we don't have to make them over and over again.
  ## -- apply to other methods -- this is the real key, as it's easier to do this than to rework the math for each case.
  ## -- allow other distributions (non-Dirichlet), another attractive generalization


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








