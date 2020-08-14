### grid-based estimation
###

# TODO: this could be made not specific to Dirichlet. pass a function and parameters
# is not specific to plurality.
unnormalized_pivot_prob_from_tie_grid_and_alpha <- function(tie_grid, alpha){
  sum(gtools::ddirichlet(as.matrix(tie_grid), alpha))
}

# now for each type of pivot probability, a tie_grid:

plurality_tie_grid <- function(increment = .025, k = 4){
  # this one is base R
  stopifnot(k >= 3)
  # beware using a small increment and high k because that will take a long time.
  increments <- round(1/increment, 0)
  effective_increment <- 1/increments
  standard_seq <- seq(0 + effective_increment/2, 1 - effective_increment/2, length = increments)

  to_expand <- as.data.frame(matrix(standard_seq, ncol = k - 2, nrow = length(standard_seq), byrow = F))

  the_grid <- expand.grid(to_expand)
  the_grid$V0 <- the_grid$V1
  the_grid$last <- 1 - apply(the_grid, 1, sum)
  the_grid <- the_grid[which(the_grid$last >= effective_increment/2),]
  # slight fudge to deal with floating point stuff, which was making slightly differently shaped grids depending on the increment.
  the_grid <- the_grid[which(the_grid$V0 > apply(the_grid, 1, max) - .00001),]
  the_grid[,c("V0", "V1", names(the_grid)[which(!names(the_grid) %in% c("V0", "V1", "last"))], "last")]

}

positional_tie_grid <- function(increment = .025, positional_s = .5){
  # this one uses dplyr.
  # can probably make this a little more efficient by pruning progressively, e.g. expand ab and ba, take out ab + ba > 1, then expand further.
  increments <- round(1/increment, 0)
  effective_increment <- 1/increments
  standard_seq <- seq(0 + effective_increment/2, 1 - effective_increment/2, length = increments)

  expand_grid(
    ab = standard_seq,
    ba = standard_seq,
    ca = standard_seq,
    cb = standard_seq
  ) %>%
    mutate(
      ac = .5*(1 - positional_s*ba + (positional_s-2)*ab - (1 + positional_s)*ca - (1 - positional_s)*cb),
      bc = 1 - ab - ac - ba - ca - cb,
      score_a = ab + ac + positional_s*(ba + ca),
      score_c = ca + cb + positional_s*(ac + bc)
    ) %>%
    # we require bc and ac to be within the bounds of the main variables so that we don't get blowups when there are small alpha values.
    filter(bc >= effective_increment/2 & bc <= 1 - effective_increment/2 & ac >= effective_increment/2 & ac <= 1  - effective_increment/2 & score_c < score_a) %>% # filtering on ac not necessary I think.
    select(ab, ac, ba, bc, ca, cb)
}

# example:
# ptg5 <- plurality_tie_grid(k = 5, incr = .01)

TIE_GRID_STORE <- list()


# NB: this could be made not specific to Dirichlet.
#' @export
plurality_pivot_probs_grid_based <- function(alpha_vec, increment = .025, cand_names = NULL, sep = "", boost_low_alphas = T){
  if(min(alpha_vec) < 1 & !boost_low_alphas){cat("Note: alpha_vec contains values below 1 -- beware infinite density near edges, and consider boosting low alphas.")}
  if(is.null(cand_names)){
    cand_names <- names(alpha_vec)
  }
  if(is.null(cand_names)){
    stop("Must provide either cand_names argument or named alpha_vec.")
  }
  k <- length(alpha_vec)
  this_ptg_name <- paste0("k_", k, "_increment_", increment)
  if(this_ptg_name %in% names(TIE_GRID_STORE)){
    ptg <- TIE_GRID_STORE[[this_ptg_name]]
  }else{
    TIE_GRID_STORE[[this_ptg_name]] <<- ptg <- plurality_tie_grid(increment, k)
  }
  effective_increment <- 1/round(1/increment, 0)
  out <- list()
  for(i in 1:(length(cand_names)-1)){
    for(j in (i+1):length(cand_names)){
      this_alpha <- c(alpha_vec[c(i,j)], alpha_vec[-c(i,j)])
      if(min(this_alpha[-c(1,2)]) < 1 & boost_low_alphas){
        cat("Boosting low alphas.\n")
        to_boost <- this_alpha[-c(1,2)]
        to_boost[to_boost < 1] <- 1
        this_alpha <- c(this_alpha[1:2], to_boost)
      }
      out[[paste0(cand_names[i], sep, cand_names[j])]] = unnormalized_pivot_prob_from_tie_grid_and_alpha(ptg, this_alpha)*(effective_increment)^(k - 2)
    }
  }

  out

}

# example:
# alpha_vec <- c(.4, .35, .20, .05)*20
# plurality_pivot_probs_grid_based(alpha_vec, increment = .01)
# compare with
# sims <- gtools::rdirichlet(10000000, alpha = alpha_vec)
# pps <- pivotal.probabilities(sims)


# this is repetitive but not easy to separate out
positional_pivot_probs_grid_based <- function(alpha_vec, positional_s = .5, increment = .025, cand_names = c("a", "b", "c"), sep = "", boost_low_alphas = T){
  if(length(alpha_vec) != 6){
    stop("For positional methods we expect three candidates and six shares.")
  }
  this_ptg_name <- paste0("positional_s_", positional_s, "_increment_", increment)
  if(this_ptg_name %in% names(TIE_GRID_STORE)){
    ptg <- TIE_GRID_STORE[[this_ptg_name]]
  }else{
    TIE_GRID_STORE[[this_ptg_name]] <<- ptg <- positional_tie_grid(increment, positional_s)
  }
  # TODO: explain why we need to divide by 2 -- I think it's because of Borda
  effective_increment <- 1/round(1/increment, 0)
  if(min(alpha_vec) < 1){
    if(!boost_low_alphas){
      cat("Some elements of alpha_vec are below 1. Consider boosting.\n")
      }else{
      alpha_vec[alpha_vec < 1] <- 1
    }
  }
  out <- list(
    unnormalized_pivot_prob_from_tie_grid_and_alpha(ptg, alpha_vec)*effective_increment^4/2,
    unnormalized_pivot_prob_from_tie_grid_and_alpha(ptg, alpha_vec[c(2,1,5,6,3,4)])*effective_increment^4/2,
    unnormalized_pivot_prob_from_tie_grid_and_alpha(ptg, alpha_vec[c(4,3,6,5,1,2)])*effective_increment^4/2
  )

  names(out) <- c(paste0(cand_names[1], sep, cand_names[2]),
                  paste0(cand_names[1], sep, cand_names[3]),
                  paste0(cand_names[2], sep, cand_names[3]))

  out

}

# example:
# alpha_vec <- c(.2, .1, .1, .1, .25, .25)*20
# pps_gb <- positional_pivot_probs_grid_based(alpha_vec, positional_s = .5, increment = .02)
# compare with
# sims <- gtools::rdirichlet(500000, alpha_vec)
# pps <- positional_piv_probs_simulation(sims)

# pps_gb_1 <- positional_pivot_probs_grid_based(alpha_vec, positional_s = 0, increment = .02)
# sims <- gtools::rdirichlet(500000, c(sum(alpha_vec[1:2]), sum(alpha_vec[3:4]), sum(alpha_vec[5:6])))
# colnames(sims) <- c("a", "b", "c")
# pps_1 <- pivotal.probabilities(sims)

# consistently twice as large.

# pps_gb_2 <- positional_pivot_probs_grid_based(alpha_vec, positional_s = .5, increment = .03)
# compare with pps
# rbind(pps_gb_2 %>% unlist(), pps)

# pps_gb_3 <- positional_pivot_probs_grid_based(alpha_vec, positional_s = 1, increment = .02)
# sims <- gtools::rdirichlet(500000, alpha_vec)
# pps_3 <- positional_piv_probs_simulation(sims, s = 1)
# rbind(pps_gb_3 %>% unlist(), pps_3)



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








