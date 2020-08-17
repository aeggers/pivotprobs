## general approach

cand_a_win_region_vertices_from_win_conditions <- function(win_conditions){

  positivity_part <- -diag(ncol(win_conditions))

  a1 <- rbind(-win_conditions, positivity_part)

  cand_a_win_region_H <- makeH(a1 = a1, b1 = rep(0, nrow(a1)), a2 = rep(1, ncol(win_conditions)), b2 = 1)

  out <- rcdd::scdd(cand_a_win_region_H)

  out$output[, -c(1,2)]
}

# test with plurality
wrv3 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0), c(1, 0, -1)))
wrv4 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0,0), c(1, 0, -1,0), c(1, 0, 0, -1)))
# borda count
s <- 1/2
bc_wcs <- rbind(c(1-s, 1, s - 1, -1, s, -s),
            c(1, 1-s, s, -s, s - 1, -s))
wrv_bc <- cand_a_win_region_vertices_from_win_conditions(bc_wcs)


simplices_to_integrate_from_win_region_vertices <- function(wrv, binding_constraint = c(1,-1,0,0)){

  # which vertices in wrv have an ab tie?
  # TODO: worry about floating point here
  vertices_with_ab_tie <- which(as.vector(wrv %*% matrix(binding_constraint, ncol = 1)) == 0)

  # get matrix describing convex hull of a's win region:
  # one row per simplex, one column per vertex
  tch <- geometry::convhulln(wrv[,-ncol(wrv)])

  # we want to select the simplices on the relevant facets.
  # these simplices will *only* have vertices where a ties b.
  # so pivot longer, making the simplex the group
  data.frame(tch) %>%
    mutate(simplex = 1:nrow(.)) %>%
    pivot_longer(cols = starts_with("X"), values_to = "vertex", names_to = "name") %>%
    select(-name) %>%
    group_by(simplex) %>%
    # and select simplices where all vertices are in vertices_with_ab_tie
    filter(sum(!vertex %in% vertices_with_ab_tie) == 0) -> simplex_vertex

  wrv_df <- data.frame(vertex = 1:nrow(wrv), wrv[,-ncol(wrv)])

  simplex_vertex %>%
    left_join(wrv_df, by = "vertex") %>%
    select(-vertex) %>% as.matrix() -> simplices

  #  )%>%
  #  nest() -> nested_simplices
  matrix_list <- list()
  ss <- unique(simplices[,"simplex"])
  for(s in ss){
    matrix_list[[as.character(s)]] = t(simplices[which(simplices[,"simplex"] == s), -1])
  }

  array(matrix_list %>% unlist(), dim = c(nrow(matrix_list[[1]]), ncol(matrix_list[[1]]), length(ss)))
}

# test
simplices_to_integrate_from_win_region_vertices(wrv4, binding_constraint = c(1,-1,0,0))

bc_sti <- simplices_to_integrate_from_win_region_vertices(wrv_bc, binding_constraint = bc_wcs[1,])

SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = bc_sti, alpha = c(8, 6, 5, 6, 4, 9), maxEvals = 100000L, absError = .1)

SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = bc_sti, alpha = c(8, 6, 5, 6, 4, 9), maxEvals = 100000L, absError = .1, integRule = 2L) # that doesn't help!

# okay this is worrying -- taking a long time. I wonder if there is some way to speed this up. maybe not.

# good: we have a procedure. bad: it takes forever.


# combining these a bit
plurality_simplices_to_integrate <- function(k){
  win_conditions <- cbind(rep(1, k-1), -diag(k-1))

  wrv_k <- cand_a_win_region_vertices_from_win_conditions(win_conditions)

  simplices_to_integrate_from_win_region_vertices(wrv_k, binding_constraint = win_conditions[1, ])
}

plurality_simplices_to_integrate(3)
plurality_simplices_to_integrate(4)
plurality_simplices_to_integrate(5)
plurality_simplices_to_integrate(6)

# a little unclear what we want here:
#TODO: work this out a bit. will it always be the case that we can re-label the candidates and get the right answer by focusing on a vs b? with Dirichlet, this is just reshuffling the alpha vector. If it were e.g. normal distribution, we can transform the mean and variance-covariance matrix. so I think it's okay to focus on this approach.

# so next we put things together:
# supply the conditions for a to win (with the condition we want to satisfy with equality being the first condition); get the S array from that (using above); then iterate through Dirichlet parameters to integrate this for each candidate. (generalize to other belief distributions later.)
# SimplicialCubature::adaptIntegrateSimplex

integrate_over_simplices <- function(f, S, tol = .0001, ...){
  SimplicialCubature::adaptIntegrateSimplex(f, S, tol, ...)
}

dirichlet_for_integration <- function(x, alpha){
  gtools::ddirichlet(c(as.numeric(x), 1 - sum(x)), alpha)
}

mvnorm_for_integration <- function(x, mu, sigma){
  mvtnorm::dmvnorm(c(as.numeric(x), 1 - sum(x)), mean = mu, sigma = sigma)
}

## now knit together for a voting method

# TODO: generalize a bit? other methods, other distributions?
plurality_pivot_probs_dirichlet <- function(alpha, cand_names = NULL, sep = "", tol = .0001, ...){

  k <- length(alpha)

  if(is.null(cand_names)){
    cand_names <- names(alpha)
    if(is.null(cand_names)){
      cand_names <- letters[1:k]
    }
  }

  S <- plurality_simplices_to_integrate(k)

  out <- list()

  for(i in 1:(k - 1)){
    for(j in (i+1):k){
      out[[paste0(cand_names[i], sep, cand_names[j])]] = SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = S, alpha = c(alpha[c(i,j)], alpha[-c(i,j)]), tol = tol, ...) #
    }
  }

  out

}

# let's test this for three candidates
# one piv prob:
SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = cbind(c(.5, .5), c(1/3, 1/3)), alpha = c(10, 8, 4))
SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = plurality_simplices_to_integrate(3), alpha = c(10, 8, 4))

# it's fine if it's just one array.
SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = array(cbind(c(.5, .5), c(1/3, 1/3)), dim = c(2,2,1)), alpha = c(10, 8, 4))

p3 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5))
p4 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4))
p5 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3)) # takes a long time

# what if we make it less precise?

system.time(p1 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3), absError = .1)) # less than a second
system.time(p01 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3), absError = .01))
system.time(p001 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3), absError = .001)) # four seconds.

system.time(p1 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3), tol = .1)) # 1.5 seconds
system.time(p01 <- plurality_pivot_probs_dirichlet(alpha = c(10, 8, 5, 4, 3), tol = .001))  # 9 seconds.

# relative seems better than absolute error.

# the speed is definitely a problem. disappointing.

# TODO: next steps: write it up for other methods, see if the performance is bad. it might not be -- maybe there aren't as many facets for e.g. Borda?

