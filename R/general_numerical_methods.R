## general approach

dirichlet_for_integration <- function(x, alpha){
  gtools::ddirichlet(as.vector(x), alpha)/sqrt(length(alpha)) # correction for dimensionality
}

# could integrate other functions

cand_a_win_region_vertices_from_win_conditions <- function(win_conditions){

  # returns a matrix of vertices for candidate a's win region (one vertex per row) given conditions under which candidate a wins.

  # win_conditions states conditions for the ballot shares such that candidate a wins. conditions are stated in terms of coefficients $\beta$ on the ballot shares $\mathbf{v}$ such that a wins when $\mathbf{v} \beta \geq 0$ e.g. three-candidate plurality we would have rbind(c(1, -1, 0), c(1, 0, -1))

  # we combine this with conditions that place the ballot shares on the unit simplex:
    #   inequalities: v_a > 0, v_b > 0, . . .
    #   equation (a2 and b2 arguments to rcdd::makeH): \sum{\mathbf{v}}  = 1

  # and we reverse the sign in the inequalities so that we can apply these to the rcdd::makeH function
  positivity_part <- -diag(ncol(win_conditions))

  a1 <- rbind(-win_conditions, positivity_part)

  cand_a_win_region_H <- rcdd::makeH(a1 = a1, b1 = rep(0, nrow(a1)), a2 = rep(1, ncol(win_conditions)), b2 = 1)

  out <- rcdd::scdd(cand_a_win_region_H)

  out$output[, -c(1,2)] ## I am not sure what the first two columns are for.
}

test <- F

if(test){
  # test with plurality
  wrv3 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0), c(1, 0, -1)))
  wrv4 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0,0), c(1, 0, -1,0), c(1, 0, 0, -1)))

  # borda count: seems correct.
  s <- 1/2
  bc_wcs <- rbind(c(1-s, 1, s - 1, -1, s, -s),
                  c(1, 1-s, s, -s, s - 1, -1))
  wrv_bc <- cand_a_win_region_vertices_from_win_conditions(bc_wcs)
}



simplices_to_integrate_from_win_region_vertices <- function(wrv, binding_constraint, epsilon = 1.0e-15){

  ## returns an array (the S array for SimplicialCubature::adaptIntegrateSimplex) describing the simplices that cover the facet(s) of the win region described by wrv where binding_constraint applies.

  # identify the vertices where binding_constraint applies
  vertices_with_ab_tie <- which(abs(as.vector(wrv %*% matrix(binding_constraint, ncol = 1))) < epsilon)

  if(length(vertices_with_ab_tie) == 0){stop("binding_constraint does not hold at any vertices in the win_region_vertices.")}

  # get matrix describing convex hull of a's win region -- crucially, this is triangulated, i.e. in terms of simplices, so it's a good start for forming the S matrix.
  # one row per simplex, one column per vertex
  # we drop a dimension when we do this, but it's fine -- we bring it back later
  tch <- geometry::convhulln(wrv[,-ncol(wrv)])

  # we want to select the simplices on the relevant facets.
  # these simplices will *only* have vertices where a ties b, i.e. where the binding_constraint holds.
  # so pivot longer, making the simplex the group
  data.frame(tch) %>%
    mutate(simplex = 1:nrow(.)) %>%
    pivot_longer(cols = starts_with("X"), values_to = "vertex", names_to = "name") %>%
    select(-name) %>%
    group_by(simplex) %>%
    # and select simplices where all vertices are in vertices_with_ab_tie
    filter(sum(!vertex %in% vertices_with_ab_tie) == 0) -> simplex_vertex

  wrv_df <- data.frame(vertex = 1:nrow(wrv), wrv) # [,-ncol(wrv)])

  simplex_vertex %>%
    left_join(wrv_df, by = "vertex") %>%
    select(-vertex) %>% as.matrix() -> simplices

  # in the S array we need for SimplicialCubature, the columns are the vertices of the simplices over which we integrate. (note that's not the usual way R works.)
  matrix_list <- list()
  ss <- unique(simplices[,"simplex"])
  for(s in ss){
    matrix_list[[as.character(s)]] = t(simplices[which(simplices[,"simplex"] == s), -1])
  }

  array(matrix_list %>% unlist(), dim = c(nrow(matrix_list[[1]]), ncol(matrix_list[[1]]), length(ss)))
}

test <- F

if(test){
  # test
  simplices_to_integrate_from_win_region_vertices(wrv4, binding_constraint = c(1,-1,0,0))

  bc_sti <- simplices_to_integrate_from_win_region_vertices(wrv_bc, binding_constraint = bc_wcs[1,])

  system.time(out <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = bc_sti, alpha = c(8, 6, 5, 6, 4, 9), maxEvals = 100000L, absError = .1))  # 5 seconds to do a single piv prob -- that's not good.
}


# good: we have a procedure. bad: it takes forever.


# combining these a bit

S_array_from_win_conditions <- function(win_conditions){
  wrv <- cand_a_win_region_vertices_from_win_conditions(win_conditions)
  simplices_to_integrate_from_win_region_vertices(wrv, binding_constraint = win_conditions[1, ])
}


plurality_simplices_to_integrate <- function(k){
  win_conditions <- cbind(rep(1, k-1), -diag(k-1))

  wrv_k <- cand_a_win_region_vertices_from_win_conditions(win_conditions)

  simplices_to_integrate_from_win_region_vertices(wrv_k, binding_constraint = win_conditions[1, ])
}

test <- F
if(test){
  plurality_simplices_to_integrate(3)
  plurality_simplices_to_integrate(4)
  plurality_simplices_to_integrate(5)
  plurality_simplices_to_integrate(6)
}



# a little unclear what we want here:
#TODO: work this out a bit. will it always be the case that we can re-label the candidates and get the right answer by focusing on a vs b? with Dirichlet, this is just reshuffling the alpha vector. If it were e.g. normal distribution, we can transform the mean and variance-covariance matrix. so I think it's okay to focus on this approach.

# so next we put things together:
# supply the conditions for a to win (with the condition we want to satisfy with equality being the first condition); get the S array from that (using above); then iterate through Dirichlet parameters to integrate this for each candidate. (generalize to other belief distributions later.)
# SimplicialCubature::adaptIntegrateSimplex


# this would need some kind of normalization I guess?
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

test <- F
if(test){
  alpha <- c(8,6,4,2)
  system.time(out <- plurality_pivot_probs_dirichlet(alpha))
  out_sim <- plurality_pivot_probs_simulation_based(alpha = alpha)
  out_sim_2 <- plurality_pivot_probs_simulation_based(alpha = alpha, N = 1000000)
  out_sim_3 <- plurality_pivot_probs_simulation_based(alpha = alpha, N = 10000000) # okay that takes a long time
  sc_vec <- c(out$ab$integral, out$ac$integral, out$ad$integral, out$bc$integral, out$bd$integral, out$cd$integral)
  sc_vec/unlist(out_sim)
  sc_vec/unlist(out_sim_2)
  sc_vec/unlist(out_sim_3)
}



ordinal_shuffle_dirichlet_pivot_probs <- function(alpha, S, cand_names = NULL, sep = "", ...){

  # for now, constrained to three candidates
  stopifnot(length(alpha) == 6)

  if(is.null(cand_names)){
    cand_names <- names(alpha)
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
  }

  out <- list()

  out[[paste0(cand_names[1], sep, cand_names[2])]] = SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = S, alpha = alpha, ...)

  out[[paste0(cand_names[1], sep, cand_names[3])]] = SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = S, alpha = alpha[c(2,1,5,6,3,4)], ...)

  out[[paste0(cand_names[2], sep, cand_names[3])]] = SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = S, alpha = alpha[c(4,3,6,5,1,2)], ...)

  out

}

positional_pivot_probs_general <- function(s = 1/2, alpha, cand_names = NULL, sep = "", tol = .01, ...){

  wcs <- rbind(c(1-s, 1, s - 1, -1, s, -s),
                  c(1, 1-s, s, -s, s - 1, -1))
  wrv <- cand_a_win_region_vertices_from_win_conditions(wcs)

  S <- simplices_to_integrate_from_win_region_vertices(wrv, binding_constraint = wcs[1,])

  ordinal_shuffle_dirichlet_pivot_probs(alpha, S, cand_names = cand_names, sep = sep, tol = tol, ...)

}

if(test){
  alpha_vec <- c(8,7,4,3,5,8)
  system.time(bc_out <- positional_pivot_probs_general(alpha = alpha_vec, maxEvals = 100000))
  # and then compare with direct MC?
  bc_sim <- positional_piv_probs_simulation(alpha_vec = alpha_vec)
  bc_sim_2 <- positional_piv_probs_simulation(alpha_vec = alpha_vec, N = 1000000)
  bc_sim_3 <- positional_piv_probs_simulation(alpha_vec = alpha_vec, N = 10000000)

  c(bc_out$ab$integral, bc_out$ac$integral, bc_out$bc$integral)/c(bc_sim$ab, bc_sim$ac, bc_sim$bc)
  c(bc_out$ab$integral, bc_out$ac$integral, bc_out$bc$integral)/c(bc_sim_2$ab, bc_sim_2$ac, bc_sim_2$bc)
  c(bc_out$ab$integral, bc_out$ac$integral, bc_out$bc$integral)/c(bc_sim_3$ab, bc_sim_3$ac, bc_sim_3$bc)

  # so they don't even really agree. WTF.
  # let's test this for three candidates
  # one piv prob:
  SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = cbind(c(.5, .5), c(1/3, 1/3)), alpha = c(10, 8, 4))
  SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = plurality_simplices_to_integrate(3), alpha = c(10, 8, 4))

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
}


# TODO: next steps: write it up for other methods, see if the performance is bad. it might not be -- maybe there aren't as many facets for e.g. Borda?

