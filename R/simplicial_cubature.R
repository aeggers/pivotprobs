simplicial_cubature_based_plurality_pivot_probs <- function(alpha = NULL, mu = NULL, precision = NULL, sigma = NULL, cand_names = NULL, sep = "", report_issues = F, full_output = F, ...){

  if(!is.null(alpha) | (!is.null(mu) & !is.null(precision)) & is.null(sigma)){
    # this is Dirichlet
    distribution = "dirichlet"
    if(is.null(alpha)){alpha <- mu*precision}
    k <- length(alpha)
  }else if(!is.null(mu) & !is.null(sigma)){
    # this is logistic normal
    distribution = "logisticnormal"
    k <- length(mu)
  }

  if(is.null(cand_names)){cand_names <- letters[1:k]}

  S <- plurality_simplices_to_integrate(k)

  out <- list()
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      these_indices <- c(i,j,(1:k)[-c(i,j)])
      this_name <- paste0(cand_names[i], sep, cand_names[j])
      if(distribution == "dirichlet"){
        out[[this_name]] = SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = S, alpha = alpha[these_indices], ...)
      }else if(distribution == "logisticnormal"){
        out[[this_name]] = SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = S, mu = mu[these_indices], sigma = sigma[these_indices, these_indices], ...)
      }
    }
  }

  if(report_issues){
    msgs <- out %>% map("message")
    if(length(unique(msgs)) != 1){
      cat("Messages from SimplicialCubature::adaptIntegrateSimplex: \n", paste(msg, collapse = "\n"))
    }
  }

  if(full_output){
    out
  }else{
    out %>% map("integral")
  }

}


cand_a_win_region_vertices_from_win_conditions <- function(win_conditions){

  # returns a matrix of vertices for candidate a's win region (one vertex per row) given conditions under which candidate a wins.

  # win_conditions states conditions for the ballot shares such that candidate a wins. conditions are stated in terms of coefficients $\beta$ on the ballot shares $\mathbf{v}$ such that a wins when $\mathbf{v} \beta \geq 0$ e.g. three-candidate plurality we would have rbind(c(1, -1, 0), c(1, 0, -1))

  # we combine this with conditions that place the ballot shares on the unit simplex:
  #   inequalities: v_a > 0, v_b > 0, . . .
  #   equation (a2 and b2 arguments to rcdd::makeH): \sum{\mathbf{v}}  = 1

  # and we reverse the sign in the inequalities so that we can apply these to the rcdd::makeH function
  positivity_part <- -diag(ncol(win_conditions))

  a1 <- rbind(-win_conditions, positivity_part)
  # a1_char <- matrix(as.character(a1), nrow = nrow(a1), ncol = ncol(a1))

  # represents the polyhedron described by the win_conditions and conditions on vote shares.
  cand_a_win_region_H <- rcdd::makeH(a1 = a1, b1 = rep(0, nrow(a1)), a2 = rep(1, ncol(win_conditions)), b2 = 1)

  # gets the convex hull of points from H representation
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

  # wrv is already the convex hull of the win region.
  # but convexhulln() triangulates this, i.e. it describes the convex hull in terms of simplices
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


# combining these

S_array_from_win_conditions <- function(win_conditions){
  wrv <- cand_a_win_region_vertices_from_win_conditions(win_conditions)
  simplices_to_integrate_from_win_region_vertices(wrv, binding_constraint = win_conditions[1, ])
}


plurality_simplices_to_integrate <- function(k){
  S_array_from_win_conditions(cbind(rep(1, k-1), -diag(k-1)))
}

# not used
simplicial_cubature_plurality_pivot_probs_dirichlet <- function(alpha, cand_names = NULL, sep = "", full_output = F, ...){

  out <- list()
  k <- length(alpha)
  S <- plurality_simplices_to_integrate(k)
  if(is.null(cand_names)){cand_names <- letters[1:k]}

  all_indices <- 1:k

  for(i in 1:(k-1)){
    for(j in (i+1):k){
      these_indices <- c(i,j,all_indices[-c(i,j)])
      out[[paste0(cand_names[i], sep, cand_names[j])]] = SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = S, alpha = alpha[these_indices], ...)
    }
  }

  if(report_issues){
    msgs <- out %>% map("message")
    if(length(unique(msgs)) != 1){
      cat("Messages from SimplicialCubature::adaptIntegrateSimplex: \n", paste(msg, collapse = "\n"))
    }
  }

  if(full_output){
    out
  }else{
    out %>% map("integral")
  }

}

# not used
simplicial_cubature_plurality_pivot_probs_logisticnormal <- function(mu, sigma, cand_names = NULL, sep = "", full_output = F, report_issues = T, normalizing_factor = 1, ...){

  out <- list()
  k <- length(mu) + 1 # key difference from dirichlet here.
  S <- plurality_simplices_to_integrate(k)
  if(is.null(cand_names)){cand_names <- letters[1:(k+1)]}

  all_indices <- 1:k

  for(i in 1:(k-1)){
    for(j in (i+1):k){
      these_indices <- c(i,j,all_indices[-c(i,j)])
      the_name <- paste0(cand_names[i], sep, cand_names[j])
      out[[the_name]] = SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = S[these_indices,,], mu = mu, sigma = sigma, ...) # note that I spin around the S rather than the mu and sigma
      out[[the_name]]$integral <- out[[the_name]]$integral/normalizing_factor
    }
  }

  if(report_issues){
    msgs <- out %>% map("message")
    if(length(unique(msgs)) != 1){
      cat("Messages from SimplicialCubature::adaptIntegrateSimplex: \n", paste(msg, collapse = "\n"))
    }
  }

  if(full_output){
    out
  }else{
    out %>% map("integral")
  }

}


## integrand functions

dirichlet_for_integration <- function(x, alpha){
  out <- gtools::ddirichlet(as.vector(x), alpha)/sqrt(length(alpha)) # normalizing constant so that is sums to 1 on the unit simplex
  ifelse(is.nan(out) | is.infinite(out), 0, out)
}

logisticnormal_for_integration <- function(x, mu, sigma){
  out <- dlogisticnormal(as.vector(x), mu, sigma)/sqrt(length(mu)) # normalizing constant so that is sums to 1 on the unit simplex
  ifelse(is.nan(out) | is.infinite(out), 0, out)
}


# #to see need for this correction:
# naive_dirichlet_fn <- function(x, alpha){
#   gtools::ddirichlet(as.vector(x), alpha)
# }

# SimplicialCubature::adaptIntegrateSimplex(naive_dirichlet_fn, S = diag(3), alpha = c(5,4,3))  # returns sqrt(3)
# SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = diag(3), alpha = c(5,4,3)) returns 1

# SimplicialCubature::adaptIntegrateSimplex(naive_dirichlet_fn, S = diag(5), alpha = c(5,4,3, 2, 2), tol = .1, maxEvals = 100000)  # returns sqrt(5)
# SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = diag(5), alpha = c(5,4,3, 2,2), tol = .1, maxEvals = 100000) returns 1

# this was not a good idea.
# mvnorm_for_integration <- function(x, mu, sigma){
#   mvtnorm::dmvnorm(as.numeric(x), mean = mu, sigma = sigma)
# }

# SimplicialCubature::adaptIntegrateSimplex(mvnorm_for_integration, S = diag(3), mu = c(.4, .35, .25), sigma = diag(3)*.02)

