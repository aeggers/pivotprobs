# simplicial cubature utils

make_inequalities_into_near_equalities <- function(inequality_mat, to_make_near_equality = c(NULL), n = 5000){
  a_mat <- inequality_mat[ ,-ncol(inequality_mat)]
  b_vec <- inequality_mat[ ,ncol(inequality_mat)]

  b_vec[to_make_near_equality] <- 0
  a_mat <- rbind(a_mat, -a_mat[to_make_near_equality, ]) # add some \leq conditions
  b_vec <- c(b_vec, rep(-1/n, length(to_make_near_equality)))

  cbind(a_mat, b_vec)
}

vertices_of_integration_region_from_inequalities <- function(inequality_mat){

  # This is a general function taking in a matrix representing \geq statements on the ballot shares (where the last column is the RHS and the others are the coefficients on the LHS) and outputting the vertices of the area we need to integrate
  a_mat <- inequality_mat[ ,-ncol(inequality_mat)]
  b_vec <- inequality_mat[ ,ncol(inequality_mat)]
  # a_mat and b_vec are stated in terms of \geq statements.
  # a_mat v \geq b_vec where v and b_vec are B-length column vectors
  # and a_mat has dimensions K \times B (K conditions).
  # to_make_near_equality is a vector of indices of conditions we want to convert to near-equalities

  # we add the 2 simplex conditions
  # first is that each vote share has to be non-negative
  positivity_part <- diag(ncol(a_mat))

  # adjust a_mat and b_vec for any near inequalities we are making
  # the \geq statement should be zero

  # makeH represents the polyhedron described by the inequalities and the simplex conditions. it is just a reformatting of the conditions you put in -- no computation.
  # a1 and b1 describe \leq inequalities, so I make LHS and RHS of my \geq inequalities (a_mat and b_vec) negative.
  # a2 and b2 describe an exact equality - here, just keeping everything on the simplex
  the_Hrep <- rcdd::makeH(a1 = -rbind(a_mat, positivity_part), b1 = c(-b_vec, rep(0, nrow(positivity_part))), a2 = matrix(1, ncol = ncol(a_mat), nrow = 1), b2 = 1)

  # gets the convex hull of points from H representation
  out <- rcdd::scdd(the_Hrep)

  out$output[, -c(1,2)] ## I am not sure what the first two columns are for.
}


S_array_from_vertices_of_integration_region <- function(voir){

  # simplified version of simplices_to_integrate_from_win_region_vertices, which subsets to the facets where a constraint binds.

  simplex_mat <- geometry::delaunayn(voir[, -ncol(voir)])

  data.frame(simplex_mat) %>%
    mutate(simplex = 1:nrow(.)) %>%
    pivot_longer(cols = starts_with("X"), values_to = "vertex", names_to = "name") %>%
    select(-name) -> simplex_vertex

  voir_df <- data.frame(vertex = 1:nrow(voir), voir)

  simplex_vertex %>%
    left_join(voir_df, by = "vertex") %>%
    select(-vertex) %>% as.matrix() -> simplices

  # in the S array we need for SimplicialCubature, the columns are the vertices of the simplices over which we integrate. (note that's not the usual way R works.)
  matrix_list <- list()
  ss <- unique(simplices[,"simplex"])
  for(s in ss){
    matrix_list[[as.character(s)]] = t(simplices[which(simplices[,"simplex"] == s), -1])
  }

  array(matrix_list %>% unlist(), dim = c(nrow(matrix_list[[1]]), ncol(matrix_list[[1]]), length(ss)))

}

plurality_inequality_matrix <- function(k, n = 5000){
  cbind(1, -diag(k - 1), rep(1/n, k-1))
}

S_array_from_inequalities <- function(inequality_mat){

  inequality_mat %>%
    vertices_of_integration_region_from_inequalities() %>%
    S_array_from_vertices_of_integration_region()

}

test <- F
if(test){
  k <- 4
  n <- 5000
  # plurality outright win conditions: v_a - v_k > 1/n forall k
  im <- plurality_inequality_matrix(4, n = n)

  voir_a <- vertices_of_integration_region_from_inequalities(im)
  S_a <- S_array_from_inequalities(im)
  S_ab <- S_array_from_inequalities(im %>% make_inequalities_into_near_equalities(1))

  voir_abc <- vertices_of_integration_region_from_inequalities(im %>% make_inequalities_into_near_equalities(c(1,2)))

  im_abc <- im %>% make_inequalities_into_near_equalities(c(1,2), n = 6000)
  voir_abc <- vertices_of_integration_region_from_inequalities(im_abc)
  S_abc <- S_array_from_inequalities(im_abc)
  out <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, S = S_abc, alpha = c(10, 8, 4, 3), maxEvals = 100000)
  out$integral
}


# so now I need to cycle through. for plurality (or any positional method)
# each combination of conditions made into a near equality is a pivot event.

# I should probably do this generically for positional methods? but hard to extend given number of ballots. so I start with plurality.

# all I need is this: how do I get a list of all combinations of indices
# c(), c(1), c(2), c(3), ..., c(1,2), c(1,2,3)
# k <- 4
# combn(k, 0)
# combn(k, 1)
# combn(k, 2)
# combn(k, 3)


exact_plurality_pivot_probabilities <- function(n = 5000, alpha = NULL, mu = NULL, sigma = NULL, precision = NULL, cand_names = NULL, sep = "", ...){

  if(!is.null(alpha) | (!is.null(mu) & !is.null(precision)) & is.null(sigma)){
    # this is Dirichlet
    distribution = "dirichlet"
    if(is.null(alpha)){alpha <- mu*precision}
    k <- length(alpha)
  }else if(!is.null(mu) & !is.null(sigma)){
    # this is logistic normal
    distribution = "logisticnormal"
    k <- length(mu)
  }else{
    stop("Please pass parameters that allow me to determine the intended distribution over election outcomes.")
  }

  if(is.null(cand_names)){cand_names <- letters[1:k]}

  im <- plurality_inequality_matrix(k, n = n)
  stopifnot(nrow(im) == k-1)

  generic_W_mat <- matrix(0, nrow = k, ncol = k, dimnames = list(cand_names, cand_names))

  out <- list()

  # we cycle over candidates
  for(i in 1:k){
    # we will shuffle the parameters and cand_names so that this candidate is first
    these_indices <- c(i,(1:k)[-i])
    these_cand_names <- cand_names[these_indices]
    # for each number of candidates who could be nearly tied with this one
    for(m in 0:(k-1)){
      # we enumerate the indices of the constraints corresponding to the possible candidates who could be in a combination of suze m
      near_equality_index_mat <- combn(k-1, m)
      for(j in 1:ncol(near_equality_index_mat)){
        # nei is a vector of those indices
        nei <- near_equality_index_mat[,j]
        # im is the set of inequalities that results when we add the near-equalities
        this_im <- im %>% make_inequalities_into_near_equalities(nei, n)
        # this_S is the space defined by this_im divided into an array of simplices
        this_S <- S_array_from_inequalities(this_im)
        # we name the pivot event: e.g. "a" means a wins outright no matter what extra ballot is submitted, "a_bc" means b and c are each less than one vote behind
        this_name <- paste0(c(these_cand_names[1], paste0(sort(these_cand_names[1 + nei]), collapse = sep)), collapse = "_")
        if(distribution == "dirichlet"){
          out[[this_name]] = SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[these_indices], ...)
        }else if(distribution == "logisticnormal"){
          out[[this_name]] = SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[these_indices], sigma = sigma[these_indices, these_indices], ...)
        }
        this_W_mat <- generic_W_mat
        for(kand in (nei + 1)){
          this_W_mat[these_cand_names[kand], these_cand_names[kand]] <- 1
        }
        this_W_mat[i, which(apply(this_W_mat, 2, sum) == 0)] <- 1
        out[[this_name]]$W_mat <- this_W_mat
      }
    }
  }
  out
}

P_mat_from_eppp <- function(out){
  # loopy way
  P <- out[[names(out)[1]]]$W_mat*out[[names(out)[1]]]$integral
  for(j in 2:length(names(out))){
    P <- P + out[[names(out)[j]]]$W_mat*out[[names(out)[j]]]$integral
  }
  P
  # probably a dplyr way. but couldn't work it out
  #tibble(name = names(out),
   #      pie = out %>% map("integral"),
    #     W = out %>% map("W_mat")) -> tpW
}

test <- F
if(test){
  out <- exact_plurality_pivot_probabilities(alpha = c(10, 9, 6, 5), tol = .001, maxEvals = 100000)
  P <- P_mat_from_eppp(out)

  # logistic normal case -- takes a bit longer?
  this_sigma <- diag(4)*c(.5, .25, .1, 0)
  out_ln <- exact_plurality_pivot_probabilities(mu = c(1,.5, .25, 0), sigma = this_sigma, tol = .001, maxEvals = 100000)
  P <- P_mat_from_eppp(out_ln)


  out <- exact_plurality_pivot_probabilities(alpha = c(10, 9, 6), tol = .001, maxEvals = 100000)
  P <- P_mat_from_eppp(out)

}



# OK that's cool!
# so what do I have: I have a way to get all the pivot probabilities given parameters. I have a justification of pivot probabilities.
# I can write this up I guess.
# then from there can consider approximations:
  # collapse adjacent pivot probs (have them straddle or not?)
  # approximate pivot probs by (1/n) times integral along the line.
