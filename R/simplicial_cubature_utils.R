# simplicial cubature utils

# This is the key function -- I think it is self-contained but should check
S_array_from_inequalities_and_conditions <- function(inequality_mat, rows_to_alter = c(NULL), drop_dimension = F, limits = c(0,NULL), epsilon = 1.0e-10, qhull_options = NULL){

  if(drop_dimension & length(rows_to_alter) > 1){
    warning("Results are unreliable when dropping a dimension on a compound pivot event (i.e. when multiple rows_to_alter).")
  }

  # reformat the inequality matrix  and add the simplex constraints ####
  a1_mat <- inequality_mat[ ,-ncol(inequality_mat)]
  b1_vec <- inequality_mat[ ,ncol(inequality_mat)]
  # the simplex inequality constraints
  positivity_part <- diag(ncol(a1_mat))
  a1_mat <- -rbind(a1_mat, positivity_part) # and make negative, as we provide \geq statements but require \leq statements
  b1_vec <- -c(b1_vec, rep(0, nrow(positivity_part)))

  # the simplex equality condition
  a2_mat <- matrix(1, nrow = 1, ncol = ncol(a1_mat))
  b2_vec <- c(1)

  if(!drop_dimension & length(rows_to_alter) >= 1){
    if(is.null(limits[2])){
      stop("If not dropping a dimension, you must provide upper and lower limits, e.g. 0 and 1/n.")
    }
    if(limits[1] == limits[2]){
      stop("If not dropping a dimension, the upper and lower limits must not be equal.")
    }
    # we are going to convert each specified inequality condition into two near-equality conditions to capture a near tie.
    # the current conditions are a1[rta,] v \geq b1[rta]
    # we need to convert these to a1[rta,] v \geq limits[1]
    # and add one that is a1[rta,] v \leq b1[rta]
    b1_vec[rows_to_alter] <- -limits[1] # recycling
    a1_mat <- rbind(a1_mat, -a1_mat[rows_to_alter, ]) # add some \leq conditions
    b1_vec <- c(b1_vec, rep(limits[2], length(rows_to_alter)))
  }

  # get the H representation -- this is just a reformatting
  the_Hrep <- rcdd::makeH(a1 = a1_mat, b1 = b1_vec, a2 = a2_mat, b2 = b2_vec)

  # get the vertices of the convex hull from H representation
  all_vertices <- rcdd::scdd(the_Hrep)$output[,-c(1,2)]

  # if user passes a set of conditions that cannot be met
  if(nrow(all_vertices) == 0){
    warning("The supplied conditions are not met at any vertices.")
    return (NULL)
  }

  # get the tesselated/triangulated convex hull (tch)
  # this is a matrix with one row per triangle
  # if a row is c(3, 1, 2), it means rows 3, 1, and 2 of the input make one of the triangles.
  # not the same as geometry::convhulln(), which I had used before, though I don't see the difference
  if(!(drop_dimension & length(rows_to_alter) >= 1)){
    # we are integrating over a D-1 dimensional space, e.g. for plurality with 3 candidates we have 2-dimensional areas to integrate over. we triangulate and return an array of these triangles.
    N <- ncol(all_vertices) - 1
    QzQx <- ifelse(N >= 4, "Qx", "Qz")
    tch <- geometry::delaunayn(all_vertices[,-ncol(all_vertices)], options = paste("Qcc Qc Qt", QzQx, qhull_options, sep = " ")) #
    vertices_satisfying_conditions <- 1:nrow(all_vertices)
  }else{
    # we are integrating on at most D-2 dimensional facets, e.g. for plurality with 3 candidates we have a line or even a point. we triangulate the convex hull, look for facets that meet the conditions, and return an array of those triangles.
    tch <- geometry::convhulln(all_vertices[,-ncol(all_vertices)], options = paste("Tv", qhull_options, sep = " "))
    # check which vertices satisfy the conditions
    # all_vertices is V by B, relevant part of inequality_mat is rows_to_alter by V+1
    condition_matrix <- inequality_mat[rows_to_alter,-ncol(inequality_mat)]
    if(!is.matrix(condition_matrix)){
      condition_matrix <- matrix(condition_matrix, nrow = length(rows_to_alter), ncol = ncol(inequality_mat) - 1)
    }
    deviations <- all_vertices %*% t(condition_matrix) - matrix(inequality_mat[rows_to_alter,ncol(inequality_mat)], nrow = nrow(all_vertices), ncol = length(rows_to_alter), byrow = T) # one column per condition that is supposed to be satisfied with equality
    vertices_satisfying_conditions <- which(apply(abs(deviations) < epsilon, 1, all))

    if(length(vertices_satisfying_conditions) == 0){
      # the conditions don't apply at of these vertices
      stop("The conditions supplied are not met at any vertex of the convex hull of the inequality matrix supplied.")
    }

  }

  #### locate the relevant simplices and organize in an S array ####

  # pivot longer, making the simplex the group
  data.frame(tch) %>%
    mutate(simplex = 1:nrow(.)) %>%
    pivot_longer(cols = starts_with("X"), values_to = "vertex", names_to = "name") %>%
    select(-name) %>%
    group_by(simplex) %>%
    # and select simplices satisfying conditions (if relevant)
    filter(sum(!vertex %in% vertices_satisfying_conditions) == 0) -> simplex_vertex

  # if we are dropping a dimension on a compound event, sometimes we don't get any simplices
  if(nrow(simplex_vertex) == 0){
  # This happens when the conditions don't apply on any facet of the convex hull.
  # e.g. if k=3 a three-way tie happens at a point; at k = 4 a three-way tie happens on a line; at k = 5  it happens on a space spanned by 4 points.
  # when k=3 or k=4 we can just return the result as a matrix and it works. but for k=5, integrating on this gives us 0. so in the main code we do not compute compound pivot event probabilities when we drop a dimension.
    this <- all_vertices[vertices_satisfying_conditions,]
    if(!is.matrix(this)){this <- matrix(this, nrow = 1)}
    return(t(this))
  }

  # otherwise we need to organize into an S array -- a bit tedious and messy code.
  all_vertices_df <- data.frame(vertex = 1:nrow(all_vertices), all_vertices)

  simplex_vertex %>%
    left_join(all_vertices_df, by = "vertex") %>%
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
  n <- 1000
  win_3 <- rbind(c(1,-1,0,1/n), c(1,0,-1, 1/n))
  S_array_from_inequalities_and_conditions(win_3)
  S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(0, 1/n))
  S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(2), limits = c(0, 1/n))
  S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1, 2), limits = c(0, 1/n))
  # this is straddling (i.e. merging pivot events) but not dropping a dimension
  S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(-1/(2*n), 1/(2*n)))
  # dropping a dimension -- just a line (though not centered, so not merging pivot probs)
  S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(0, 1/n), drop_dimension = T)
  # dropping a dimension and centered (merging pivot probs)
  win_3a <- rbind(c(1,-1,0,0), c(1,0,-1, 0))
  S_array_from_inequalities_and_conditions(win_3a, rows_to_alter = c(1), limits = c(0, 1/n), drop_dimension = T)
  S_array_from_inequalities_and_conditions(win_3a, rows_to_alter = c(2), limits = c(0, 1/n), drop_dimension = T)
  # a point.
  S_array_from_inequalities_and_conditions(win_3a, rows_to_alter = c(1,2), limits = c(0, 1/n), drop_dimension = T)
  win_4a <- rbind(c(1,-1,0,0,0), c(1,0,-1, 0,0), c(1,0,0, -1,0))
  S_array_from_inequalities_and_conditions(win_4a, limits = c(0, 1/n))
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1), limits = c(0, 1/n), drop_dimension = F)
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1), limits = c(-1/n, 1/n), drop_dimension = F)
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1,2), limits = c(0, 1/n), drop_dimension = F)
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1,2), limits = c(-1/n, 1/n), drop_dimension = F)
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1), limits = c(0, 1/n), drop_dimension = T)
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1,2), limits = c(0, 1/n), drop_dimension = T)
  # a point
  S_array_from_inequalities_and_conditions(win_4a, rows_to_alter = c(1,2,3), limits = c(0, 1/n), drop_dimension = T)

  n <- 5000
  win_3 <- rbind(c(1,-1,0,1/n), c(1,0,-1, 1/n))
  alpha3 <- c(10, 6, 4)
  this_S_1 <- S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(0, 1/n), drop_dimension = F)
  out_1 <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, this_S_1, alpha = alpha3)
  this_S_2 <- S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(-1/(2*n), 1/(2*n)), drop_dimension = F)
  out_2 <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, this_S_2, alpha = alpha3)
  this_S_3 <- S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(1/(2*n), 1/(2*n)), drop_dimension = T)
  out_3 <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, this_S_3, alpha = alpha3)
  this_S_4 <- S_array_from_inequalities_and_conditions(win_3, rows_to_alter = c(1), limits = c(0, 1/(2*n)), drop_dimension = T)
  out_4 <- SimplicialCubature::adaptIntegrateSimplex(dirichlet_for_integration, this_S_4, alpha = alpha3)
  c(out_1$integral, out_2$integral, out_3$integral/(sqrt(2)*n), out_4$integral/(sqrt(2)*n))
}




# this also used
P_mat_from_eppp <- function(out){
  # thought about a dplyr way but couldn't work it out.
  # here is a loopy way
  P <- out[[names(out)[1]]]$P*out[[names(out)[1]]]$integral
  for(j in 2:length(names(out))){
    the_integral <- out[[names(out)[j]]]$integral
    if(is.null(the_integral)){next} # "total" for example doesn't have an integral (it has "seconds_elapsed")
    if(is.na(the_integral)){next}
    P <- P + out[[names(out)[j]]]$P*the_integral
  }
  P
}





## rest is obsolete


plurality_win_conditions <- function(k, n = 5000){
  cbind(1, -diag(k - 1), rep(1/n, k-1))
}

# and this is implementing it in plurality
plurality_event_probabilities <- function(skip_non_pivot_events = F, merge_adjacent_pivot_events = F, drop_dimension = F, n = 5000, alpha = NULL, mu = NULL, sigma = NULL, precision = NULL, cand_names = NULL, sep = "", store_time = T, ...){

  if(merge_adjacent_pivot_events & drop_dimension){
    limits <- c(1/(2*n), 1/(2*n))
  } else if(merge_adjacent_pivot_events){
    limits <- c(-1/(2*n), 1/(2*n))
  } else if(drop_dimension){
    limits <- c(0, 0)
  } else{
    limits <- c(0, 1/n)
  }

  # it is plurality specific to get the number of candidates from alpha
  time_start <- Sys.time()
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

  # make matrix of inequalities for candidate a winning outright, given number of candidates k and electorate size n
  im <- plurality_win_conditions(k, n = n) # plurality specific
  stopifnot(nrow(im) == k-1)

  # W_mat indicates which candidate (rows) wins given this event and an additional ballot (columns)
  # we fill it in below
  # this part is plurality specific
  generic_W_mat <- matrix(0, nrow = k, ncol = k, dimnames = list(cand_names, cand_names))

  # storage -- one entry per election event
  out <- list()

  # we cycle over candidates
  for(i in 1:k){
    # we will shuffle the parameters and cand_names so that this candidate is first, i.e. candidate a
    these_indices <- c(i,(1:k)[-i]) # in plurality we can use this both for candidates and for alpha/mu/sigma.
    these_cand_names <- cand_names[these_indices]
    # for each number of candidates who could be nearly tied with this one
    for(m in 0:(k-1)){
      # you can skip the non-pivot events
      if(skip_non_pivot_events & m == 0){next}
      # we enumerate the indices of the constraints corresponding to the possible candidates who could be in a combination of suze m
      matrix_of_rows_to_alter <- combn(k-1, m)
      for(j in 1:ncol(matrix_of_rows_to_alter)){
        # rta is a vector of those indices, saying which constraint will bind
        rta <- matrix_of_rows_to_alter[,j]

        this_time_start <- Sys.time()
        # we name the pivot event: e.g. "a" means a wins outright no matter what extra ballot is submitted, "a_bc" means b and c are each less than one vote behind
        this_name <- paste0(c(these_cand_names[1], paste0(sort(these_cand_names[1 + rta]), collapse = sep)), collapse = "_")
        if(merge_adjacent_pivot_events){
          # first we check if this pivot event is already covered

          # strip this name
          altered_sep <- sep
          if(sep == ""){altered_sep = "STRING NOONE WOULD USE"}
          this_stripped_name <- this_name %>% str_replace_all("_", "") %>% str_replace_all(altered_sep, "") %>% str_split("") %>% `[[`(1) %>% sort() %>% paste0(collapse = "")
          stored_stripped_names <- names(out) %>% str_replace_all("_", "") %>% str_replace_all(altered_sep, "") %>% str_split("") %>% map(sort) %>% map(str_c, collapse = "") %>% unlist()
          to_grab <- which(stored_stripped_names == this_stripped_name)[1]
          if(!is.na(to_grab)){
            out[[this_name]] <- out[[names(out)[to_grab]]]
            time_diff <- Sys.time() - this_time_start
            units(time_diff) <- "secs"
            out[[this_name]]$seconds_elapsed <- as.double(time_diff)
            next
          }
        }
        # this is the crucial function: getting the S array (the array of simplices over which top integrate) from conditions, with options
        this_S <- S_array_from_inequalities_and_conditions(im, rows_to_alter = rta, limits = limits, drop_dimension = drop_dimension)
        if(distribution == "dirichlet"){
          if(class(this_S) == "matrix" && ncol(this_S) == 1){
            out[[this_name]] <- list(integral = gtools::ddirichlet(as.vector(this_S), alpha[these_indices])/sqrt(length(alpha)), functionEvaluations = 1, message = "OK")
          }else{
            out[[this_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[these_indices], ...)
          }
        }else if(distribution == "logisticnormal"){
          if(class(this_S) == "matrix" && ncol(this_S) == 1){
            out[[this_name]] <- list(integral = dlogisticnormal(as.vector(this_S), mu[these_indices], sigma[these_indices, these_indices])/sqrt(length(mu)), functionEvaluations = 1, message = "OK")
          }else{
            out[[this_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[these_indices], sigma = sigma[these_indices, these_indices], ...)
          }
        }
        # if we are integrating at equality rather than near-equality (i.e. integrating on facets rather than in thin subspaces), then we need to thicken out the space. if there were m equalities we multiply by (1/n)^m.
        if(drop_dimension){
          out[[this_name]]$integral <- out[[this_name]]$integral*((1/(sqrt(2)*n))^m)
        }
        # this part is plurality specific
        this_W_mat <- generic_W_mat
        for(kand in (rta + 1)){
          this_W_mat[these_cand_names[kand], these_cand_names[kand]] <- 1
        }
        this_W_mat[i, which(apply(this_W_mat, 2, sum) == 0)] <- 1
        out[[this_name]]$W_mat <- this_W_mat

        if(store_time){
          time_diff <- Sys.time() - this_time_start
          units(time_diff) <- "secs"
          out[[this_name]]$seconds_elapsed <- as.double(time_diff)
        }

      }
    }
  }
  if(store_time){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    out[["total"]]$seconds_elapsed <- as.double(time_diff)
  }
  out
}









plurality_event_probabilities_original <- function(n = 5000, alpha = NULL, mu = NULL, sigma = NULL, precision = NULL, cand_names = NULL, sep = "", store_time = T, ...){

  time_start <- Sys.time()
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

  # make matrix of inequalities for candidate a winning outright, given number of candidates k and electorate size n
  im <- plurality_inequality_matrix(k, n = n)
  stopifnot(nrow(im) == k-1)

  # W_mat indicates which candidate (rows) wins given this event and an additional ballot (columns)
  generic_W_mat <- matrix(0, nrow = k, ncol = k, dimnames = list(cand_names, cand_names))

  # storage -- one entry per election event
  out <- list()

  # we cycle over candidates
  for(i in 1:k){
    # we will shuffle the parameters and cand_names so that this candidate is first, i.e. candidate a
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
  if(store_time){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    out$seconds_elapsed <- as.double(time_diff)
  }
  out
}




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

  # if voir is just two points in R3, we just return the transpose of it:
  if(nrow(voir) < ncol(voir)){
    return(t(voir))
  }

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






vertices_of_integration_region_from_inequalities_and_equalities <- function(inequality_mat, indices_to_alter = c(NULL), drop_dimension = F, limits = c(0,0)){
  # this was an attempt to get the S array from a matrix of inequalities in a way that allows us to specify inequalities to alter.
  # we can specify inequalities to turn into equalities (drop_dimension). then we get vertices of the integration region.
  # this works fine if we don't drop_dimension. But if we do, then we can't tesselate the resulting vertices. delaunayn() and convhulln() can't figure out that if we have a plane in 3d, we should switch to finding triangles rather than tetrahedrons. (It only works on convex surfaces.) Seems like a simple issue, but I couldn't find a solution, and I don't think there is a simple general fix: if k = 3 we're talking about a line so no triangulation is needed, if k = 4 we have a plane and unless there are only three vertices we need to apply some rules.
  # so we'll go back to the previous approach.
  # the problem is that

  # indices_to_alter indicates the rows of inequality mat that we want to convert to a near inequality (drop_dimension = F) or an equality (drop_dimension = T).
  # in the case where we do not drop a dimension, `limits` specifies that v_a - v_b \geq limits[1] and v_a - v_b \leq limits[2]
  # in the case where we *do* drop a dimension, we convert each inequality v_a - v_b \geq 1/n into v_a - v_b = limits[1]
  # limits could in principle be a matrix when we have several indices to alter, but I can't think of why we would want to use different limits for different conditions

  # the inequalities (leaving out for now the simplex conditions, i.e. v_1 \geq 0, v_2 \geq 0, etc
  a1_mat <- inequality_mat[ ,-ncol(inequality_mat)]
  b1_vec <- inequality_mat[ ,ncol(inequality_mat)]

  # the equalities (initially just the simplex condition: \sum{v} = 1)
  a2_mat <- matrix(1, nrow = 1, ncol = ncol(a1_mat))
  b2_vec <- c(1)

  if(!drop_dimension){
    # we are going to convert each specified inequality condition into two near-equality conditions to capture a near tie. .
    # the current condition is a1[ita,] v \geq b1[ita]
    # we need to convert that one to a1[ita,] v \geq 0
    # and add one that is a1[ita,] v \leq b1[ita]
    b1_vec[indices_to_alter] <- limits[1] # recycling
    a1_mat <- rbind(a1_mat, -a1_mat[indices_to_alter, ]) # add some \leq conditions
    b1_vec <- c(b1_vec, rep(-limits[2], length(indices_to_alter)))
  }else if(drop_dimension){
    if(length(indices_to_alter) == 0){
      stop("You can't drop a dimension unless you specify indices_to_alter.")
    }
    # we are going to convert each specified inequality condition into an equality condition
    # stick the current inequalities into the equalities matrix
    a2_mat <- rbind(a2_mat, a1_mat[indices_to_alter,])
    b2_vec <- c(b2_vec, rep(limits[1], length(indices_to_alter)))
    # take them out of the inequalities matrix
    a1_mat <- matrix(a1_mat[-indices_to_alter,], ncol = ncol(a1_mat))
    b1_vec <- b1_vec[-indices_to_alter]
  }

  # add the simplex inequalities, and convert \geq statements to \leq statements
  positivity_part <- diag(ncol(a1_mat))
  a1_mat <- -rbind(a1_mat, positivity_part)
  b1_vec <- -c(b1_vec, rep(0, nrow(positivity_part)))

  # get the H representation -- this is just a reformatting
  the_Hrep <- rcdd::makeH(a1 = a1_mat, b1 = b1_vec, a2 = a2_mat, b2 = b2_vec)

  # get the vertices of the convex hull from H representation
  out <- rcdd::scdd(the_Hrep)

  matrix(out$output[, -c(1,2)], nrow = nrow(out$output)) ## I don't need the first two columns
  # matrix part is because it gets turned into a vector if just one row

  # this is my function
  # S_array_from_vertices_of_integration_region(voir)

}

S_array_from_inequalities2 <- function(inequality_mat){

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
