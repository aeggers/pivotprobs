
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
    tidyr::pivot_longer(cols = starts_with("X"), values_to = "vertex", names_to = "name") %>%
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

