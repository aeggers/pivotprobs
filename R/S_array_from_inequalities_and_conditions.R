#' Converts system of inequalities into an array of simplices for integration
#'
#' The \code{SimplicialCubature::adaptIntegrateSimplex()} function requires an array
#' of vertices \code{S} over which to integrate. This function produces
#' this array from
#' a matrix of inequalities and optional arguments specifiying alterations to
#' these inequalities.
#'
#' For estimating pivot event probabilities, the general procedure is to specify
#' a matrix of conditions under which a given candidates wins
#' (\code{inequality_mat}) and
#' then specify which of these conditions (\code{rows_to_alter}) should
#' be met with near equality (\code{drop_dimension=F}) or
#' exact equality (\code{drop_dimension=T}). The function could be used for
#' other purposes, however.
#'
#' @param inequality_mat An R x C+1 matrix. Each row is an inequality of the
#' form Ax > b. The coefficients A come from the first C columns,
#' the constant b comes from the last column, and the x's are (in this case)
#' the shares of each ballot. For example, if there were only two ballots
#' and the sole condition is \code{x_1 - x_2 > .1}, \code{inequality_mat}
#' would be \code{matrix(c(1, -1, .1), nrow = 1)}.
#' @param rows_to_alter The indices of any rows in \code{inequality_mat} that you
#' wish to set to equality or near-equality. To generate the \code{S} array for
#' the plurality pivot event in which candidate 1 ties candidate 2,
#' specify the row in  \code{inequality_mat} that indicates that candidate 1
#' must win more than candidate 2.
#' @param drop_dimension If \code{F} and \code{rows_to_alter} are supplied,
#' then the inequality in \code{inequality_mat} is converted to a near inequality
#' based on the \code{limits} and the result represents a narrow hypervolume
#' over which to integrate. If \code{T} and \code{rows_to_alter} are supplied,
#' then the inequality in \code{inequality_mat} is converted to a near inequality
#' at \code{b=limits[1]} and the result represents a hyperplane over which to
#' integrate.
#' @param limits Determine the range within which the inequality
#' in \code{rows_to_alter} is required to hold. If the inequality
#' to be altered was \code{v_1 - v_2 > 0} and \code{drop_dimension=F}, the new conditions
#' are \code{v_1 - v_2 > limits[1]} and \code{v_1 - v_2 < limits[2]}.
#' If \code{drop_dimension=T}, the new condition is effectively
#' \code{v_1 - v_2 = limits[1]}.
#' @param epsilon Tolerance for determining which vertices of the convex hull
#' defined by \code{inequality_mat} satisfy the altered conditions when
#' \code{drop_dimension} is specified.
#' @param qhull_options These are passed to \code{geometry::delaunayn()}
#' and \code{geometry::convexhulln()}. See QHull documentation.
#'
#' @return An array of matrices, each representing one simplex. Note that
#' each *column* of each matrix is one vertex. (You might expect rows to
#' indicate vertices.)
#'
#' @examples
#' # plurality conditions with 3 candidates
#' n <- 1000
#' inequality_mat <- rbind(c(1,-1,0,1/n),
#'                         c(1,0,-1,1/n))
#' # simplices covering region where candidate 1 wins
#' S_array_from_inequalities_and_conditions(inequality_mat)
#' # simplices covering region where candidate 1 finishes just ahead of candidate 2
#' S_array_from_inequalities_and_conditions(inequality_mat,
#'        rows_to_alter = c(1), limits = c(0, 1/n))
#' # simplices covering facet where candidates 1 and 2 tie
#' S_array_from_inequalities_and_conditions(inequality_mat,
#'        rows_to_alter = c(1), drop_dimension = T, limits = c(0, 1/n))
#'
#'
#' @export
S_array_from_inequalities_and_conditions <- function(inequality_mat, rows_to_alter = c(NULL), drop_dimension = F, limits = c(0,NULL), epsilon = 1.0e-10, qhull_options = NULL){

  if(drop_dimension & length(rows_to_alter) > 1){
    warning("Results are unreliable when dropping a dimension on a compound pivot event (i.e. when multiple rows_to_alter).")
  }

  if(length(rows_to_alter) >= 1){
    # the inequality mat specifies conditions where the guy wins outright, so the rows_to_alter constant is typically 1/n.
    # when we alter conditions we go from limits[1] to limits[2]
    inequality_mat[rows_to_alter, ncol(inequality_mat)] <- limits[1]
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
  # it would make sense to apply an equality condition when drop_dimension=T,
  # but this doesn't work with the tesselation.

  if(!drop_dimension & length(rows_to_alter) >= 1){ #i.e. if this is a pivot event and we are not dropping a dimension
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
    a1_mat <- rbind(a1_mat, -a1_mat[rows_to_alter, ]) # add some \leq conditions
    b1_vec <- c(b1_vec, rep(limits[2], length(rows_to_alter)))
  }

  # get the H representation -- this is just a reformatting
  the_Hrep <- rcdd::makeH(a1 = a1_mat, b1 = b1_vec, a2 = a2_mat, b2 = b2_vec)

  # if(save_it){save(the_Hrep, file = paste0("~/Dropbox/research/strategic_voting/pivotprobs_paper/data/the_Hrep.RData"))}

  # rational arithmetic! very fun.
  the_Vrep <- rcdd::d2q(the_Hrep) %>% rcdd::scdd()
  # get the vertices of the convex hull from H representation
  all_vertices <- the_Vrep$output[,-c(1,2)] %>% rcdd::q2d()

  # if user passes a set of conditions that cannot be met
  if(nrow(all_vertices) == 0){
    warning("The supplied conditions are not met at any vertices.")
    return (NULL)
  }

  # if(save_it){save(all_vertices, file = paste0("~/Dropbox/research/strategic_voting/pivotprobs_paper/data/all_vertices.RData"))}

  # get the tesselated/triangulated convex hull (tch)
  # this is a matrix with one row per triangle
  # if a row is c(3, 1, 2), it means rows 3, 1, and 2 of the input make one of the triangles.
  if(!(drop_dimension & length(rows_to_alter) >= 1)){
    # i.e. if we are not dropping a dimension, or this is not a pivot event
    # we are integrating over a D-1 dimensional space, e.g. for plurality with 3 candidates we have 2-dimensional areas to integrate over. we triangulate and return an array of these triangles.
    N <- ncol(all_vertices) - 1
    QzQx <- ifelse(N >= 4, "Qx", "Qz")
    tch <- geometry::delaunayn(all_vertices[,-ncol(all_vertices)], options = paste("Qcc Qc Qt", QzQx, qhull_options, sep = " ")) #
    vertices_satisfying_conditions <- 1:nrow(all_vertices)
  }else{
    # we are integrating on at most D-2 dimensional facets, e.g. for plurality with 3 candidates we have a line or even a point. we triangulate the convex hull, look for facets that meet the conditions, and return an array of those triangles.
    tch <- geometry::convhulln(all_vertices[,-ncol(all_vertices)], options = paste("Tv", qhull_options, sep = " "))
    # check which vertices satisfy the conditions with equality
    # all_vertices is V by B, relevant part of inequality_mat is rows_to_alter by V+1
    condition_matrix <- inequality_mat[rows_to_alter,-ncol(inequality_mat), drop = F]
    #if(!is.matrix(condition_matrix)){
    #  condition_matrix <- matrix(condition_matrix, nrow = length(rows_to_alter), ncol = ncol(inequality_mat) - 1)
    #}
    # we use matrix multiplication to check the conditions --
    deviations <- all_vertices %*% t(condition_matrix) - matrix(inequality_mat[rows_to_alter,ncol(inequality_mat)], nrow = nrow(all_vertices), ncol = length(rows_to_alter), byrow = T) # one column per condition that is supposed to be satisfied with equality. I think I should have used drop=F here.
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





