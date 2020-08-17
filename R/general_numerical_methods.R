## general approach

cand_a_win_region_vertices_from_win_conditions <- function(win_conditions){

  positivity_part <- -diag(ncol(win_conditions))

  a1 <- rbind(-win_conditions, positivity_part)

  cand_a_win_region_H <- makeH(a1 = a1, b1 = rep(0, nrow(a1)), a2 = rep(1, ncol(win_conditions)), b2 = 1)

  out <- rcdd::scdd(cand_a_win_region_H)

  out$output[, -c(1,2)]
}

# test:
wrv3 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0), c(1, 0, -1)))
wrv4 <- cand_a_win_region_vertices_from_win_conditions(rbind(c(1, -1, 0,0), c(1, 0, -1,0), c(1, 0, 0, -1)))

simplices_to_integrate_from_win_region_vertices <- function(wrv, a_beats_b_condition = c(1,-1,0,0)){

  # which vertices in wrv have an ab tie?
  # TODO: worry about floating point here
  vertices_with_ab_tie <- which(as.vector(wrv %*% matrix(a_beats_b_condition, ncol = 1)) == 0)

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
    select(-vertex) %>%
    nest() -> nested_simplices

  # bit messy here
  ns <- nested_simplices$data
  array(ns %>% unlist(), dim = c(nrow(ns[[1]]), ncol(ns[[1]]), length(ns)))
}

# test
simplices_to_integrate_from_win_region_vertices(wrv4, a_beats_b_condition = c(1,-1,0,0))

# a little unclear what we want here:
#TODO: work this out a bit. will it always be the case that we can re-label the candidates and get the right answer by focusing on a vs b? with Dirichlet, this is just reshuffling the alpha vector. If it were e.g. normal distribution, we can transform the mean and variance-covariance matrix. so I think it's okay to focus on this approach.

# so next we put things together:
# supply the conditions for a to win (with the condition we want to satisfy with equality being the first condition); get the S array from that (using above); then iterate through Dirichlet parameters to integrate this for each candidate. (generalize to other belief distributions later.)
# SimplicialCubature::adaptIntegrateSimplex
