#' Compute IRV winners for simulated elections
#'
#' Takes a data frame of simulated ordinal elections (one election per row,
#' one ballot ordering per column) and returns a data fram of IRV winners.
#'
#' @param sims A tibble of simulations, with an `id` column and
#' a (named) column for each distinct ballot (e.g. ABCD, ABDC,...)
#' @return A tibble with columns `id` and `winner`
#' @examples
#' sims <- simulate_ordinal_results_from_dirichlet(k = 3, n = 10)
#' irv_winners(sims)

irv_winners <- function(sims){
  stopifnot(length(unique(colnames(sims))) == ncol(sims))
  while(!near(max(as.matrix(sims)[1,-which(names(sims) == "id")], na.rm = T), 1)){
    sims <- one_round_of_irv(sims)
  }
  sims %>%
    pivot_longer(cols = -id, names_to = "winner", values_to = "won") %>%
    filter(near(won, 1)) %>%
    select(id, winner)
}

one_round_of_irv <- function(sims){
  sims %>%
    mutate(loser = get_loser_from_tibble_of_sims(.)) %>%
    # now we group so that we only have to drop and condense once per losing candidate (rather than once per row)
    group_by(loser) %>%
    nest() %>%
    mutate(reduced = map(data, drop_candidate_and_condense, loser)) %>%
    ungroup() %>%
    select(reduced) %>%
    unnest(cols = c(reduced)) -> out

  out %>%
    colSums(na.rm = T) -> col_sums

  col_sums[col_sums > 0] %>%
    names() %>%
    sort() -> non_zero_cols

  out %>%
    select(all_of(non_zero_cols)) %>%
    arrange(id) %>%
    relocate(id)
}


