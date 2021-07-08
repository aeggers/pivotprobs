#' Compute pivot probabilities and generate P-matrices for four-candidate IRV elections
#'
#' @param sims A tibble of simulations, with an `id` column and
#' a column for each distinct ballot
#' @param n The number of voters.
#' @param reporting The level of reporting detail (0: none; 1:
#' type of pivot event; 2: name of pivot event)
#' @return A list of pivot events, each of which is a list with elements
#' `integral` and `P`
#' @examples
#' sims <- simulate_ordinal_results_from_dirichlet(k = 4, n = 10000)
#' out <- irv_pivot_probs_four_cands(sims)
#' out %>% combine_P_matrices()
#' @export
irv_pivot_probs_four_cands <- function(sims, n = 1000, reporting = 1) {
  container <- new.env()

  container$PP_LIBRARY <- list() ## warning: writes to global variable
  if (reporting >= 1) {
    cat("Round 0: ")
  }
  container$round_0 <- sims %>% round_0_pivot_probs(reporting = reporting)
  if (reporting >= 1) {
    cat("done.\nRound 1: ")
  }
  container$round_1 <- sims %>% round_1_pivot_probs(reporting = reporting)
  if (reporting >= 1) {
    cat("done.\nRound 2: ")
  }
  container$round_2 <- sims %>% last_round_pivot_probs(reporting = reporting)
  if (reporting >= 1) {
    cat("done.\n")
  }
  c(container$round_0, container$round_1, container$round_2)
}
