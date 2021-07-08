#' Generate simulated election results from Dirichlt distribution
#'
#' @param k The number of candidates
#' @param n The number of simulated elections to draw
#' @param alpha Optional vector of Dirichlet parameters of length `factorial(k)`.
#' If omitted, random parameters will be drawn.
#' @return A tibble with a column `id` and a column for each unique ballot type.
#' @examples
#' simulate_ordinal_results_from_dirichlet(n = 10)
#' @export
simulate_ordinal_results_from_dirichlet <- function(k = 4, n = 10000, alpha = NULL){
  stopifnot(k >= 3)
  if(is.null(alpha)){alpha <- sample(1:12, size = factorial(k), replace = T)}
  stopifnot(length(alpha) == factorial(k))
  sims <- gtools::rdirichlet(n, alpha = alpha)
  colnames(sims) <- gtools::permutations(n = k, r = k, v = LETTERS[1:k]) %>%
    apply(1, paste, collapse = "")
  sims %>%
    as_tibble() %>%
    mutate(id = 1:n()) %>%
    relocate(id)
}
