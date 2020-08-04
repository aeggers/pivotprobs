
positional_piv_probs_simulation <- function(sims, tol = .01, s = .5){

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  # this is the brute force -- nothing fancy
  list(
    "AB" = mean(score_a > score_b & score_a - score_b < tol & score_b > score_c)/tol,
    "AC" = mean(score_a > score_c & score_a - score_c < tol & score_c > score_b)/tol,
    "BC" = mean(score_b > score_c & score_b - score_c < tol & score_c > score_a)/tol
  )

}
