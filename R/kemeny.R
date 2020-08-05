##

kemeny_piv_probs_simulation <- function(sims, tol = .01){

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  a_over_b <- sims[,1] + sims[,2] + sims[,5]
  a_over_c <- sims[,1] + sims[,2] + sims[,3]
  b_over_c <- sims[,3] + sims[,4] + sims[,1]

  forward_cycle <- a_over_b > .5 & a_over_c < .5 & b_over_c > .5
  reverse_cycle <- a_over_b < .5 & a_over_c > .5 & b_over_c < .5

  list(
    # outright condorcet winner probs
    # e.g. A and B both beat C and A barely beats B
    "AB" = mean(a_over_b > .5 & a_over_c > .5 & b_over_c > .5 & a_over_b - .5 < tol)/tol,
    "AC" = mean(a_over_b > .5 & a_over_c > .5 & b_over_c < .5 & a_over_c - .5 < tol)/tol,
    "BC" = mean(a_over_b < .5 & a_over_c < .5 & b_over_c > .5 & b_over_c - .5 < tol)/tol,
    # cycle probs
    # e.g. forward cycle (A beats B, B beats C, C beats A); C's loss is the worst, and A's loss to C is slightly better than B's loss to A
    "AB_cycle" = mean(
      (# forward cycle, and
        forward_cycle &
          # c's loss to b worse than b's loss to a
          b_over_c > a_over_b &
          # b's loss to a worse than a's loss to c
          a_over_b > (1 - a_over_c) &
          # but just barely
          a_over_b - (1 - a_over_c) < tol) |
        (  # OR reverse cycle
          reverse_cycle &
            # c's loss to a worse than a's loss to b
            a_over_c > (1 - a_over_b) &
            # b's loss to c worse than a's loss to b
            (1 - b_over_c) > (1 - a_over_b) &
            # but just barely
            (1 - b_over_c) - (1 - a_over_b) < tol))/tol,
    "AC_cycle" = mean(
      (# forward cycle, and
        forward_cycle &
          # b's loss to a worse than c's loss to b
          a_over_b > b_over_c &
          # c's loss to b worse than a's loss to c
          b_over_c > (1 - a_over_c) &
          # but just barely
          b_over_c - (1 - a_over_c) < tol) |
        (  # OR reverse cycle
          reverse_cycle &
            # b's loss to c worse than a's loss to b
            (1 - b_over_c) > (1 - a_over_b) &
            # a's loss to b worse than c's loss to a
            (1 - a_over_b) > a_over_c &
            # but just barely
            (1 - a_over_b) - a_over_c < tol))/tol,
    "BC_cycle" = mean(
      (# forward cycle, and
        forward_cycle &
          # a's loss to c worse than b's loss to a
          a_over_b > b_over_c &
          # b's loss to a worse than c's loss to b
          b_over_c > (1 - a_over_c) &
          # but just barely
          b_over_c - (1 - a_over_c) < tol) |
        (  # OR reverse cycle
          reverse_cycle &
            # a's loss to b worse than c's loss to a
            (1 - a_over_b) > a_over_c &
            # c's loss to a worse than b's loss to c
            a_over_c > (1 - b_over_c) &
            # but just barely
            a_over_c - (1 - b_over_c) < tol))/tol
  )

}
