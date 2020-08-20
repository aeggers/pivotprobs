## condorcet

# a function to compute the probability of a decisive Condorcet tie between each pair of candidates
condorcet_pivot_probs_from_sims <- function(sims, tol = .01, cand_names = NULL, sep = "", kemeny = T){

  if(is.null(cand_names)){
    cand_names <- letters[1:3]
  }

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  a_vs_b <- apply(sims[,c(1,2,5)], 1, sum)
  a_vs_c <- apply(sims[,c(1,2,3)], 1, sum)
  b_vs_c <- apply(sims[,c(1,3,4)], 1, sum)

  normalizer <- tol/(sqrt(6)/2) # a little unsure about this.

  out <- list(
    mean(a_vs_c > 1/2 & b_vs_c > 1/2 & abs(a_vs_b - 1/2) < tol/2)/normalizer,
    mean(a_vs_b > 1/2 & b_vs_c < 1/2 & abs(a_vs_c - 1/2) < tol/2)/normalizer,
    mean(a_vs_b < 1/2 & a_vs_c < 1/2 & abs(b_vs_c - 1/2) < tol/2)/normalizer
    )

  out_names <- c(paste0(cand_names[1], sep, cand_names[2]), paste0(cand_names[1], sep, cand_names[3]), paste0(cand_names[2], sep, cand_names[3]))

  names(out) <- out_names

  if(kemeny){
    forward_cycle <- a_vs_b > .5 & a_vs_c < .5 & b_vs_c > .5
    reverse_cycle <- a_vs_b < .5 & a_vs_c > .5 & b_vs_c < .5

    out[[paste0(out_names[1], "_forward")]] <- mean(
      # forward cycle, and
        forward_cycle &
          # c's loss to b worse than b's loss to a
          b_vs_c > a_vs_b &
          # b's loss to a worse than a's loss to c
          a_vs_b > (1 - a_vs_c) &
          # but just barely
          a_vs_b - (1 - a_vs_c) < tol)/normalizer

    out[[paste0(out_names[1], "_reverse")]] <- mean(
        # reverse cycle
          reverse_cycle &
            # c's loss to a worse than a's loss to b
            a_vs_c > (1 - a_vs_b) &
            # b's loss to c worse than a's loss to b
            (1 - b_vs_c) > (1 - a_vs_b) &
            # but just barely
            (1 - b_vs_c) - (1 - a_vs_b) < tol)/normalizer

    out[[paste0(out_names[2], "_forward")]] <- mean(
      # forward cycle, and
        forward_cycle &
          # b's loss to a worse than c's loss to b
          a_vs_b > b_vs_c &
          # c's loss to b worse than a's loss to c
          b_vs_c > (1 - a_vs_c) &
          # but just barely
          b_vs_c - (1 - a_vs_c) < tol)/normalizer

    out[[paste0(out_names[2], "_reverse")]] <- mean(
        # reverse cycle
          reverse_cycle &
            # b's loss to c worse than a's loss to b
            (1 - b_vs_c) > (1 - a_vs_b) &
            # a's loss to b worse than c's loss to a
            (1 - a_vs_b) > a_vs_c &
            # but just barely
            (1 - a_vs_b) - a_vs_c < tol)/normalizer

    out[[paste0(out_names[3], "_forward")]] <- mean(
      # forward cycle, and
        forward_cycle &
          # a's loss to c worse than b's loss to a
          a_vs_b > b_vs_c &
          # b's loss to a worse than c's loss to b
          b_vs_c > (1 - a_vs_c) &
          # but just barely
          b_vs_c - (1 - a_vs_c) < tol)/normalizer

    out[[paste0(out_names[3], "_reverse")]] <- mean(
        # OR reverse cycle
          reverse_cycle &
            # a's loss to b worse than c's loss to a
            (1 - a_vs_b) > a_vs_c &
            # c's loss to a worse than b's loss to c
            a_vs_c > (1 - b_vs_c) &
            # but just barely
            a_vs_c - (1 - b_vs_c) < tol)/normalizer

  }

  out

}
