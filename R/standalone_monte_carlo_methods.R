# standalone monte carlo methods

#' @export
plurality_pivot_probs_from_sims <- function(sims = NULL, n = 1000, tol = .01, cand_names = NULL, sep = ""){

  out <- list()
  if(is.null(cand_names)){cand_names <- letters[1:ncol(sims)]}
  for(i in 1:(ncol(sims)-1)){
    for(j in (i+1):ncol(sims)){
      out[[paste0(cand_names[i], sep, cand_names[j])]] = ab_plurality_tie_for_first_from_sims(cbind(sims[,c(i,j)], sims[,-c(i,j)]), tol = tol, n = n)
    }
  }
  out
}

ab_plurality_tie_for_first_from_sims <- function(sims, tol = .01, n = 1000){
  row_max <- apply(sims, 1, max)
  mean((sims[,1] == row_max | sims[,2] == row_max) & abs(sims[,1] - sims[,2]) < tol/2)/(tol*n) # Had divided by sqrt because: "Normalization here because the width of a channel where a and b is separated by tol holding fixed others' votes is sqrt(2)*tol"
}


#' @export
positional_pivot_probs_froms_sims <- function(sims, tol = .01, n = 1000, s = .5, cand_names = NULL, sep = ""){

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  if(is.null(cand_names)){
    cand_names <- names(sims)
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
  }

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  normalization_factor <- n*tol

  # this is direct Monte Carlo -- nothing fancy
  out <- list(
    mean(score_a > score_b & score_a - score_b < tol & score_b > score_c)/normalization_factor,
    mean(score_a > score_c & score_a - score_c < tol & score_c > score_b)/normalization_factor,
    mean(score_b > score_c & score_b - score_c < tol & score_c > score_a)/normalization_factor
  )

  names(out) <- c(paste0(cand_names[1], sep, cand_names[2]), paste0(cand_names[1], sep, cand_names[3]), paste0(cand_names[2], sep, cand_names[3]))

  out

}


#' @export
condorcet_pivot_probs_from_sims <- function(sims, n = 1000, tol = .01, cand_names = NULL, sep = "", kemeny = T){

  if(is.null(cand_names)){
    cand_names <- letters[1:3]
  }

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  a_vs_b <- apply(sims[,c(1,2,5)], 1, sum) # the share a got against b
  a_vs_c <- apply(sims[,c(1,2,3)], 1, sum) # the share a got against c
  b_vs_c <- apply(sims[,c(1,3,4)], 1, sum) # the share b got against c

  normalizer <- tol*n

  out <- list(
    mean(a_vs_c > 1/2 & b_vs_c > 1/2 & abs(a_vs_b - 1/2) < tol/4)/normalizer,
    mean(a_vs_b > 1/2 & b_vs_c < 1/2 & abs(a_vs_c - 1/2) < tol/4)/normalizer,
    mean(a_vs_b < 1/2 & a_vs_c < 1/2 & abs(b_vs_c - 1/2) < tol/4)/normalizer
  )
  # why divide by 4?
  # we want cases where (v_a + v_{ca}) - (v_b + v_{cb}) \in (-tol/2, tol/2)
  # ie 2(v_a + v_{ca}) - 1 > -tol/2 &  2(v_a + v_{ca}) - 1 < tol/2
  # ie v_a + v_{ca} > 1/2 - tol/4 & v_a + v_{ca} < 1/2 + tol/4
  # ie v_a + v_{ca} - 1/2 \in (-tol/4, tol/4)

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
        a_vs_c < (1 - a_vs_b) &
        # b's loss to a worse than c's loss to b
        (1 - a_vs_b) < (1 - b_vs_c) &
        # but just barely
        (1 - b_vs_c) - (1 - a_vs_b) < tol)/normalizer

    out[[paste0(out_names[3], "_reverse")]] <- mean(
      # OR reverse cycle
      reverse_cycle &
        # a's loss to b worse than c's loss to a
        a_vs_b < (1 - a_vs_c) &
        # c's loss to a worse than b's loss to c
        (1 - a_vs_c) < b_vs_c &
        # but just barely
        b_vs_c - (1 - a_vs_c) < tol)/normalizer

  }

  out

}
