#' Faster Monte Carlo estimation of pivot probabilities
#'
#' These methods take simulated election results and return a list of
#' pivot events, each with a pivot probability (\code{integral}) and
#' P matrix (\code{P}). They produce identical results to
#' \code{election_event_probs()} (or nearly so) but they are faster because
#' they are written in a less general way.
#'
#' For plurality elections, any number of candidates can be specified.
#' The other voting systems can only handle three candidates.
#'
#' Results have been validated against \code{election_event_probs()}.
#'
#' @param sims A matrix of simulated election results, with one column per
#' ballot type. Must be 6 columns for the ordinal methods (IRV, Kemeny-Young, positional).
#' @param n Size of electorate
#' @param tol Window within which two candidates are considered to be tied for the purpose of
#' simulation. This is normalized away. Wider window means lower variance but more bias.
#' @param cand_names Names of the candidates.
#' @param sep Separation between candidate names.
#' @param kemeny For Condorcet method, compute Kemeny-Young pivot event probabilities? If \code{F}, only handles event in which the Condorcet winner
#' is decided.
#'
#' @examples
#' sims <- gtools::rdirichlet(100000, alpha = c(9,7,3,4,4,6))
#' plurality_pivot_probs_from_sims(sims)
#' positional_pivot_probs_from_sims(sims, s=.5)
#' irv_pivot_probs_from_sims(sims)
#' condorcet_pivot_probs_from_sims(sims)
#'
#'
#'@name standalone_monte_carlo_methods
NULL

#' @rdname standalone_monte_carlo_methods
#' @export
plurality_pivot_probs_from_sims <- function(sims = NULL, n = 1000, window = .01, cand_names = NULL, sep = "_", method = "density", merge = F){

  out <- list()
  if(is.null(cand_names)){cand_names <- letters[1:ncol(sims)]}
  for(i in 1:(ncol(sims)-1)){
    for(j in (i+1):ncol(sims)){
      pp <- ab_plurality_tie_for_first_from_sims(cbind(sims[,c(i,j), drop = F], sims[,-c(i,j), drop = F]), method = method, n = n, merge = merge, window = window)
      out[[paste0(cand_names[i], sep, cand_names[j])]] <- list(
        integral = pp[1],
        P = P_matrix_from_indices(i,j,ncol(sims))
      )
      out[[paste0(cand_names[j], sep, cand_names[i])]] <- list(
        integral = pp[2],
        P = P_matrix_from_indices(j,i,ncol(sims))
      )
    }
  }
  out
}

P_matrix_from_indices <- function(i,j,k){
  out <- matrix(0, nrow = k, ncol = k)
  out[i,] <- 1
  out[c(i,j),j] <- c(0,1)
  out
}

ab_plurality_tie_for_first_from_sims <- function(sims, method = "density", n = 1000, merge = F, window = .01){

  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    if(method == "density"){
      cond <- rep(T, nrow(sims))
      sum_12 <- sims[,1] + sims[,2]
      for(j in 3:ncol(sims)){
        cond <- cond & (sims[,j] < sum_12/2)
      }
    }else if(method == "naive_density"){
      rank_mat <- apply(-sims, 1, rank) %>% t() # not sure why I have to do the transpose
      cond <- rank_mat[,1] + rank_mat[,2] == 3
    }
    mean(cond)*ks::kde(x = (sims[,1] - sims[,2])[cond], eval.points = limits)$estimate*(1/n)
  }else if(method == "rectangular"){
    row_max <- apply(sims, 1, max)
    pp <- mean((sims[,1] == row_max | sims[,2] == row_max) & abs(sims[,1] - sims[,2]) < window/2)/(window*n)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for plurality pivot prob estimation: ", method, "\n")
  }

}


#' @rdname standalone_monte_carlo_methods
#' @export
positional_pivot_probs_from_sims <- function(sims, window = .01, n = 1000, s = .5, cand_names = NULL, sep = "_", method = "density", merge = F){

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

  # assemble the scores
  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  # go through the events
  out <- list()
  ab_P <- rbind(c(1,1,s,0,1,1-s),
        c(0,0,1-s,1,0,s),
        0)

  ab <- positional_pivot_probs_12_from_scores(score_a, score_b, score_c, method = method, n = n, merge = merge, window = window)

  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = ab[1],
    P = ab_P
  )

  out[[paste0(cand_names[2], sep, cand_names[1])]] <-
    list(integral = ab[2],
         P = ab_P[c(2,1,3), c(3,4,1,2,6,5)])

  ac <- positional_pivot_probs_12_from_scores(score_a, score_c, score_b, method = method, n = n, merge = merge, window = window)

  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = ac[1],
    P = ab_P[c(1,3,2), c(2,1,5,6,3,4)]
  )

  out[[paste0(cand_names[3], sep, cand_names[1])]] <-
    list(integral = ac[2],
         P = ab_P[c(2,3,1), c(4,3, 6,5, 1,2)])

  bc <- positional_pivot_probs_12_from_scores(score_b, score_c, score_a, method = method, n = n, merge = merge, window = window)

  out[[paste0(cand_names[2], sep, cand_names[3])]] <-
    list(integral = bc[1],
         P = ab_P[c(3,1,2), c(5,6,2,1,4,3)])

  out[[paste0(cand_names[3], sep, cand_names[2])]] <-
    list(integral = bc[2],
         P = ab_P[c(3,2,1), c(6,5,4,3,2,1)])

  out

}

positional_pivot_probs_12_from_scores <- function(score_1, score_2, score_3, method = "density", n = 1000, merge = F, window = .01){
  # returns a pair of results, the first for 1 being ahead of 2, the second for 2 being ahead of 1 (these are the same for merge = T or method = "rectangular")
  # methods: "density", "naive_density", "rectangular"
  # easily extended to more candidates.

  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    if(method == "density"){
      cond <- score_3 < .5*(score_1 + score_2)
    }else if(method == "naive_density"){
      cond <- score_1 > score_3 & score_2 > score_3
    }
    mean(cond)*ks::kde(x = (score_1 - score_2)[cond], eval.points = limits)$estimate*(1/n)
  }else if(method == "rectangular"){
    pp <- mean(abs(score_1 - score_2) < window/2 & score_1 - score_3 > 1/n & score_2 - score_3 > 1/n)/(n*window)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }

}

#' @rdname standalone_monte_carlo_methods
#' @export
irv_pivot_probs_from_sims <- function(sims, tol = .01, n = 1000, s = 0, cand_names = NULL, sep = "_"){

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

  # positional scores
  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  # pairwise margins
  a_vs_b <- 2*(apply(sims[,c(1,2,5)], 1, sum) - .5) # a and b
  a_vs_c <- 2*(apply(sims[,c(1,2,3)], 1, sum) - .5) # a and c
  b_vs_c <- 2*(apply(sims[,c(1,3,4)], 1, sum) - .5) # b and c

  normalization_factor <- n*tol
  out <- list()

  # second round pivot events
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = mean(score_a > score_c & score_b > score_c & abs(a_vs_b) < tol/2)/normalization_factor,
    P = rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1),0)
  )

  out[[paste0(cand_names[1], sep, cand_names[3])]] <- out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = mean(score_a > score_b & score_c > score_b & abs(a_vs_c) < tol/2)/normalization_factor,
    P = rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1))
  )

  out[[paste0(cand_names[2], sep, cand_names[3])]] <- out[[paste0(cand_names[3], sep, cand_names[2])]] <- list(
    integral = mean(score_b > score_a & score_c > score_a & abs(b_vs_c) < tol/2)/normalization_factor,
    P = rbind(0, c(1,0,1,1,0,0), c(0,1,0,0,1,1))
  )

  # first round pivot events -- so many of these!
  ab_tie <- score_c > score_a & score_c > score_b & abs(score_a - score_b) < tol/2
  ab_P <- rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0)
  ba_P <- rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0)

  abab <- mean(ab_tie & a_vs_c > 0 & b_vs_c > 0)/normalization_factor
  # a_b|ab
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = abab,
    P = ab_P
  )
  # b_a|ba
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = abab,
    P = ba_P
  )

  # a_b|ac
  abac <- mean(ab_tie & a_vs_c > 0 & b_vs_c < 0)/normalization_factor
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = abac,
    P = ab_P[c(1,3,2),]
  )
  # b_a|ca
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = abac,
    P = ba_P[c(1,3,2),]
  )

  # a_b|cb
  abcb <- mean(ab_tie & a_vs_c < 0 & b_vs_c > 0)/normalization_factor
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = abcb,
    P = ab_P[c(3,2,1),]
  )
  # b_a|bc
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = abcb,
    P = ba_P[c(3,2,1),]
  )

  ac_tie <- score_b > score_a & score_b > score_c & abs(score_a - score_c) < tol/2
  ac_P <- rbind(c(1,1,1,1,0,0), c(0,0,0,0,1,1), 0)
  ca_P <- rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0)

  acac <- mean(ac_tie & a_vs_b > 0 & b_vs_c < 0)/normalization_factor
  # a_c|ac
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = acac,
    P = ac_P
  )
  # c_a|ca
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = acac,
    P = ca_P
  )

  # a_c|ab
  acab <- mean(ac_tie & a_vs_b > 0 & b_vs_c > 0)/normalization_factor
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = acab,
    P = ac_P[c(1,3,2),]
  )
  # c_a|ba
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = acab,
    P = ca_P[c(1,3,2),]
  )

  # a_c|bc
  acbc <- mean(ac_tie & a_vs_b < 0 & b_vs_c < 0)/normalization_factor
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = acbc,
    P = ac_P[c(3,2,1),]
  )
  # c_a|cb
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = acbc,
    P = ca_P[c(3,2,1),]
  )


  bc_tie <- score_a > score_b & score_a > score_c & abs(score_b - score_c) < tol/2
  bc_P <- rbind(0, c(1,1,1,1,0,0), c(0,0,0,0,1,1))
  cb_P <- rbind(0, c(0,0,1,1,0,0), c(1,1,0,0,1,1))

  bcbc <- mean(bc_tie & a_vs_b < 0 & a_vs_c < 0)/normalization_factor
  # b_c|bc
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = bcbc,
    P = bc_P
  )
  # c_b|cb
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = bcbc,
    P = cb_P
  )

  # b_c|ba
  bcba <- mean(bc_tie & a_vs_b < 0 & a_vs_c > 0)/normalization_factor
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = bcba,
    P = bc_P[c(1,3,2),]
  )
  # c_b|ab
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = bcba,
    P = cb_P[c(1,3,2),]
  )

  # b_c|ac
  bcac <- mean(bc_tie & a_vs_b > 0 & a_vs_c < 0)/normalization_factor
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = bcac,
    P = bc_P[c(3,2,1),]
  )
  # c_b|ca
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = bcac,
    P = cb_P[c(3,2,1),]
  )

  out

}

#' @rdname standalone_monte_carlo_methods
#' @export
condorcet_pivot_probs_from_sims <- function(sims, n = 1000, tol = .01, cand_names = NULL, sep = "_", kemeny = T){

  if(is.null(cand_names)){
    if(kemeny & (cand_names %>% sort() %>% paste(collapse = "") != "abc")){
      warning("Kemeny-Young pivot event names will use a, b, c, not the supplied candidate names.")
    }
    cand_names <- letters[1:3]
  }

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  a_vs_b <- apply(sims[,c(1,2,5)], 1, sum) # the share a got against b
  a_vs_c <- apply(sims[,c(1,2,3)], 1, sum) # the share a got against c
  b_vs_c <- apply(sims[,c(1,3,4)], 1, sum) # the share b got against c

  normalizer <- tol*n

  out <- list()

  a_b <- mean(a_vs_c > 1/2 & b_vs_c > 1/2 & abs(a_vs_b - 1/2) < tol/4)/normalizer
  a_b_P <- rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1),0)
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = a_b,
    P = a_b_P
  )
  out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = a_b,
    P = a_b_P
  )

  a_c <- mean(a_vs_b > 1/2 & b_vs_c < 1/2 & abs(a_vs_c - 1/2) < tol/4)/normalizer
  a_c_P <- rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1))
  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = a_c,
    P = a_c_P
  )
  out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = a_c,
    P = a_c_P
  )

  b_c <- mean(a_vs_b < 1/2 & a_vs_c < 1/2 & abs(b_vs_c - 1/2) < tol/4)/normalizer
  b_c_P <- rbind(0, c(1,0,1,1,0,0), c(0,1,0,0,1,1))
  out[[paste0(cand_names[2], sep, cand_names[3])]] <- list(
    integral = b_c,
    P = b_c_P
  )
  out[[paste0(cand_names[3], sep, cand_names[2])]] <- list(
    integral = b_c,
    P = b_c_P
  )
  # why divide by 4?
  # we want cases where (v_a + v_{ca}) - (v_b + v_{cb}) \in (-tol/2, tol/2)
  # ie 2(v_a + v_{ca}) - 1 > -tol/2 &  2(v_a + v_{ca}) - 1 < tol/2
  # ie v_a + v_{ca} > 1/2 - tol/4 & v_a + v_{ca} < 1/2 + tol/4
  # ie v_a + v_{ca} - 1/2 \in (-tol/4, tol/4)

  if(kemeny){
    forward_cycle <- a_vs_b > .5 & a_vs_c < .5 & b_vs_c > .5
    reverse_cycle <- a_vs_b < .5 & a_vs_c > .5 & b_vs_c < .5

    ab_forward <- mean(
      # forward cycle, and
      forward_cycle &
        # c's loss to b worse than b's loss to a
        b_vs_c > a_vs_b &
        # b's loss to a about the same as a's loss to c
        abs(a_vs_b - (1 - a_vs_c)) < tol/2)/normalizer
    out[["ac_ba|abca"]] <- list(
      integral = ab_forward,
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(0,1,0), c(1,0,0), c(0,1,0))
    )
    out[["ba_ac|abca"]] <- list(
      integral = ab_forward,
      P = cbind(c(1,0,0), c(1,0,0), c(0,1,0), c(0,1,0), c(0,1,0), c(0,1,0))
    )

    ab_reverse <- mean(
      # reverse cycle
      reverse_cycle &
        # c's loss to a worse than a's loss to b
        a_vs_c > (1 - a_vs_b) &
        # b's loss to c about the same as a's loss to b
        abs((1 - b_vs_c) - (1 - a_vs_b)) < tol/2)/normalizer
    out[["ab_bc|bacb"]] <- list(
      integral = ab_reverse, # b needs a losing to b but b beating c: bac bca
      P = cbind(c(1,0,0), c(1,0,0), c(0,1,0), c(0,1,0), c(1,0,0), c(1,0,0))
    )
    out[["bc_ab|bacb"]] <- list(
      integral = ab_reverse,  # a needs b losing to c but a beating b: acb, cab
      P = cbind(c(0,1,0), c(1,0,0), c(0,1,0), c(0,1,0), c(1,0,0), c(0,1,0))
    )

    ac_forward <- mean(
      # forward cycle, and
      forward_cycle &
        # b's loss to a worse than c's loss to b
        a_vs_b > b_vs_c &
        # c's loss to b about the same as a's loss to c
        abs(b_vs_c - (1 - a_vs_c)) < tol/2)/normalizer
    out[["ac_cb|cabc"]] <- list(
      integral = ac_forward, # c needs to beat both: cba cab
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(1,0,0), c(0,0,1), c(0,0,1))
    )
    out[["cb_ac|cabc"]] <- list(
      integral = ac_forward,  # a needs c to finish last: abc, bac
      P = cbind(c(1,0,0), c(0,0,1), c(1,0,0), c(0,0,1), c(0,0,1), c(0,0,1))
    )

    ac_reverse <- mean(
      # reverse cycle
      reverse_cycle &
        # b's loss to c worse than a's loss to b
        (1 - b_vs_c) > (1 - a_vs_b) &
        # a's loss to b about the same as c's loss to a
        abs((1 - a_vs_b) - a_vs_c) < tol/2)/normalizer
    out[["ab_ca|acba"]] <- list(
      integral = ac_reverse, # c needs a to finish last: cba, bca
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(0,0,1), c(1,0,0), c(0,0,1))
    )
    out[["ca_ab|acba"]] <- list(
      integral = ac_reverse,  # a needs to beat both: abc, acb
      P = cbind(c(1,0,0), c(1,0,0), c(0,0,1), c(0,0,1), c(0,0,1), c(0,0,1))
    )

    bc_forward <- mean(
      # forward cycle, and
      forward_cycle &
        # a's loss to c worse than b's loss to a
        a_vs_c < (1 - a_vs_b) &
        # b's loss to a about the same as c's loss to b
        abs((1 - a_vs_b) - (1 - b_vs_c)) <  tol/2)/normalizer
    out[["ba_cb|bcab"]] <- list(
      integral = bc_forward, # c needs b last: acb, cab
      P = cbind(c(0,1,0), c(0,0,1), c(0,1,0), c(0,1,0), c(0,0,1), c(0,1,0))
    )
    out[["cb_ba|bcab"]] <- list(
      integral = bc_forward,  # b needs to beat both: bca, bac
      P = cbind(c(0,0,1), c(0,0,1), c(0,1,0), c(0,1,0), c(0,0,1), c(0,0,1))
    )

    bc_reverse <- mean(
      # OR reverse cycle
      reverse_cycle &
        # a's loss to b worse than c's loss to a
        a_vs_b < (1 - a_vs_c) &
        # c's loss to a about the same as b's loss to c
        abs((1 - a_vs_c) - b_vs_c) < tol/2)/normalizer
    out[["bc_ca|cbac"]] <- list(
      integral = bc_reverse, # c needs to beat both: cba, cab
      P = cbind(c(0,1,0), c(0,1,0), c(0,1,0), c(0,1,0), c(0,0,1), c(0,0,1))
    )
    out[["ca_bc|cbac"]] <- list(
      integral = bc_reverse,  # b needs c to finish last: bac, abc
      P = cbind(c(0,1,0), c(0,0,1), c(0,1,0), c(0,0,1), c(0,0,1), c(0,0,1))
    )
  }

  out

}
