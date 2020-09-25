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
#' @param window Window within which two candidates are considered to be tied for the purpose of
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

density_estimate <- function(x, bw_divisor = 1, eval.points = c(0)){
  if(length(x) <= 1){return(rep(0, length(eval.points)))} # can't get a bandwidth with only 1 point
  bw <- ks::hpi(x, binned = T)/bw_divisor
  ks::kde(x = x, h = bw, eval.points = eval.points)$estimate
}

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
        P = plurality_P_matrix_from_indices(i,j,ncol(sims))
      )
      out[[paste0(cand_names[j], sep, cand_names[i])]] <- list(
        integral = pp[2],
        P = plurality_P_matrix_from_indices(j,i,ncol(sims))
      )
    }
  }
  out
}

plurality_P_matrix_from_indices <- function(i,j,k){
  out <- matrix(0, nrow = k, ncol = k)
  out[i,] <- 1
  out[c(i,j),j] <- c(0,1)
  out
}

ab_plurality_tie_for_first_from_sims <- function(sims, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 1){

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
        cond <- cond & (sims[,j] < sum_12/2 - 1/n)
      }
    }else if(method == "naive_density"){
      cond <- rep(T, nrow(sims))
      for(j in 3:ncol(sims)){
        cond <- cond & (sims[,j] < sims[,1] - 1/n & sims[,j] < sims[,2] - 1/n)
      }
    }
    the_density <- density_estimate(x = (sims[,1] - sims[,2])[cond], eval.points = limits, bw_divisor = bw_divisor)
    mean(cond)*the_density*(1/n)
    # mean(cond)*ks::kde(x = (sims[,1] - sims[,2])[cond], eval.points = limits)$estimate*(1/n)
  }else if(method == "rectangular"){
    cond <- rep(T, nrow(sims))
    for(j in 3:ncol(sims)){
      cond <- cond & (sims[,j] < sims[,1] - 1/n & sims[,j] < sims[,2] - 1/n)
    }
    pp <- mean(cond & abs(sims[,1] - sims[,2]) < window/2)/(window*n)
    # row_max <- apply(sims, 1, max)
    # pp <- mean((sims[,1] == row_max | sims[,2] == row_max) & abs(sims[,1] - sims[,2]) < window/2)/(window*n)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for plurality pivot prob estimation: ", method, "\n")
  }

}


#' @rdname standalone_monte_carlo_methods
#' @export
positional_pivot_probs_from_sims <- function(sims, window = .01, n = 1000, s = .5, cand_names = NULL, sep = "_", method = "density", merge = F, bw_divisor = 1){

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

  ab <- positional_pivot_probs_12_from_scores(score_a, score_b, score_c, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)

  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = ab[1],
    P = ab_P
  )

  out[[paste0(cand_names[2], sep, cand_names[1])]] <-
    list(integral = ab[2],
         P = ab_P[c(2,1,3), c(3,4,1,2,6,5)])

  ac <- positional_pivot_probs_12_from_scores(score_a, score_c, score_b, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)

  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = ac[1],
    P = ab_P[c(1,3,2), c(2,1,5,6,3,4)]
  )

  out[[paste0(cand_names[3], sep, cand_names[1])]] <-
    list(integral = ac[2],
         P = ab_P[c(2,3,1), c(4,3, 6,5, 1,2)])

  bc <- positional_pivot_probs_12_from_scores(score_b, score_c, score_a, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)

  out[[paste0(cand_names[2], sep, cand_names[3])]] <-
    list(integral = bc[1],
         P = ab_P[c(3,1,2), c(5,6,2,1,4,3)])

  out[[paste0(cand_names[3], sep, cand_names[2])]] <-
    list(integral = bc[2],
         P = ab_P[c(3,2,1), c(6,5,4,3,2,1)])

  out

}

positional_pivot_probs_12_from_scores <- function(score_1, score_2, score_3, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 2){
  # returns a pair of results, the first for 1 being ahead of 2, the second for 2 being ahead of 1 (these are the same for merge = T or method = "rectangular")
  # methods: "density", "naive_density", "rectangular"
  # easily extended to more candidates.
  # bw_divisor > 1 undersmooths the density estimate

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
    the_density <- density_estimate(x = (score_1 - score_2)[cond], eval.points = limits, bw_divisor = bw_divisor)
    mean(cond)*the_density*(1/n)
  }else if(method == "rectangular"){
    pp <- mean(abs(score_1 - score_2) < window/2 & score_1 - score_3 > 1/n & score_2 - score_3 > 1/n)/(n*window)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }

}

#' @rdname standalone_monte_carlo_methods
#' @export
irv_pivot_probs_from_sims <- function(sims, window = .01, n = 1000, s = 0, cand_names = NULL, sep = "_", method = "density", merge = F, bw_divisor = 1){

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

  out <- list()

  # second round pivot events
  # a_b and b_a
  ab_P <- rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1),0)
  ab_pp <- irv_second_round_pivot_prob_ab(score_a, score_b, score_c, a_vs_b, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = ab_pp[1],
    P = ab_P
  )
  out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = ab_pp[2],
    P = ab_P
  )

  # a_c and c_a
  ac_P <- rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1))
  ac_pp <- irv_second_round_pivot_prob_ab(score_a, score_c, score_b, a_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = ac_pp[1],
    P = ac_P
  )
  out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = ac_pp[2],
    P = ac_P
  )

  # b_c and c_b
  bc_P <- rbind(0, c(1,0,1,1,0,0), c(0,1,0,0,1,1))
  bc_pp <- irv_second_round_pivot_prob_ab(score_b, score_c, score_a, b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[2], sep, cand_names[3])]] <- list(
    integral = bc_pp[1],
    P = bc_P
  )
  out[[paste0(cand_names[3], sep, cand_names[2])]] <- list(
    integral = bc_pp[2],
    P = bc_P
  )

  # first-round pivot events -- one set for each pair of candidates

  # a and b
  ab_pps <- irv_first_round_pivot_probs_ab(score_a, score_b, score_c, a_vs_c, b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  ab_P <- rbind(c(1,1,0,0,1,1), c(0,0,1,1,0,0), 0)
  ba_P <- rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), 0)

  # a_b|ab
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = ab_pps[["i_j|ij"]][1],
    P = ab_P
  )
  # b_a|ba
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = ab_pps[["i_j|ij"]][2],
    P = ba_P
  )

  # a_b|ac
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = ab_pps[["i_j|ik"]][1],
    P = ab_P[c(1,3,2),]
  )
  # b_a|ca
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = ab_pps[["i_j|ik"]][2],
    P = ba_P[c(1,3,2),]
  )

  # a_b|cb
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = ab_pps[["i_j|kj"]][1],
    P = ab_P[c(3,2,1),]
  )
  # b_a|bc
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = ab_pps[["i_j|kj"]][2],
    P = ba_P[c(3,2,1),]
  )


  # a and c
  ac_pps <- irv_first_round_pivot_probs_ab(score_a, score_c, score_b, a_vs_b, -b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  ac_P <- rbind(c(1,1,1,1,0,0), 0, c(0,0,0,0,1,1))
  ca_P <- rbind(c(1,1,0,0,0,0), 0, c(0,0,1,1,1,1))

  # a_c|ac
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = ac_pps[["i_j|ij"]][1],
    P = ac_P
  )
  # b_a|ba
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = ac_pps[["i_j|ij"]][2],
    P = ca_P
  )

  # a_c|ab
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = ac_pps[["i_j|ik"]][1],
    P = ac_P[c(1,3,2),]
  )
  # c_a|ba
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = ac_pps[["i_j|ik"]][2],
    P = ca_P[c(1,3,2),]
  )

  # a_c|bc
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = ac_pps[["i_j|kj"]][1],
    P = ac_P[c(2,1,3),]
  )
  # c_a|cb
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = ac_pps[["i_j|kj"]][2],
    P = ca_P[c(2,1,3),]
  )


  # b and c
  bc_pps <- irv_first_round_pivot_probs_ab(score_b, score_c, score_a, -a_vs_b, -a_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  bc_P <- rbind(0, c(1,1,1,1,0,0), c(0,0,0,0,1,1))
  cb_P <- rbind(0, c(0,0,1,1,0,0), c(1,1,0,0,1,1))

  # b_c|bc
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = bc_pps[["i_j|ij"]][1],
    P = bc_P
  )
  # c_b|cb
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = bc_pps[["i_j|ij"]][2],
    P = cb_P
  )

  # b_c|ba
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = bc_pps[["i_j|ik"]][1],
    P = bc_P[c(3,2,1),]
  )
  # c_b|ab
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = bc_pps[["i_j|ik"]][2],
    P = bc_P[c(3,2,1),]
  )

  # b_c|ac
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = bc_pps[["i_j|kj"]][1],
    P = bc_P[c(2,1,3),]
  )
  # c_b|ca
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = bc_pps[["i_j|kj"]][2],
    P = cb_P[c(2,1,3),]
  )

  out
}

irv_second_round_pivot_prob_ab <- function(score_a, score_b, score_c, a_vs_b, s, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 1){
  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    if(method == "density"){
      cond_1 <- score_a - (a_vs_b/2)*(1 - s/2) - score_c > 1/n
      cond_2 <- score_b + (a_vs_b/2)*(1 - s/2) - score_c > 1/n
    }else if(method == "naive_density"){
      cond_1 <- score_a - score_c > 1/n
      cond_2 <- score_b - score_c > 1/n
    }
    cond <- cond_1 & cond_2
    the_density <- density_estimate(x = a_vs_b[cond], eval.points = limits, bw_divisor = bw_divisor)
    mean(cond)*the_density*(1/n)
  }else if(method == "rectangular"){
    pp <- mean(score_a - score_c > 1/n & score_b - score_c > 1/n & abs(a_vs_b) < window/2)/(n*window)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }
}

irv_first_round_pivot_probs_ab <- function(score_a, score_b, score_c, a_vs_c, b_vs_c, s, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 1){
  out <- list()
  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    delta <- score_a - score_b
    if(method == "density"){
      cond_1 <- score_c - score_a > 1/n - (delta/2)*(1 - s*(1-s)/2)
      cond_2 <- score_c - score_b > 1/n + (delta/2)*(1 - s*(1-s)/2)
      diff_4 <- a_vs_c - (delta/2)*(1 + s)
      diff_5 <- b_vs_c + (delta/2)*(1 + s)
      a_b_ab_cond <- cond_1 & cond_2 & diff_4 > 1/n & diff_5 > 1/n
      a_b_cb_cond <- cond_1 & cond_2 & diff_4 < -1/n & diff_5 > 1/n
      a_b_ac_cond <- cond_1 & cond_2 & diff_4 > 1/n & diff_5 < -1/n
    }else if(method == "naive_density"){
      c_first_cond <- score_c - score_a > 1/n & score_c - score_b > 1/n
      a_b_ab_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c > 1/n
      a_b_cb_cond <- c_first_cond & a_vs_c < -1/n & b_vs_c > 1/n
      a_b_ac_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c < -1/n
    }
    out[["i_j|ij"]] <- mean(a_b_ab_cond)*density_estimate(x = delta[a_b_ab_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
    out[["i_j|kj"]] <- mean(a_b_cb_cond)*density_estimate(x = delta[a_b_cb_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
    out[["i_j|ik"]] <- mean(a_b_ac_cond)*density_estimate(x = delta[a_b_ac_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
  }else if(method == "rectangular"){
    c_first_cond <- score_c - score_a > 1/n & score_c - score_b > 1/n
    a_b_ab_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c > 1/n
    a_b_cb_cond <- c_first_cond & a_vs_c < -1/n & b_vs_c > 1/n
    a_b_ac_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c < -1/n
    ab_tie_cond <- abs(score_a - score_b) < window/2
    out[["i_j|ij"]] <- rep(mean(a_b_ab_cond & ab_tie_cond)/(n*window), 2)
    out[["i_j|kj"]] <- rep(mean(a_b_cb_cond & ab_tie_cond)/(n*window), 2)
    out[["i_j|ik"]] <- rep(mean(a_b_ac_cond & ab_tie_cond)/(n*window), 2)
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }
  out
}


#' @rdname standalone_monte_carlo_methods
#' @export
condorcet_pivot_probs_from_sims <- function(sims, n = 1000, window = .01, cand_names = NULL, sep = "_", kemeny = T){

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

  normalizer <- window*n

  out <- list()

  a_b <- mean(a_vs_c > 1/2 & b_vs_c > 1/2 & abs(a_vs_b - 1/2) < window/4)/normalizer
  a_b_P <- rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1),0)
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = a_b,
    P = a_b_P
  )
  out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = a_b,
    P = a_b_P
  )

  a_c <- mean(a_vs_b > 1/2 & b_vs_c < 1/2 & abs(a_vs_c - 1/2) < window/4)/normalizer
  a_c_P <- rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1))
  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = a_c,
    P = a_c_P
  )
  out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = a_c,
    P = a_c_P
  )

  b_c <- mean(a_vs_b < 1/2 & a_vs_c < 1/2 & abs(b_vs_c - 1/2) < window/4)/normalizer
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
  # we want cases where (v_a + v_{ca}) - (v_b + v_{cb}) \in (-window/2, window/2)
  # ie 2(v_a + v_{ca}) - 1 > -window/2 &  2(v_a + v_{ca}) - 1 < window/2
  # ie v_a + v_{ca} > 1/2 - window/4 & v_a + v_{ca} < 1/2 + window/4
  # ie v_a + v_{ca} - 1/2 \in (-window/4, window/4)

  if(kemeny){
    forward_cycle <- a_vs_b > .5 & a_vs_c < .5 & b_vs_c > .5
    reverse_cycle <- a_vs_b < .5 & a_vs_c > .5 & b_vs_c < .5

    ab_forward <- mean(
      # forward cycle, and
      forward_cycle &
        # c's loss to b worse than b's loss to a
        b_vs_c > a_vs_b &
        # b's loss to a about the same as a's loss to c
        abs(a_vs_b - (1 - a_vs_c)) < window/2)/normalizer
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
        abs((1 - b_vs_c) - (1 - a_vs_b)) < window/2)/normalizer
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
        abs(b_vs_c - (1 - a_vs_c)) < window/2)/normalizer
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
        abs((1 - a_vs_b) - a_vs_c) < window/2)/normalizer
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
        abs((1 - a_vs_b) - (1 - b_vs_c)) <  window/2)/normalizer
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
        abs((1 - a_vs_c) - b_vs_c) < window/2)/normalizer
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
