## functions for recovering events in each system

plurality_event_list <- function(n = 1000, k = 4, max_pivot_event_degree = 2){

  if(!max_pivot_event_degree %in% 1:2){
    stop("max_pivot_event_degree needs to be either 1 or 2. The event_probabilities_from_event_list() function only permutes three candidates, so we cannot handle pivot events of degree 3 or higher, e.g. four-way ties.")
  }

  el <- list()

  el[["params"]] <- list(
    n = n,
    ordinal = F
  )

  win_conditions <- cbind(1, -diag(k-1), 1/n) # i beats j, k, ...
  win_P <- matrix(0, nrow = k, ncol = k)
  win_P[1, ] <- 1

  el[["i_"]] <- list(
    conditions = win_conditions,
    rows_to_alter = c(),
    P = win_P
  )

  this_P <- win_P
  this_P[1:2, 2] <- c(0,1)

  el[["i_j"]] <- list(
    conditions = el[["i_"]]$conditions,
    rows_to_alter = c(1),                      # i barely beats j
    scaling_factor = sqrt(2),
    P = this_P,
    adjacent_events = "j_i"
  )

  if(max_pivot_event_degree == 2){

    this_P <- win_P
    this_P[1:2, 2] <- c(0,1)
    this_P[1:3, 3] <- c(0,0,1)

    el[["i_jk"]] <- list(
      conditions = el[["i_"]]$conditions,
      rows_to_alter = c(1,2),                    # i barely beats j and k
      P = this_P,
      adjacent_events = c("j_ik", "k_ij") # note "k_ij" and "k_ji" both get stored as "k_ij"
    )
  }

  el
}

positional_event_list <- function(n = 1000, s = .5){

  el <- list()

  el[["params"]] <- list(
    n = n,
    ordinal = T
  )

  el[["i_"]] <- list(
    conditions = rbind(c(1-s, 1, s-1, -1, s, -s, 1/n),  # i scores above j
                     c(1, 1-s, s, -s, s-1, -1, 1/n)),   # i scores above k
    rows_to_alter = c(),
    P = rbind(rep(1, 6), 0, 0, 0)
    )

  el[["i_j"]] <- list(
    conditions = el[["i_"]]$conditions,
    rows_to_alter = c(1),                               # i barely scores above j
    scaling_factor = 2*sqrt(1 - s + s^2),
    P = rbind(c(1,1,s,0,1,1-s),
            c(0,0,1-s,1,1,s),
            0),
    adjacent_events = "j_i"
    )

  # not doing second-degree pivot events (a point), but could add -- would need to get scaling factor

  el
}

irv_event_list <- function(n = 1000, s = 0){

  # first round is positional, so could be Borda Count, plurality, anti-plurality, etc

  irv_events <- list()

  irv_events[["params"]] <- list(
    n = n,
    ordinal = T
  )

  # positional scaling factor
  psf <- 2*sqrt(1 - s + s^2)

  # i wins over j in the second round
  irv_events[["i__j"]] <- list(
    conditions = rbind(c(1, 1-s,s,-s,s-1,-1, 1/n),  # i beats k in first round
                       c(s,-s,1,1-s,-1,s-1,1/n),    # so does j
                       c(1,1,-1,-1,1,-1, 1/n)), # i beats j in second round
    rows_to_alter = c(),
    P = rbind(rep(1, 6), 0, 0),
    adjacent_events = "j__i"
  )

  # second round pivot event
  irv_events[["i_j"]] <- list(
    conditions = irv_events[["i__j"]]$conditions,
    rows_to_alter = c(3),                        # i barely beats j in second round
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0),
    adjacent_events = "j_i")

  # first round pivot events
  irv_events[["i_j|ij"]] <- list(
    conditions = rbind(c(-1, s-1,-s,s,1-s,1, 1/n),    # k beats i in first round
                       c(-s,s,-1,s-1,1,1-s, 1/n),      # k beats j in first round
                       c(1-s,1,s-1,-1,s,-s, 1/n),     # i beats j in first round
                       c(1,1,1,-1,-1,-1, 1/n),    # i beats k in second round
                       c(1,-1,1,1,-1,-1, 1/n)),   # j beats k in second round
    rows_to_alter = c(3),                         # i barely beats j in first round
    scaling_factor = psf,
    P = rbind(c(1,1,s,0,1,1-s), c(0,0,1-s,1,0,s), 0),
    adjacent_events = "j_i|ji")

  i_j_ik_conditions <- irv_events[["i_j|ij"]]$conditions
  i_j_ik_conditions[5,1:6] <- -i_j_ik_conditions[5,1:6] # reverse the j k ordering

  irv_events[["i_j|ik"]] <- list(
    conditions = i_j_ik_conditions,
    rows_to_alter = c(3),
    scaling_factor = psf,
    P = rbind(c(1,1,s,0,1,1-s), 0, c(0,0,1-s,1,0,s)),
    adjacent_events = "j_i|ki")

  i_j_kj_conditions <- irv_events[["i_j|ij"]]$conditions
  i_j_kj_conditions[4,1:6] <- -i_j_kj_conditions[4,1:6]  # reverse the i k ordering

  irv_events[["i_j|kj"]] <- list(
    conditions = i_j_kj_conditions,
    rows_to_alter = c(3),
    scaling_factor = psf,
    P = rbind(0, c(0,0,1-s,1,0,s), c(1,1,s,0,1,1-s)),
    adjacent_events = "j_i|jk")

  irv_events
}

kemeny_young_event_list <- function(n = 1000){

  el <- list()

  el[["params"]] <- list(
    n = n,
    ordinal = T
  )

  # i is the condorcet winner, with j behind.
  el[["i__j"]] <- list(
    conditions = rbind(c(1, 1, -1, -1, 1, -1, 1/n),   # i pairwise beats j
                       c(1, 1, 1, -1, -1, -1, 1/n),   # i pairwise beats k
                       c(1, -1, 1, 1, -1, -1, 1/n)),  # j pairwise beats k
    rows_to_alter = c(),
    P = rbind(rep(1, 6), 0, 0)
  )

  # i is barely the condorcet winner, over j.
  el[["i_j"]] <- list(
    conditions = el[["i__j"]]$conditions,
    rows_to_alter = c(1),  # i barely beats j
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0),
    adjacent_events = "j_i"
  )

  # i wins in a cycle
  el[["ijki"]] = list(
    conditions = rbind(c(1,1,-1,-1,1,-1, 1/n),    # i pairwise beats j
                       c(1,-1,1,1,-1,-1, 1/n),    # j pairwise beats k
                       c(-1, -1, -1, 1, 1, 1, 1/n), # k pairwise beats i
                       c(1,1,0,-1,0,-1,1/n),      # i's loss to k better than j's loss to i
                       c(0,-1, 1, 1,-1, 0, 1/n)   # j's loss to i better than k's loss to j
    ),
    rows_to_alter = c(),
    P = rbind(rep(1, 6), 0, 0)
  )

  # i barely wins in a cycle, having a slightly better loss than j
  el[["ik_ji|ijki"]] = list(
    conditions = el[["ijki"]]$conditions,
    rows_to_alter = c(4), # i's loss to k barely better than j's loss to i
    scaling_factor = 2,  # is it sqrt(number of tallies involved)?
    # P matrix is tricky.
    # a ballot that improves both i relative to k and j relative to i (jik) or makes both worse (kij) does not change the outcome. you only elect j by putting jki or kji.
    P = rbind(c(1,1,1,0,1,0),
              c(0,0,0,1,0,1),
              0),
    adjacent_events = "ji_ik|ijki"
  )

  # reverse the i j comparison
  ji_ik_conditions <- el[["ijki"]]$conditions
  ji_ik_conditions[4,1:6] <- -ji_ik_conditions[4,1:6]

  el[["ji_ik|ijki"]] = list(
    conditions = ji_ik_conditions,
    rows_to_alter = c(4), # i's loss to k barely better than j's loss to i
    scaling_factor = 2,  # is it sqrt(number of tallies involved)?
    # P matrix is tricky.
    # a ballot that improves both i relative to k and j relative to i (jik) or makes both worse (kij) does not change the outcome. you only elect j by putting jki or kji.
    P = rbind(c(1,1,1,0,1,0),
              c(0,0,0,1,0,1),
              0),
    adjacent_events = "ik_ji|ijki"
  )

  el

}
