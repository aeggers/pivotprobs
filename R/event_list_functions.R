#' Create election list for analysis
#'
#' An election list represents an election
#' in a minimal way: electorate size, ballot format, and
#' rules for determining the winner, electorate size.
#' It is passed to the \code{election_event_probabilities()}
#' function to compute the probability of each election event.
#'
#' An election list contains the electorate size \code{n},
#' an \code{ordinal}
#' flag indicating whether this is an ordinal or single-ballot
#' voting
#' system, and a list of election \code{events}. For each event,
#' the election list indicates \itemize{\item conditions under which the event takes place, and \item which candidate wins given that this
#' event takes place as a function of which additional ballot is submitted.}
#'
#' These functions produce plain lists, so functions can be written
#' for further election methods.
#'
#' @param n The size of the electorate. This (slightly) affects the
#' conditions for election events in the election list
#' (e.g. the region in which one
#' candidate is more than one vote ahead of another).
#' \code{n} is also used when computing pivot event probabilities
#' from an election list in \code{election_event_probabilities()},
#' where higher \code{n} means lower probability.
#' @param k The number of candidates in plurality.
#' @param s The value of a second-ranking in a three-candidate
#' positional election, where a top ranking is worth 1 and a bottom ranking is worth 0.
#'
#' @name election_list_functions
NULL

#' @rdname election_list_functions
#' @export
plurality_election <- function(n = 1000, k = 4, max_pivot_event_degree = 2){

  if(!max_pivot_event_degree %in% 1:2){
    stop("max_pivot_event_degree needs to be either 1 or 2. The event_probabilities_from_event_list() function only permutes three candidates, so we cannot handle pivot events of degree 3 or higher, e.g. four-way ties.")
  }

  el <- list()

  el$n <- n
  el$ordinal <- F
  el$events <- list()

  win_P <- matrix(0, nrow = k, ncol = k)
  win_P[1, ] <- 1

  el[["events"]][["i_"]] <- list(
    win_conditions = cbind(1, -diag(k-1), 1/n), # i beats j, k, ...
    tie_condition_rows = c(),
    P = win_P
  )

  this_P <- win_P
  this_P[1:2, 2] <- c(0,1)

  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i_"]]$win_conditions,
    tie_condition_rows = c(1),                      # i barely beats j
    scaling_factor = sqrt(2),
    P = this_P,
    adjacent_events = "j_i"
  )

  if(max_pivot_event_degree == 2){

    this_P <- win_P
    this_P[1:2, 2] <- c(0,1)
    this_P[1:3, 3] <- c(0,0,1)

    el[["events"]][["i_jk"]] <- list(
      win_conditions = el[["events"]][["i_"]]$win_conditions,
      tie_condition_rows = c(1,2),                    # i barely beats j and k
      P = this_P,
      adjacent_events = c("j_ik", "k_ij") # note "k_ij" and "k_ji" both get stored as "k_ij"
    )
  }

  el
}

#' @rdname election_list_functions
#' @export
positional_election <- function(n = 1000, s = .5){

  el <- list()

  el$n <- n
  el$ordinal <- T
  el$events <- list()

  el[["events"]][["i_"]] <- list(
    win_conditions = rbind(c(1-s, 1, s-1, -1, s, -s, 1/n),  # i scores above j
                     c(1, 1-s, s, -s, s-1, -1, 1/n)),   # i scores above k
    tie_condition_rows = c(),
    P = rbind(rep(1, 6), 0, 0, 0)
    )

  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i_"]]$win_conditions,
    tie_condition_rows = c(1),                               # i barely scores above j
    scaling_factor = 2*sqrt(1 - s + s^2),
    P = rbind(c(1,1,s,0,1,1-s),
            c(0,0,1-s,1,1,s),
            0),
    adjacent_events = "j_i"
    )

  # not doing second-degree pivot events (a point), but could add -- would need to get scaling factor

  el
}

#' @rdname election_list_functions
#' @export
irv_election <- function(n = 1000, s = 0){

  # first round is positional, so could be Borda Count, plurality, anti-plurality, etc

  el <- list()

  el$n <- n
  el$ordinal <- T
  el$events <- list()

  # positional scaling factor
  psf <- 2*sqrt(1 - s + s^2)

  # i wins over j in the second round
  el[["events"]][["i__j"]] <- list(
    win_conditions = rbind(c(1, 1-s,s,-s,s-1,-1, 1/n),  # i beats k in first round
                       c(s,-s,1,1-s,-1,s-1,1/n),    # so does j
                       c(1,1,-1,-1,1,-1, 1/n)), # i beats j in second round
    tie_condition_rows = c(),
    P = rbind(rep(1, 6), 0, 0),
    adjacent_events = "j__i"
  )

  # second round pivot event
  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i__j"]]$win_conditions,
    tie_condition_rows = c(3),                        # i barely beats j in second round
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0),
    adjacent_events = "j_i")

  # first round pivot events
  el[["events"]][["i_j|ij"]] <- list(
    win_conditions = rbind(c(-1, s-1,-s,s,1-s,1, 1/n),    # k beats i in first round
                       c(-s,s,-1,s-1,1,1-s, 1/n),      # k beats j in first round
                       c(1-s,1,s-1,-1,s,-s, 1/n),     # i beats j in first round
                       c(1,1,1,-1,-1,-1, 1/n),    # i beats k in second round
                       c(1,-1,1,1,-1,-1, 1/n)),   # j beats k in second round
    tie_condition_rows = c(3),                         # i barely beats j in first round
    scaling_factor = psf,
    P = rbind(c(1,1,s,0,1,1-s), c(0,0,1-s,1,0,s), 0),
    adjacent_events = "j_i|ji")

  i_j_ik_conditions <- el[["events"]][["i_j|ij"]]$win_conditions
  i_j_ik_conditions[5,1:6] <- -i_j_ik_conditions[5,1:6] # reverse the j k ordering

  el[["events"]][["i_j|ik"]] <- list(
    win_conditions = i_j_ik_conditions,
    tie_condition_rows = c(3),
    scaling_factor = psf,
    P = rbind(c(1,1,s,0,1,1-s), 0, c(0,0,1-s,1,0,s)),
    adjacent_events = "j_i|ki")

  i_j_kj_conditions <- el[["events"]][["i_j|ij"]]$win_conditions
  i_j_kj_conditions[4,1:6] <- -i_j_kj_conditions[4,1:6]  # reverse the i k ordering

  el[["events"]][["i_j|kj"]] <- list(
    win_conditions = i_j_kj_conditions,
    tie_condition_rows = c(3),
    scaling_factor = psf,
    P = rbind(0, c(0,0,1-s,1,0,s), c(1,1,s,0,1,1-s)),
    adjacent_events = "j_i|jk")

  el
}

#' @rdname election_list_functions
#' @export
kemeny_young_election <- function(n = 1000){

  el <- list()

  el$n <- n
  el$ordinal <- T
  el$events <- list()

  # i is the condorcet winner, with j behind.
  el[["events"]][["i__j"]] <- list(
    win_conditions = rbind(c(1, 1, -1, -1, 1, -1, 1/n),   # i pairwise beats j
                       c(1, 1, 1, -1, -1, -1, 1/n),   # i pairwise beats k
                       c(1, -1, 1, 1, -1, -1, 1/n)),  # j pairwise beats k
    tie_condition_rows = c(),
    P = rbind(rep(1, 6), 0, 0)
  )

  # i is barely the condorcet winner, over j.
  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i__j"]]$win_conditions,
    tie_condition_rows = c(1),  # i barely beats j
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1), 0),
    adjacent_events = "j_i"
  )

  # i wins in a cycle
  el[["events"]][["ijki"]] = list(
    win_conditions = rbind(c(1,1,-1,-1,1,-1, 1/n),    # i pairwise beats j
                       c(1,-1,1,1,-1,-1, 1/n),    # j pairwise beats k
                       c(-1, -1, -1, 1, 1, 1, 1/n), # k pairwise beats i
                       c(1,1,0,-1,0,-1,1/n),      # i's loss to k better than j's loss to i
                       c(0,-1, 1, 1,-1, 0, 1/n)   # j's loss to i better than k's loss to j
    ),
    tie_condition_rows = c(),
    P = rbind(rep(1, 6), 0, 0)
  )

  # i barely wins in a cycle, having a slightly better loss than j
  el[["events"]][["ik_ji|ijki"]] = list(
    win_conditions = el[["events"]][["ijki"]]$win_conditions,
    tie_condition_rows = c(4), # i's loss to k barely better than j's loss to i
    scaling_factor = 2,  # is it sqrt(number of tallies involved)?
    # P matrix is tricky.
    # a ballot that improves both i relative to k and j relative to i (jik) or makes both worse (kij) does not change the outcome. you only elect j by putting jki or kji.
    P = rbind(c(1,1,1,0,1,0),
              c(0,0,0,1,0,1),
              0),
    adjacent_events = "ji_ik|ijki"
  )

  # reverse the i j comparison
  ji_ik_conditions <- el[["events"]][["ijki"]]$win_conditions
  ji_ik_conditions[4,1:6] <- -ji_ik_conditions[4,1:6]

  el[["events"]][["ji_ik|ijki"]] = list(
    win_conditions = ji_ik_conditions,
    tie_condition_rows = c(4), # i's loss to k barely better than j's loss to i
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
