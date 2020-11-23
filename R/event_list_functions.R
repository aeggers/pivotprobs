#' Create election list for analysis
#'
#' An election list represents an election
#' in a minimal way: it contains the electorate size, the ballot
#' format, and a list of "election events", i.e. situations in which
#' a single ballot can affect the outcome in a given way.
#' It also contains a label for the voting system, e.g. "positional".
#' An election list is passed to the
#' \code{election_event_probs()}
#' function to compute the probability of each election event given
#' a model of voting outcomes.
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
#' The election methods provided are:
#' \itemize{
#' \item \code{plurality_election()}, where each voter votes for
#' one of \code{k>=3} candidates and the candidate with the most votes wins.
#' \item \code{positional_election()}, where each voter ranks the
#' three candidates and the candidate with the highest score wins, with
#' each top ranking worth 1 point, each bottom ranking worth 0 points,
#' and each middle ranking worth \code{s} points.
#' For now we can only handle 3 candidates.
#' (I have written code to handle more, but the QHull
#' functions seem to have trouble with the geometry.)
#' \item \code{irv_election()} (instant-runoff), where each voter
#' ranks the three candidates and the winner is determined by a process
#' of elimination. The ballots are initially tallied using a
#' positional method; the most commonly used option is plurality,
#' i.e. \code{s=0}.) The lowest-scoring candidate is eliminated.
#' Given three candidates, the winner is then determined according
#' to which of the remaining candidates is ranked higher on more
#' ballots.
#' \item \code{kemeny_young_election()} (maximin), where any
#' candidate who defeats all others in a pairwise comparison (i.e. is ranked higher on more ballots)
#' is the winner and, if there is no such candidate, the winner is the
#' candidate who is defeated by the smallest margin.
#' }
#'
#' @examples
#' plurality_election(k = 3)
#' positional_election(s = .5) # borda count
#' irv_election()
#' kemeny_young_election()
#'
#' @param n The size of the electorate. This (slightly) affects the
#' conditions for election events in the election list
#' (e.g. the region in which one
#' candidate is more than one vote ahead of another).
#' \code{n} is also used when computing pivot event probabilities
#' from an election list in \code{election_event_probs()},
#' where higher \code{n} means lower probability.
#' @param k The number of candidates in plurality.
#' @param max_pivot_event_degree If 1, we consider two-way (near) ties
#' in plurality. If 2, we also consider three-way (near) ties.
#' @param s The value of a second-ranking in a three-candidate
#' positional election, where a top ranking is worth 1 and a bottom ranking is worth 0.
#'
#' @name election_list_functions
NULL

#' @rdname election_list_functions
#' @export
plurality_election <- function(n = 1000, k = 4, max_pivot_event_degree = 1){

  if(!max_pivot_event_degree %in% 1:2){
    stop("max_pivot_event_degree needs to be either 1 or 2. The event_probabilities_from_event_list() function only permutes three candidates, so we cannot handle pivot events of degree 3 or higher, e.g. four-way ties.")
  }

  if(k < 3){stop("There must be at least three candidates in a plurality election.")}

  el <- list()

  el$n <- n
  el$ordinal <- F
  el$system <- "plurality"
  el$k <- k

  el$events <- list()

  win_P <- matrix(0, nrow = k, ncol = k+1)
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
pr_election <- function(n = 1000, M = 3){
  # this one works differently
  if(M < 1){stop("M must be an integer 1 or greater")}
  el <- list()
  el$n <- n
  el$M <- M
  el$ordinal <- F
  el$system <- "dhondt"

  el$events <- list()

  dgp <- dhondt_gridpoints(M = M)

  event_names <- dgp %>% pull(event) %>% unique()

  adjacent_event_names <- event_names %>% str_split("_") %>% lapply(function(x){paste0(x[2], "_", x[1])}) %>% unlist()

  for(i in 1:length(event_names)){
    event_name <- event_names[i]
    endpoints <- dgp %>% filter(event == event_name) %>% select(a, b, c)

    el$events[[event_name]] = list(
      endpoints = endpoints,
      P = P_matrix_from_PR_event_name_and_M(event_name, M),
      adjacent_events = adjacent_event_names[i]
    )
    el$events[[adjacent_event_names[i]]] = list(
      endpoints = endpoints,
      P = P_matrix_from_PR_event_name_and_M(event_name, M, reverse = F),
      adjacent_events = event_name
    )

  }

  el

}

#' @rdname election_list_functions
#' @export
positional_election <- function(n = 1000, s = .5){

  if(s > 1 | s < 0){stop("s must be between 0 and 1.")}

  el <- list()

  el$n <- n
  el$ordinal <- T
  el$system <- "positional"
  el$s <- s

  el$events <- list()

  el[["events"]][["i_"]] <- list(
    win_conditions = rbind(c(1-s, 1, s-1, -1, s, -s, 1/n),  # i scores above j
                     c(1, 1-s, s, -s, s-1, -1, 1/n)),   # i scores above k
    tie_condition_rows = c(),
    P = rbind(rep(1, 7), 0, 0, 0)
    )

  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i_"]]$win_conditions,
    tie_condition_rows = c(1),                               # i barely scores above j
    scaling_factor = 2*sqrt(1 - s + s^2),
    P = rbind(c(1,1,s,0,1,1-s,1),
            c(0,0,1-s,1,0,s,0),
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
  el$system <- "irv"
  el$s <- s

  el$events <- list()

  # positional scaling factor
  psf <- 2*sqrt(1 - s + s^2)

  # i wins over j in the second round
  el[["events"]][["i__j"]] <- list(
    win_conditions = rbind(c(1, 1-s,s,-s,s-1,-1, 1/n),  # i beats k in first round
                       c(s,-s,1,1-s,-1,s-1,1/n),    # so does j
                       c(1,1,-1,-1,1,-1, 1/n)), # i beats j in second round
    tie_condition_rows = c(),
    P = rbind(rep(1, 7), 0, 0),
    adjacent_events = "j__i"
  )

  # second round pivot event
  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i__j"]]$win_conditions,
    tie_condition_rows = c(3),                        # i barely beats j in second round
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0,1), c(0,0,1,1,0,1,0), 0),
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
    P = rbind(c(1,1,s,0,1,1-s,1), c(0,0,1-s,1,0,s,0), 0),
    adjacent_events = "j_i|ji")

  i_j_ik_conditions <- el[["events"]][["i_j|ij"]]$win_conditions
  i_j_ik_conditions[5,1:6] <- -i_j_ik_conditions[5,1:6] # reverse the j k ordering

  el[["events"]][["i_j|ik"]] <- list(
    win_conditions = i_j_ik_conditions,
    tie_condition_rows = c(3),
    scaling_factor = psf,
    P = rbind(c(1,1,s,0,1,1-s,1), 0, c(0,0,1-s,1,0,s,0)),
    adjacent_events = "j_i|ki")

  i_j_kj_conditions <- el[["events"]][["i_j|ij"]]$win_conditions
  i_j_kj_conditions[4,1:6] <- -i_j_kj_conditions[4,1:6]  # reverse the i k ordering

  el[["events"]][["i_j|kj"]] <- list(
    win_conditions = i_j_kj_conditions,
    tie_condition_rows = c(3),
    scaling_factor = psf,
    P = rbind(0, c(0,0,1-s,1,0,s,0), c(1,1,s,0,1,1-s,1)),
    adjacent_events = "j_i|jk")

  el
}

#' @rdname election_list_functions
#' @export
kemeny_young_election <- function(n = 1000){

  el <- list()

  el$n <- n
  el$ordinal <- T
  el$system <- "kemeny_young"

  el$events <- list()

  # i is the condorcet winner by more than one vote over j, with k the condorcet loser.
  el[["events"]][["i__j"]] <- list(
    win_conditions = rbind(c(1, 1, -1, -1, 1, -1, 1/n),   # i pairwise beats j
                       c(1, 1, 1, -1, -1, -1, 1/n),   # i pairwise beats k
                       c(1, -1, 1, 1, -1, -1, 1/n)),  # j pairwise beats k
    tie_condition_rows = c(),
    P = rbind(rep(1, 7), 0, 0)
  )

  # i is the condorcet winner by less than one vote over j
  el[["events"]][["i_j"]] <- list(
    win_conditions = el[["events"]][["i__j"]]$win_conditions,
    tie_condition_rows = c(1),  # i barely beats j
    scaling_factor = sqrt(6),
    P = rbind(c(1,1,0,0,1,0,1), c(0,0,1,1,0,1,0), 0),
    adjacent_events = "j_i"
  )

  # i wins in a cycle where i > j > k > i, with j next and k last
  el[["events"]][["ik__ji|ijki"]] = list(
    win_conditions = rbind(c(1,1,-1,-1,1,-1, 1/n),    # i pairwise beats j
                       c(1,-1,1,1,-1,-1, 1/n),    # j pairwise beats k
                       c(-1, -1, -1, 1, 1, 1, 1/n), # k pairwise beats i
                       c(1,1,0,-1,0,-1,1/n),      # i's loss to k better than j's loss to i
                       c(0,-1, 1, 1,-1, 0, 1/n)   # j's loss to i better than k's loss to j
    ),
    tie_condition_rows = c(),
    P = rbind(rep(1, 7), 0, 0)
  )

  # i barely wins in a cycle,
  # with i's loss to k slightly better than j's loss to i
  el[["events"]][["ik_ji|ijki"]] = list(
    win_conditions = el[["events"]][["ik__ji|ijki"]]$win_conditions,
    tie_condition_rows = c(4), # i's loss to k barely better than j's loss to i
    scaling_factor = 2,  # sqrt(number of tallies involved)
    # a ballot that improves both i relative to k and j relative to i (jik) or makes both worse (kij) does not change the outcome. you only elect j by putting jki or kji.
    P = rbind(c(1,1,1,0,1,0,1),
              c(0,0,0,1,0,1,0),
              0),
    adjacent_events = "ji_ik|ijki"  # j's loss to i slightly better than i's loss to k
  )

  # here is the adjacent pivot event:
  # j's loss to i slightly better than i's loss to k
  ji_ik_conditions <- el[["events"]][["ik__ji|ijki"]]$win_conditions
  ji_ik_conditions[4,1:6] <- -ji_ik_conditions[4,1:6]
  ji_ik_conditions[5,1:6] <- c(1,0,1,0,-1,-1) # i's loss to k better than k's loss to j

  el[["events"]][["ji_ik|ijki"]] = list(
    win_conditions = ji_ik_conditions,
    tie_condition_rows = c(4), # j's loss to i barely better than i's loss to k
    scaling_factor = 2,  # sqrt(number of tallies involved)
    # you only elect i by putting i above k and j below i, so ikj or ijk
    P = rbind(c(1,1,0,0,0,0,0),
              c(0,0,1,1,1,1,1),
              0),
    adjacent_events = "ik_ji|ijki"
  )

  el

}
