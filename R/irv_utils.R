condense_mat <- function(mat){
  # replace NA with zero
  mat[is.na(mat)] <- 0
  cns <- colnames(mat)
  ucns <- unique(colnames(mat))
  new_mat <- matrix(0, nrow = nrow(mat), ncol = length(ucns))
  colnames(new_mat) <- ucns
  for(j in 1:ncol(mat)){
    new_mat[,cns[j]] <- new_mat[,cns[j]] + mat[,j]
  }
  new_mat
}

drop_candidate_and_condense_matrix <- function(sims, cand_to_drop = "A"){
  # drop
  colnames(sims) <- str_replace(colnames(sims), cand_to_drop, "")
  # condense
  condense_mat(sims)
}

drop_candidate_and_condense <- function(tb, cand_to_drop = "A"){
  ids <- tb$id

  tb %>%
    select(-id) %>%
    as.matrix() %>%
    drop_candidate_and_condense_matrix(cand_to_drop) %>%
    as_tibble() %>%
    mutate(id = ids) %>%
    relocate(id)
}

rank_mat <- function(frs){
  out <- matrix(1, nrow = nrow(frs), ncol = ncol(frs))
  colnames(out) <- colnames(frs)
  for(j in 1:(ncol(out)-1)){ # j is the column of the candidate we are ranking
    for(k in 1:ncol(frs)){ # k is the column of the candidate we are comparing to
      if(j == k){next}
      out[,j] <- out[,j] + as.integer(frs[,j] < frs[,k])
    }
  }
  out[,ncol(out)] <- sum(1:ncol(out)) - apply(out[,1:(ncol(out) - 1)], 1, sum)
  out
}

# test: much faster than ranking, and gets it right.
# frs <- gtools::rdirichlet(100000, alpha = rep(5, 3))
# system.time(rank_slow <- apply(frs, 1, rank))
# system.time(rank_fast <- rank_mat(frs))


get_loser2 <- function(sims){
  # optimized I think
  loser <- rep("", nrow(sims))
  the_min <- rep(1, nrow(sims))
  for(j in 1:ncol(sims)){
    this_is_lower <- sims[,j] > 0 & sims[,j] < the_min
    the_min[this_is_lower] <-  sims[this_is_lower, j]
    loser[this_is_lower] <- colnames(sims)[j]
  }
  loser
}

get_first_rank_shares <- function(sims){
  colnames(sims) <- str_extract(colnames(sims), "^.")
  condense_mat(sims)
}

get_loser_from_tibble_of_sims <- function(df){
  df %>%
    select(-id) %>%
    as.matrix() %>%
    get_first_rank_shares() %>%
    get_loser2()
}

get_winner_for_each_elimination <- function(sims){

  cands <- names(sims)[2] %>% str_split("") %>% .[[1]] %>% sort()

  tibble(dropped_cand = cands, sims = list(sims)) %>%
    mutate(sims = map2(sims, dropped_cand, drop_candidate_and_condense)) %>%
    mutate(winners = map(sims, irv_winners)) %>%
    select(-sims) %>%
    unnest(cols = c(winners)) %>%
    pivot_wider(names_from = dropped_cand, values_from = winner)

}

parse_event_name <- function(event){
  cands <<- event %>% str_replace("_", "") %>% str_split("") %>% .[[1]]
  w <<- cands[1]; x <<- cands[2]; y <<- cands[3]; z <<- cands[4]
}

round_0_irv_pivot_prob <- function(sims, n = 1000, event = "AB_CD", wfee = NULL, noisy = F){
  # wfee stands for "winner for each elimination"
  # assign candidate names to variables
  parse_event_name(event)
  # get first rank shares and first rank ranks
  frs <- PP_LIBRARY[["frs"]]
  if(is.null(frs)){frs <- get_first_rank_shares(sims %>% select(-id) %>% as.matrix()); PP_LIBRARY[["frs"]] <<- frs}
  ranks <- PP_LIBRARY[["ranks"]]
  if(is.null(ranks)){ranks <- rank_mat(frs); PP_LIBRARY[["ranks"]] <<- ranks}
  # condition 1: w and x in last two places
  cond1 <- ranks[,x] + ranks[,w] == 2*ncol(frs) - 1
  if(noisy){cat("cond1 passed by ", sum(cond1, na.rm = T), " cases.\n")}
  if(is.null(wfee)){
    # who would win given each
    wfee <- get_winner_for_each_elimination(sims)
  }
  cond2 <- wfee[,x] == y
  if(noisy){cat("cond2 passed by ", sum(cond2, na.rm = T), " cases.\n")}
  cond3 <- wfee[,w] == z
  if(noisy){cat("cond3 passed by ", sum(cond3, na.rm = T), " cases.\n")}
  cond <- cond1 & cond2 & cond3
  if(sum(cond, na.rm = T) == 0){return(0)}
  # now compute density
  w_vs_x <- frs[,w] - frs[,x]
  the_density <- try(density_estimate(x = w_vs_x[cond], eval.points = c(0)), silent = T)
  if(class(the_density) == "try-error"){
    msg <- geterrmessage()
    cat("Error: ", msg, "-- setting density to zero.\n")
    the_density <- 0
  }
  mean(cond)*the_density*(1/n)
}

make_P_from_event_name <- function(event, cands = c("A", "B", "C", "D")){
  parse_event_name(event)
  P <- make_empty_P(cands)
  ballots <- colnames(P)
  z_wins <- str_detect(ballots, str_c("^", x))
  P[y, !z_wins] <- 1
  P[z, z_wins] <- 1
  P
}

round_0_pivot_probs <- function(sims, n = 1000, reporting = 1){
  out <- list()
  PP_LIBRARY <<- list() # clearing out by default
  wfee <- get_winner_for_each_elimination(sims)
  cands <- names(sims)[2] %>% str_split("") %>% .[[1]]
  for(w in cands){
    for(x in cands){
      if(w >= x){next}
      for(y in cands){
        if(x == y){next}
        for(z in cands){
          if(y == z | w == z){next}
          event <- str_c(w, x, "_", y, z)
          if(reporting >= 2){cat(event)}
          pp <- sims %>% round_0_irv_pivot_prob(n = n, event = event, wfee = wfee)
          this_P <- make_P_from_event_name(event, cands)
          out[[event]] <- list(integral = pp, P = this_P)
          # mirror event
          event <- str_c(x, w, "_", z, y)
          if(reporting >= 2){cat(" | ", event, "\n")}
          this_P <- make_P_from_event_name(event, cands)
          out[[event]] <- list(integral = pp, P = this_P)
        }
      }
    }
  }
  out
}

PP_LIBRARY <- list()

# initially all events had zero probability because cond1 could not be met -- was operating on frs and not ranks.

parse_1r_event_name <- function(event){
  cands <<- event %>% str_replace("\\.", "") %>% str_replace("_", "") %>% str_split("") %>% .[[1]]
  v <<- cands[1]; w <<- cands[2]; x <<- cands[3]; y <<- cands[4]; z <<- cands[5]
}

round_1_irv_pivot_prob <- function(sims, n = 1000, event = "D.AB_AC", wfee = NULL, noisy = F){
  # wfee stands for "winner for each elimination"
  # assign candidate names to variables
  parse_1r_event_name(event)
  # get first rank shares and first rank ranks
  frs <- PP_LIBRARY[["frs"]]
  if(is.null(frs)){frs <- get_first_rank_shares(sims %>% select(-id) %>% as.matrix()); PP_LIBRARY[["frs"]] <<- frs}
  ranks <- PP_LIBRARY[["ranks"]]
  if(is.null(ranks)){ranks <- rank_mat(frs); PP_LIBRARY[["ranks"]] <<- ranks}
  # condition 1: v is last in 0th round
  cond1 <- ranks[,v] == ncol(frs)
  if(noisy){cat("cond1 passed by ", sum(cond1, na.rm = T), " cases.\n")}
  # now we shift to the next round
  # condition 2: neither w nor x wins the first round
  sims1_v <- PP_LIBRARY[[str_c("sims_drop_", v)]]
  if(is.null(sims1_v)){sims1_v <- drop_candidate_and_condense(sims, cand_to_drop = v); PP_LIBRARY[[str_c("sims_drop_", v)]] <<- sims1_v}
  frs1_v <- PP_LIBRARY[[str_c("frs_drop_", v)]]
  if(is.null(frs1_v)){frs1_v <- get_first_rank_shares(sims1_v %>% select(-id) %>% as.matrix()); PP_LIBRARY[[str_c("frs_drop_", v)]] <<- frs1_v}
  ranks1_v <- PP_LIBRARY[[str_c("ranks_drop_", v)]]
  if(is.null(ranks1_v)){ranks1_v <- rank_mat(frs1_v); PP_LIBRARY[[str_c("ranks_drop_", v)]] <<- ranks1_v}
  cond2 <- ranks1_v[,w] != 1 & ranks1_v[,x] != 1
  if(noisy){cat("cond2 passed by ", sum(cond2, na.rm = T), " cases.\n")}
  if(is.null(wfee)){
    # who would win given each elimination
    wfee <- get_winner_for_each_elimination(sims1)
  }
  # cond3: eliminate x and get y
  cond3 <- wfee[,x] == y
  if(noisy){cat("cond3 passed by ", sum(cond3, na.rm = T), " cases.\n")}
  # cond4: eliminate w and get z
  cond4 <- wfee[,w] == z
  if(noisy){cat("cond4 passed by ", sum(cond4, na.rm = T), " cases.\n")}
  cond <- cond1 & cond2 & cond3 & cond4
  if(sum(cond, na.rm = T) == 0){return(0)}
  # now compute density
  w_vs_x <- frs1_v[,w] - frs1_v[,x]
  the_density <- try(density_estimate(x = w_vs_x[cond], eval.points = c(0)), silent = T)
  if(class(the_density) == "try-error"){
    msg <- geterrmessage()
    cat("Error: ", msg, "-- setting density to zero.\n")
    the_density <- 0
  }
  mean(cond)*the_density*(1/n)
}

# this could be expanded to deal with multiple previous rounds
round_1_pivot_probs <- function(sims, n = 1000, reporting = 1){
  out <- list()
  PP_LIBRARY <<- list() # clear out by default -- slightly inefficient
  cands <- names(sims)[2] %>% str_split("") %>% .[[1]]
  # v.wx_yz restrictions: v is not w, x, y, or z. w is not z, x is not y, y is not z.
  for(v in cands){
    drop_v_sims <- drop_candidate_and_condense(sims, cand_to_drop = v)
    wfee <- get_winner_for_each_elimination(drop_v_sims)
    for(w in cands){
      for(x in cands){
        if(w >= x){next}
        for(y in cands){
          if(x == y){next}
          for(z in cands){
            if(y == z | w == z | v %in% c(w, x, y, z)){next}
            event <- str_c(v, ".", w, x, "_", y, z)
            if(reporting >= 2){cat(event)}
            pp <- sims %>% round_1_irv_pivot_prob(n = n, event = event, wfee = wfee)
            this_P <- make_P_from_1r_event_name(event, cands)
            out[[event]] <- list(integral = pp, P = this_P)
            # mirror event
            event <- str_c(v, ".", x, w, "_", z, y)
            if(reporting >= 2){cat(" | ", event, "\n")}
            this_P <- make_P_from_1r_event_name(event, cands)
            out[[event]] <- list(integral = pp, P = this_P)
          }
        }
      }
    }
  }
  out
}

make_empty_P <- function(cands){
  k <- length(cands)
  P <- matrix(0, nrow = k, ncol = factorial(k))
  rownames(P) <- cands
  ballots <- gtools::permutations(n = k, r = k, v = cands) %>%
    apply(1, paste, collapse = "")
  colnames(P) <- ballots
  P
}

# v.wx_yz P: y wins unless you rank x first or v first and x second.
make_P_from_1r_event_name <- function(event, cands = c("A", "B", "C", "D")){
  parse_1r_event_name(event)
  P <- make_empty_P(cands)
  ballots <- colnames(P)
  z_wins <- str_detect(ballots, str_c("^", x)) | str_detect(ballots, str_c("^", v, x))
  P[y, !z_wins] <- 1
  P[z, z_wins] <- 1
  P
}

# last round pivot events: we just crank until we have only two
get_to_last_round <- function(sims){
  while(str_length(colnames(sims)[2]) > 2){
    sims <- one_round_of_irv(sims)
  }
  sims
}

make_P_from_last_round_event_name <- function(event, cands = c("A", "B", "C", "D")){
  event_split <- event %>% str_split("") %>% .[[1]]
  y <- event_split[1]; z <- event_split[2]

  P <- make_empty_P(cands)
  ballots <- colnames(P)
  y_before_z <- str_detect(ballots, str_c(y, ".*", z))
  P[y, y_before_z] <- 1
  P[z, !y_before_z] <- 1
  P
}

last_round_pivot_probs <- function(sims, n = 1000, reporting = 1){
  cands <- names(sims)[2] %>% str_split("") %>% .[[1]]
  sims_last <- get_to_last_round(sims)
  sims_last %>%
    pivot_longer(cols = -id) %>%
    filter(!is.na(value) & value > 0) %>%
    group_by(id) %>%
    slice(1) %>%
    separate(col = name, into = c("y", "z"), sep = 1) %>%
    mutate(y_vs_z = 2*(value - .5)) -> margins
  out <- list()
  for(y in cands){
    for(z in cands[cands > y]){
      event <- str_c(y, z)
      if(reporting >= 2){cat(event)}
      cond <- margins$y == y & margins$z == z
      the_density <- density_estimate(x = margins$y_vs_z[cond], eval.points = c(0))
      pp <- mean(cond)*the_density*(1/n)
      this_P <- make_P_from_last_round_event_name(event, cands)
      out[[event]] <- list(integral = pp, P = this_P)
      event <- str_c(z, y)
      if(reporting >= 2){cat(" | ", event, "\n")}
      this_P <- make_P_from_last_round_event_name(event, cands)
      out[[event]] <- list(integral = pp, P = this_P)
    }
  }
  out
}
