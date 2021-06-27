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
    mutate(id = ids)
}

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

one_round_of_irv <- function(sims){
  sims %>%
    mutate(loser = get_loser_from_tibble_of_sims(.)) %>%
    # now we group so that we only have to drop and condense once per losing candidate (rather than once per row)
    group_by(loser) %>%
    nest() %>%
    mutate(reduced = map(data, drop_candidate_and_condense, loser)) %>%
    ungroup() %>%
    select(reduced) %>%
    unnest(cols = c(reduced)) -> out

   out %>%
     colSums(na.rm = T) -> col_sums

  col_sums[col_sums > 0] %>%
    names() %>%
    sort() -> non_zero_cols

  out %>%
    select(all_of(non_zero_cols)) %>%
    arrange(id) %>%
    relocate(id)
}


get_irv_winners_fast <- function(sims){
  stopifnot(length(unique(colnames(sims))) == ncol(sims))
  while(!(all.equal(max(as.matrix(sims)[1,-1], na.rm = T), 1)) == T){
    sims <- one_round_of_irv(sims)
  }
  sims %>%
    pivot_longer(cols = -id, names_to = "winner", values_to = "won") %>%
    filter(near(won, 1)) %>%
    select(id, winner)
}

make_irv_sims <- function(k = 4, n = 10000, alpha = NULL){
  if(is.null(alpha)){alpha <- sample(1:12, size = factorial(k), replace = T)}
  sims <- gtools::rdirichlet(n, alpha = alpha)
  colnames(sims) <- gtools::permutations(n = k, r = k, v = LETTERS[1:k]) %>%
    apply(1, paste, collapse = "")
  sims %>%
    as_tibble() %>%
    mutate(id = 1:n()) %>%
    relocate(id)
}

# sims <- make_irv_sims(k = 4, n = 100000)
# system.time(winners <- get_irv_winners_fast(sims))

sims <- make_irv_sims(k = 4, n = 100000)
system.time(winners <- get_irv_winners_fast(sims))

sims <- make_irv_sims(k = 3, n = 100000)
system.time(winners <- get_irv_winners_fast(sims))

# so now the idea must be to
# do we have to carry around an id so that we know which race ends which way?
# the idea is: figure out cases where e.g. B would win if A were eliminated but D would win if C were eliminated, and A and C would tie for last if you averaged their scores, and then get the density.
# so yes I need to keep straight the ids.
# where will that be an issue?

# is there a general principle for all pivot events?
# -- to get event AB.AC in round 1, you need D eliminated in round 0, and A and B tied in round 1, and A wins if B eliminated and C wins if B eliminated.
# -- AB.AB. in round 1 similar but could happen with C or D eliminated in round 0
# -- to get event AB.CD in round 0, you need A and B tied in round 0, C wins if B eliminated, D wins if A eliminated.
# -- to get event AB in round 2, you need A and B in round 2, and A and B tied.
