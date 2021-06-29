## we will have a 3x6 P mat.
## we want to expand this into a 3x24 P mat, noting equivalencies.

# this is not that hard

ballots_ABC <- gtools::permutations(n = 3, r = 3, v = c("A", "B", "C")) %>%
  apply(1, paste, collapse = "")

ballots_ABD <- gtools::permutations(n = 3, r = 3, v = c("A", "B", "D")) %>%
  apply(1, paste, collapse = "")

ballots4 <- gtools::permutations(n = 4, r = 4, v = c("A", "B", "C", "D")) %>%
  apply(1, paste, collapse = "")

irv_election(n = 1000, s = 0) %>%
  election_event_probs(method = "mc", alpha = c(10,3, 4,5,2,7)) %>%
  combine_P_matrices() -> P3_1

colnames(P3_1) <- c(ballots_ABC, "nul")

irv_election(n = 1000, s = 0) %>%
  election_event_probs(method = "mc", alpha = c(11,4, 3,4,5,2)) %>%
  combine_P_matrices() -> P3_2

colnames(P3_2) <- c(ballots_ABD, "nul")

P4 <- matrix(0, nrow = 3, ncol = length(ballots4))
colnames(P4) <- ballots4
rownames(P4) <- c("A", "B", "C")

# this bit of code applies the smaller P3 to the appropriate places.
# we need to add one layer: loop across the different P3s, and weight appropriately.
# but wait: are these P matrices what we want to be combining? I think so. at least if we ignore first round pivot events
for(j in 1:6){
  this_colname <- colnames(P3_1)[j]
  target_cols <- str_replace(ballots4, "D", "") == this_colname
  P4[ , target_cols] <- P4[ , target_cols] + P3_1[,j]
}



# really much better if we don't ignore first round pivot events.
# could we look at them as independent of what happens afterward?
# like, get the utility given 3 candidates running (via simulation?)
# and then compute the probability of being able to decide between
# them (plurality). but how then to combine with the rest?
# well, let's extend this approach. could compute probability of each candidate winning given each pair in the last round. then in first round you choose who to put first by saying, "what is probability of A and B tied for second? and if they do, and I rank A first, then I get the AC pairing; if I vote B, I get the BC pairing."
# I am not thinking very well about all of this. I guess need to re-read some of my own stuff on this, get back in the mindset.

# if I can write a function that gives me the 0th round pivot probability for a set of candidates, given a set of simulated results, and gives me a P matrix given that set of candidates, then I can just fill up a list according to rules about how many we need. that's definitely the ticket.

# remembering that in pivotprobs we used various $s$'s. do we start without that?

random_alpha <- sample(1:12, size = 24, replace = T)
sims <- gtools::rdirichlet(100, alpha = random_alpha)
colnames(sims) <- gtools::permutations(n = 4, r = 4, v = c("A", "B", "C", "D")) %>%
  apply(1, paste, collapse = "")

left0 <- "A"
right0 <- "C"
win_left0 <- "B"
win_right0 <- "D"

# if left0 wins, eliminate right0, reallocate
# then see who has lowest, eliminate that one, reallocate

# so basically we're saying figure out winner

sims <- gtools::rdirichlet(10, alpha = sample(1:12, 6, T))
colnames(sims) <- gtools::permutations(n = 3, r = 3, v = c("A", "B", "C")) %>% apply(1, paste, collapse = "")

condense_mat <- function(mat){
  cns <- unique(colnames(mat))
  new_mat <- matrix(nrow = nrow(mat), ncol = length(cns))
  colnames(new_mat) <- cns
  for(j in 1:ncol(new_mat)){
    new_mat[,j] <- apply(mat[,which(colnames(mat) == cns[j]), drop = F], 1, sum)
  }
  new_mat
}

drop_candidate_and_condense <- function(sims, cand_to_drop = "A"){
  # drop
  colnames(sims) <- str_replace(colnames(sims), cand_to_drop, "")
  # condense
  condense_mat(sims)
}

get_first_rank_shares <- function(sims){
  colnames(sims) <- str_extract(colnames(sims), "^.")
  condense_mat(sims)
}

# a tidier way to do this -- possibly too slow?
drop_loser_and_collapse <- function(result){
  result %>%
    mutate(name = str_extract(name, "^.")) %>%
    group_by(name) %>%
    summarize(vs = sum(vs)) %>%
    ungroup() %>%
    filter(vs == min(vs)) %>%
    pull(name) -> loser

  result %>%
    mutate(name = str_replace(name, loser, "")) %>%
    group_by(name) %>%
    summarize(vs = sum(vs)) %>%
    ungroup()
}

get_result_row_by_row <- function(sims){

  # check that we have colnames
  stopifnot(length(unique(colnames(sims))) == ncol(sims))

  # make into nested long sims
  sims %>%
    as_tibble() %>%
    mutate(elec = 1:n()) %>%
    pivot_longer(cols = c(-elec), values_to = "vs") %>%
    group_by(elec) %>%
    nest() %>%
    rename(result = data) -> nls

  while(nls$result[[1]] %>% nrow() > 1){
    nls %>%
      mutate(result = map(result, drop_loser_and_collapse)) -> nls
  }

  nls %>%
    unnest(cols = c(result)) %>%
    select(elec, winner = name)
}

sims %>% get_result_row_by_row()

# so what are we looking to do? for 0th round pivot events, we want to:
# drop a candidate for all rows (possibly after having adjusted the shares so there is a tie), see who wins in each row.
# to make this faster, could divide based on who gets eliminated first?
# so, the first elimination is manual. second is based on shares.
# what I do in the existing code is: compute score for each in first round based on fixed ordering of v-vec, and pairwise margins.
# the thing that's slow is applying drop_loser_and_collapse to each group.
# what about: we have the sims; we ask who the loser is for each row; we divide according to the loser; then we drop the loser and collapse; repeat.

get_loser <- function(sims){
  colnames(sims)[apply(sims, 1, which.min)]
}

sims %>%
  get_first_rank_shares() %>%
  get_loser() -> losers

# tried with by, thinking about manually
out <- list()
for(loser in unique(losers)){
  out[[loser]] <- drop_candidate_and_condense(sims[losers == loser,], loser)
}

out[["C"]] %>%
  get_first_rank_shares() %>%
  get_loser() -> losers2

out2 <- list()
for(loser in unique(losers2)){
  out2[[loser]] <- drop_candidate_and_condense(out[["C"]][losers2 == loser,], loser)
}

# OK we are close to some nice use of recursion here. Need to think about the problem a bit more.
# what about we start with a tibble of sims.
# we get first rank shares to add a column of first round losers.
# we group_by the losers, drop the loser in each nested df, condense. so now we can either return this to a tibble of sims or proceed on the grouped sims.

drop_candidate_and_condense_matrix <- function(sims, cand_to_drop = "A"){
  # drop
  colnames(sims) <- str_replace(colnames(sims), cand_to_drop, "")
  # condense
  condense_mat(sims)
}

drop_candidate_and_condense <- function(tb, cand_to_drop = "A"){
  as.matrix(tb) %>%
    drop_candidate_and_condense_matrix(cand_to_drop) %>%
    as_tibble()
}

which.non.zero.min <- function(vec){
  which(vec > 0 & vec == min(vec[!vec == 0]))
}

get_loser <- function(sims){
  colnames(sims)[apply(sims, 1, which.non.zero.min)]
}

get_first_rank_shares <- function(sims){
  colnames(sims) <- str_extract(colnames(sims), "^.")
  condense_mat(sims)
}

get_loser_from_tibble_of_sims <- function(df){
  as.matrix(df) %>%
    get_first_rank_shares() %>%
    get_loser()
}

one_round_of_irv <- function(sims){
  as_tibble(sims) %>%
    mutate(loser = get_loser_from_tibble_of_sims(.)) %>%
    group_by(loser) %>%
    nest() %>%
    mutate(reduced = map(data, drop_candidate_and_condense, loser)) %>%
    mutate(loser = map(reduced, get_loser_from_tibble_of_sims)) %>%
    ungroup() %>%
    select(reduced) %>%
    unnest(cols = c(reduced)) %>%
    replace(is.na(.), 0) -> out

  out %>%
    summarize(across(where(is.double), sum)) %>%
    pivot_longer(cols = where(is.double)) %>%
    filter(value > 0) %>%
    pull(name) %>%
    sort() -> non_zero_cols

  out %>% select(all_of(non_zero_cols))
}


get_irv_winners_fast <- function(sims){
  stopifnot(length(unique(colnames(sims))) == ncol(sims))
  while(max(as.matrix(sims)[1,]) < 1){
    sims <- one_round_of_irv(sims)
  }
}

random_alpha <- sample(1:12, size = 24, replace = T)
sims <- gtools::rdirichlet(100000, alpha = random_alpha)
colnames(sims) <- gtools::permutations(n = 4, r = 4, v = c("A", "B", "C", "D")) %>%
  apply(1, paste, collapse = "")

winners <- get_irv_winners_fast(sims)





sims <- gtools::rdirichlet(10, alpha = sample(1:12, 6, T))
colnames(sims) <- gtools::permutations(n = 3, r = 3, v = c("A", "B", "C")) %>% apply(1, paste, collapse = "")



sims %>%
  one_round_of_irv() %>%
  one_round_of_irv() -> x

%>%
  one_round_of_irv() %>%
  one_round_of_irv()


# so we only need as many operations as there are candidates in each round.
# but to do that we need to bring it back together


#
# we take a matrix of sims.





nested %>%




nested3 %>%
  select(loser = loser2, result = reduced) %>%
  group_by(loser) %>%
  nest() -> nested4
# this is promising, come back to it.

# the idea is to reduce the number of operations -- not one per row, but one per loser in each round of counting. I think we are VERY close.

# maybe a function on a grouped dataset?
take the sims, get the first rank shares, get the loser,


# OK this is promising, come back to it.

A_lost <- sims[which(losers == "A"),]
A_lost %>%
  drop_candidate_and_condense("A") %>%
  get_first_rank_shares() %>%
  get_loser() -> losers2

by(sims, losers, drop_candidate_and_condense, losers)


random_alpha <- sample(1:12, size = 24, replace = T)
sims <- gtools::rdirichlet(1000, alpha = random_alpha)
colnames(sims) <- gtools::permutations(n = 4, r = 4, v = c("A", "B", "C", "D")) %>%
  apply(1, paste, collapse = "")

sims %>%
  drop_candidate_and_condense("A") %>%
  get_result_row_by_row() -> result_without_A
# okay that takes way too long.



nls %>%
  mutate(result = map(result, drop_loser_and_collapse)) -> nls2

get_row_loser <- function(sims){

}

do_irv_round <- function(sims){
  frs <- get_first_rank_shares(sims)

  colnames(sims) <- str_extract(colnames(sims), "^.")
  condense_mat(sims)
}
