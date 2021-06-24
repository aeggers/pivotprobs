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

colnames(P3_2) <- c(ballots_ABD, "nul", WTF)

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



P_mat_3_to_4 <-

ordinal_reduction_of <- function(shorter, longer){
  splitted_shorter <- str_split(shorter, "")
  splitted_longer <- str_split(longer, "")
  purge
}
