# I use this function.
# given a vector of candidates, return the ordered vector of ordered ballots
ordered_ballot_vector_from_candidate_vector <- function(cand_names){
  # data frame of misordered candidate names
  df <- as.data.frame(gtools::permutations(n = length(cand_names), r = length(cand_names), v = cand_names))
  for(j in 1:ncol(df)){
    df[[paste0("V", j)]] <- factor(df[[paste0("V", j)]], levels = cand_names)
  }
  odf <- psychTools::dfOrder(df, 1:ncol(df))
  odf %>% apply(1, paste, collapse = "")
}

test <- FALSE
if(test){
  ordered_ballot_vector_from_candidate_vector(c("b", "a", "c"))
  ordered_ballot_vector_from_candidate_vector(c("c", "a", "b"))
  ordered_ballot_vector_from_candidate_vector(c("c", "a", "b", "d", "e"))
}

## note in the archive I have a positional method allowing k candidates
