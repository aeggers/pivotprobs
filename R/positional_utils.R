
positional_win_conditions <- function(score_vector = NULL, cand_names = NULL, n = 5000, sep = ""){
  if(is.null(score_vector)){"You must provide a score_vector."}
  k <- length(score_vector)
  if(is.null(cand_names)){cand_names <- letters[1:k]}

  # get the ballot vector
  df <- as.data.frame(permutations(n = length(cand_names), r = length(cand_names), v = cand_names))
  ballot_vector <- apply(df, 1, paste, collapse = "")

  sv <- normalize_score_vector(score_vector)

  # this chunk should be a function
  split_ballot_list <- ballot_vector %>% str_split(sep)
  cand_names <- split_ballot_list %>% unlist() %>% unique() %>% sort()
  ballot_order_list <- split_ballot_list %>% map(order)
  # some error checking
  ballot_lengths <- lapply(ballot_order_list, length) %>% unlist() %>% unique()
  if(length(ballot_lengths) != 1){stop("Ballots provided are not all of the same length.")}
  if(ballot_lengths != length(sv)){stop("score_vector must be the same length as the ballots.")}
  # back to business
  score_mat <- t(matrix(sv[ballot_order_list %>% unlist()], ncol = length(sv), byrow = T)) # transpose makes it look like P matrix
  rownames(score_mat) <- cand_names
  colnames(score_mat) <- ballot_vector

  cbind(matrix(score_mat[1, ], nrow = nrow(score_mat) - 1, ncol = ncol(score_mat), byrow = T) - score_mat[2:nrow(score_mat), ], rep(1/n, nrow(score_mat) - 1))

}

test <- FALSE
if(test){
  out1 <- positional_event_probabilities(alpha = c(9, 6, 3,4, 5,7), tol = .1)
  P_mat_from_eppp(out1)
  t(c(10, 7, 4)) %*% P_mat_from_eppp(out1)
  out2 <- positional_event_probabilities(alpha = c(9, 6, 3,4, 5,7), skip_non_pivot_events = T, merge_adjacent_pivot_events = T, tol = .1)
  out2a <- positional_event_probabilities(alpha = c(9, 6, 4,6, 3,7), skip_non_pivot_events = T, merge_adjacent_pivot_events = T, drop_dimension = T, tol = .1)
  # these below produce an erroneous facet error
  # look at that output further, and consider rational arithmetic, or ignoring the error
  out3 <- positional_event_probabilities(alpha = sample(2:7, size = 24, replace = T), skip_non_pivot_events = T, merge_adjacent_pivot_events = T, tol = .1)
  # this also produces an erroneous facet error
  out4 <- positional_event_probabilities(alpha = rep(4, 24), skip_non_pivot_events = T, drop_dimension = T, merge_adjacent_pivot_events = T, tol = .1)

}

# this reuses a lot of code from the plurality version, but it was a tradeoff between duplicating code and aiming for abstraction that complicated too much
positional_event_probabilities <- function(score_vector = NULL, skip_non_pivot_events = F, merge_adjacent_pivot_events = F, drop_dimension = F, n = 5000, alpha = NULL, mu = NULL, sigma = NULL, precision = NULL, cand_names = NULL, sep = "", store_time = T, ...){

  if(merge_adjacent_pivot_events & drop_dimension){
    limits <- c(1/(2*n), 1/(2*n))
  } else if(merge_adjacent_pivot_events){
    limits <- c(-1/(2*n), 1/(2*n))
  } else if(drop_dimension){
    limits <- c(0, 0)
  } else{
    limits <- c(0, 1/n)
  }

  time_start <- Sys.time()
  if(!is.null(alpha) | (!is.null(mu) & !is.null(precision)) & is.null(sigma)){
    # this is Dirichlet
    distribution = "dirichlet"
    if(is.null(alpha)){alpha <- mu*precision}
    num_ballots <- length(alpha)
  }else if(!is.null(mu) & !is.null(sigma)){
    # this is logistic normal
    distribution = "logisticnormal"
    num_ballots <- length(mu)
  }else{
    stop("Please pass parameters that allow me to determine the intended distribution over election outcomes.")
  }

  cands_n_ballots <- data.frame(cands = 1:10, ballots = factorial(1:10))
  num_cands <- cands_n_ballots$cands[cands_n_ballots$ballots == num_ballots]
  if(length(num_cands) == 0){stop(paste0("The number of ballot parameters supplied (", num_ballots, ") does not correspond with an integer number of candidates between 1 and 10. Did you not provide a vector with a parameter for all possible orderings of the candidates?"))}

  if(is.null(cand_names)){cand_names <- letters[1:num_cands]}
  if(is.null(score_vector)){
    warning(paste0("No score_vector supplied. Assuming Borda count with ", num_cands, " candidates."))
    score_vector <- seq(num_cands - 1, 0, by = -1)
  }

  # make matrix of inequalities for candidate a winning outright, given number of candidates k and electorate size n
  im <- positional_win_conditions(score_vector = score_vector, cand_names = cand_names, n = n)

  # storage -- one entry per election event
  out <- list()

  # we cycle over candidates
  for(i in 1:num_cands){
    # we will shuffle the parameters and cand_names so that this candidate is first, i.e. candidate a
    these_cand_indices <- c(i,(1:num_cands)[-i]) # in plurality we can use this both for candidates and for alpha/mu/sigma.
    # but in positional methods we need to get the indices for the ballots, which we apply to alpha etc
    these_cand_names <- cand_names[these_cand_indices]
    these_ballot_indices <- ordered_ballot_vector_from_candidate_vector(these_cand_names) %>% names() %>% as.numeric()

    # for now we focus on non pivot events and single ties
    for(m in 0:1){ # could eventually go to (num_cands-1)){
      # you can skip the non-pivot events
      if(skip_non_pivot_events & m == 0){next}
      # we enumerate the indices of the constraints corresponding to the possible candidates who could be in a combination of suze m
      matrix_of_rows_to_alter <- combn(num_cands-1, m) # pointless when m can only be 0 or 1
      for(j in 1:ncol(matrix_of_rows_to_alter)){
        # rta is a vector of those indices, saying which constraint will bind
        rta <- matrix_of_rows_to_alter[,j]

        this_time_start <- Sys.time()
        # we name the pivot event: e.g. "a" means a wins outright no matter what extra ballot is submitted, "a_bc" means b and c are each less than one vote behind
        this_name <- paste0(c(these_cand_names[1], paste0(sort(these_cand_names[1 + rta]), collapse = sep)), collapse = "_")
        if(merge_adjacent_pivot_events){
          # first we check if this pivot event is already covered

          # strip this name
          altered_sep <- sep
          if(sep == ""){altered_sep = "STRING NOONE WOULD USE"}
          this_stripped_name <- this_name %>% str_replace_all("_", "") %>% str_replace_all(altered_sep, "") %>% str_split("") %>% `[[`(1) %>% sort() %>% paste0(collapse = "")
          stored_stripped_names <- names(out) %>% str_replace_all("_", "") %>% str_replace_all(altered_sep, "") %>% str_split("") %>% map(sort) %>% map(str_c, collapse = "") %>% unlist()
          to_grab <- which(stored_stripped_names == this_stripped_name)[1]
          if(!is.na(to_grab)){
            out[[this_name]] <- out[[names(out)[to_grab]]]
            time_diff <- Sys.time() - this_time_start
            units(time_diff) <- "secs"
            out[[this_name]]$seconds_elapsed <- as.double(time_diff)
            next
          }
        }
        # this is the crucial function: getting the S array (the array of simplices over which top integrate) from conditions, with options
        this_S <- S_array_from_inequalities_and_conditions(im, rows_to_alter = rta, limits = limits, drop_dimension = drop_dimension)
        if(distribution == "dirichlet"){
          if(class(this_S) == "matrix" && ncol(this_S) == 1){
            out[[this_name]] <- list(integral = gtools::ddirichlet(as.vector(this_S), alpha[these_ballot_indices])/sqrt(length(alpha)), functionEvaluations = 1, message = "OK")
          }else{
            out[[this_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[these_ballot_indices], ...)
          }
        }else if(distribution == "logisticnormal"){
          if(class(this_S) == "matrix" && ncol(this_S) == 1){
            out[[this_name]] <- list(integral = dlogisticnormal(as.vector(this_S), mu[these_ballot_indices], sigma[these_ballot_indices, these_ballot_indices])/sqrt(length(mu)), functionEvaluations = 1, message = "OK")
          }else{
            out[[this_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[these_ballot_indices], sigma = sigma[these_ballot_indices, these_ballot_indices], ...)
          }
        }
        # if we are integrating at equality rather than near-equality (i.e. integrating on facets rather than in thin subspaces), then we need to thicken out the space. if there were m equalities we multiply by (1/n)^m.
        if(drop_dimension){
          out[[this_name]]$integral <- out[[this_name]]$integral*((1/(sqrt(2)*n))^m)
        }

        ballot_vector <- colnames(im)[-ncol(im)]
        this_W_mat <- positional_P_mat_for_event(event = this_name, ballot_vector = ballot_vector, score_vector = score_vector)
        out[[this_name]]$W_mat <- this_W_mat

        if(store_time){
          time_diff <- Sys.time() - this_time_start
          units(time_diff) <- "secs"
          out[[this_name]]$seconds_elapsed <- as.double(time_diff)
        }

      }
    }
  }
  if(store_time){
    time_diff <- Sys.time() - time_start
    units(time_diff) <- "secs"
    out[["total"]]$seconds_elapsed <- as.double(time_diff)
  }
  out
}

normalize_score_vector <- function(score_vector){
  sv <- sort(score_vector, decreasing = T)
  sv <- sv - sv[length(sv)]
  sv/(sv[1])
}

positional_P_mat_for_event <- function(event = NULL, ballot_vector = NULL, score_vector = NULL, sep = ""){
  # event is an election event with an underscore separating the leader from the others. "a_", "a_b", "a_bc" are valid examples, though for now we ignore higher-order pivot events like a_bc.
  # ballot_vector needs to be in the order used for distribution parameters (e.g. alpha). Usually I use abc, acb, bac, bca, etc.
  # score vector indicates how many points you get for a top ranking, second ranking, etc. we assume this is weakly decreasing, so if you pass (1,2,3) we convert to 3,2,1. No need to normalize as that happens in the function.
  # sep is an optional character separating the trailing candidate(s). must not be "_" as this is reserved for the separation between the leader and follower(s).

  if(is.null(event)){stop("You must provide an event, e.g. 'a_' or 'a_b'.")}
  if(is.null(ballot_vector)){stop("You must provide a vector of ballot names, e.g. c('abc', 'acb', 'bac', ...).")}
  if(is.null(score_vector)){stop("You must provide a score_vector.")}

  # normalize score vector: max is 1, min is 0
  sv <- normalize_score_vector(score_vector)

  # convert ballot_vector (e.g. c("abc", "acb") into matrix saying what score each candidate gets for each ballot
  split_ballot_list <- ballot_vector %>% str_split(sep)
  cand_names <- split_ballot_list %>% unlist() %>% unique() %>% sort()
  ballot_order_list <- split_ballot_list %>% map(order)
  # some error checking
  ballot_lengths <- lapply(ballot_order_list, length) %>% unlist() %>% unique()
  if(length(ballot_lengths) != 1){stop("Ballots provided are not all of the same length.")}
  if(ballot_lengths != length(sv)){stop("score_vector must be the same length as the ballots.")}
  # back to business
  score_mat <- t(matrix(sv[ballot_order_list %>% unlist()], ncol = length(sv), byrow = T)) # transpose makes it look like P matrix
  rownames(score_mat) <- cand_names
  colnames(score_mat) <- ballot_vector

  # make the output
  P <- matrix(0, nrow = length(cand_names), ncol = length(ballot_vector), dimnames = list(cand_names, ballot_vector))

  # now detect what kind of event we have here
  if(!str_detect(event, "_")){stop("We expect an event with an underscore (_).")}
  leader <- event %>% str_extract("^\\w(?=_)")
  followers <- event %>% str_replace(".+?_", "") %>% str_split(sep) %>% unlist()
  if(length(followers) == 0){ # e.g. "a_"
    # this is not a pivot event
    P[leader,] <- 1
  }else if(length(followers) == 1){
    P[leader,] <- 1 + score_mat[leader,] - score_mat[followers,]
    P[P > 1] <- 1
    P[followers,] <- 1 - P[leader,]
  }else{
    warning("We are ignoring higher-order pivot events in making the P matrix.")
  }
  P

}

test <- FALSE
if(test){
  positional_P_mat_for_event(event = "a_b", ballot_vector = c("abc", "acb", "bac", "bca", "cab", "cba"), score_vector = c(3,2,0))
}

# given a vector of candidaates, return the ordered vector of ordered ballots
ordered_ballot_vector_from_candidate_vector <- function(cand_names){
  # data frame of misordered candidate names
  df <- as.data.frame(permutations(n = length(cand_names), r = length(cand_names), v = cand_names))
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
