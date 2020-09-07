# I think this is the way forward.

event_probabilities_from_event_list <- function(event_list = NULL, alpha = NULL, mu = NULL, precision = NULL, sigma = NULL, sims = NULL, draw_sims = F, num_sims = 100000, window_factor = 100, cand_names = NULL, ordinal = T, drop_dimension = F, merge_adjacent_pivot_events = F, skip_non_pivot_events = F, skip_compound_pivot_events = F, store_time = T, maxEvals = 100000, tol = .01, ...){
  # window_factor is a multiplier -- considered a tie if share is in the window (-(1/n2)*window_factor, (1/n2)*window_factor)

  # start the clock
  time_start <- Sys.time()

  if(is.null(event_list)){stop("You have to pass an event_list.")}

  # storage for output
  out <- list()

  if(is.null(sims)){
    # figure out the distribution from the parameters
    if(!is.null(alpha) | (!is.null(mu) & !is.null(precision)) & is.null(sigma)){
      # this is Dirichlet
      distribution <- "dirichlet"
      if(is.null(alpha)){alpha <- mu*precision}
      mu <- alpha/sum(alpha)
    }else if(!is.null(mu) & !is.null(sigma)){
      # this is logistic normal
      distribution <- "logisticnormal"
    }else{
      stop("Please pass parameters that allow me to determine the intended distribution over election outcomes.")
    }
  }else{
    distribution <- "sims"
  }

  # a label for the situation where we are going to count cases rather than integrate a function
  from_sims <- (draw_sims | !is.null(sims))

  # in that case, we need to merge_adjacent_pivot_events and drop_dimension
  if(from_sims){
    cat("This is simulation-based!\n")
    merge_adjacent_pivot_events <- T
    drop_dimension <- T
  }

  if(draw_sims & is.null(sims)){
    if(!is.numeric(num_sims)){
      stop("If draw_sims is TRUE, we need a numeric num_sims.")
    }
    if(!num_sims > 0){
      stop("num_sims needs to be a positive number.")
    }
    if(num_sims < 10000){
      warning(paste0("num_sims is only ", num_sims, "; probably want more."))
    }
    # now, generate the simulations
    if(distribution == "dirichlet"){
      sims <- gtools::rdirichlet(n = num_sims, alpha = alpha)
    }else if(distribution == "logisticnofmal"){
      sims <- rlogisticnormal(n = num_sims, mu = mu, sigma = sigma)
    }
  }

  # for the limits, we want the value of one vote in terms of vote share, i.e. 1/n. We get this from the conditions supplied in the event_list.
  conditions_1 <- event_list[[1]]$conditions
  base_limit <- conditions_1[1, ncol(conditions_1)]
  if(class(base_limit) != "numeric"){
    stop(paste0("Was expecting last column of first row of first condition to be 1/n; instead it is ", base_limit, "."))
  }
  if(base_limit <= 0){
    stop("I assumed that the last column of the first condition on the first event is 1/n, but it is", base_limit, ".")
  }

  # set limits (width of integration region) based on arguments
  if(from_sims){
    limits <- window_factor*c(-base_limit/2, base_limit/2)
    cat("Base limit was ", base_limit, ".\nLimits are ", paste(limits, collapse = ", "), ".\n", sep = "")
  }else if(merge_adjacent_pivot_events & drop_dimension){
    limits <- c(0, 0) # only first limit used -- a line at the boundary
  } else if(merge_adjacent_pivot_events){
    limits <- c(-base_limit/2, base_limit/2)  # a thin strip straddling the boundary
  } else if(drop_dimension){
    limits <- rep(base_limit/2, 2)   # only first limit used -- a line just short of the boundary
  } else{
    limits <- c(0, base_limit)     # a thin strip short of the boundary
  }

  # check number of ballots/candidates, fill in cand_names if necessary
  if(ordinal){
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
    if(length(mu) != 6){
      stop("For ordinal methods I expect parameters for 6 ballot types. You provided parameters indicating more or fewer than 6 ballot types. Maybe you meant ordinal = F?")
    }
    if(length(cand_names) != 3){
      stop("For ordinal methods I expect 3 candidates. You provided more or fewer than that.")
    }
  }else{  # single-ballot syste
    if(is.null(cand_names)){
      cand_names <- letters[1:length(mu)]
    }
    if(length(cand_names) != length(mu)){
      stop("The length of the parameter vector you provided (mu or alpha) doesn't match the length of cand_names.")
    }
  }

  # get the possible permutations of the candidates
  candidate_permutation_mat <- gtools::permutations(length(cand_names), length(cand_names), cand_names)
  colnames(candidate_permutation_mat) <- letters[9:(9 + length(cand_names) - 1)] # label the columns i,j,k,...

  # storage so we don't make the same S_array more than once
  # the S_array is the key input to SimplicialCubature::adaptIntegrateSimplex.
  S_list <- list()

  # now we permute through possible orderings of candidates
  for(p_row in 1:nrow(candidate_permutation_mat)){
    cand_vector <- candidate_permutation_mat[p_row,]  # e.g. c("c", "a", "b")
    cand_order <- rank(cand_vector)                   # (3, 1, 2)
    if(ordinal){
      ballot_order <- ordered_ballot_vector_from_candidate_vector(cand_vector) %>% names() %>% as.numeric()                                      # (5, 6, 2, 1, 4, 3)
    }else{
      ballot_order <- cand_order                      # (3, 1, 2)
    }
    for(el_index in 1:length(names(event_list))){

      # start the clock on this election event
      this_time_start <- Sys.time()

      this_P <- event_list[[names(event_list)[el_index]]]$P

      # we use the P matrix to determine if this is a pivot event
      non_pivot_event <- this_P %>% apply(1, mean) %>% max() == 1
      if(skip_non_pivot_events & non_pivot_event){next}

      # we can specify skip_compound_pivot_events, i.e. those with more than one "row to alter"; we MUST skip compound pivot events when drop_dimension is TRUE.
      if(skip_compound_pivot_events | (drop_dimension & length(event_list[[names(event_list)[el_index]]]$rows_to_alter) > 1)){next}

      # we need to convert e.g. i_j to c_a
      generic_event_name <- names(event_list)[el_index]

      # a function I use twice here. uses local variables in definition.
      turn_generic_to_specific <- function(generic_event_name){
        generic_event_name %>%
          str_replace_all("i", candidate_permutation_mat[p_row, "i"]) %>%
          str_replace_all("j", candidate_permutation_mat[p_row,"j"]) %>%
          str_replace_all("k", candidate_permutation_mat[p_row, "k"])
      }

      # turn generic event name (e.g. i_j) into name that corresponds to this permutation (e.g. c_a)
      specific_event_name <- generic_event_name %>%
        turn_generic_to_specific() %>%    # i_j => c_a
        regularize_candidate_order_after_underscore()  # c_ba => c_ab, but c_ba|cb stays same (if we ever have such an event)

      # if we have already done this one (after regularization), we skip
      # e.g. this is i_jk and we had done c_ab; this one is c_ba, which gets regularized to c_ab, so we don't need to do it again.
      if(specific_event_name %in% names(out)){
        cat(generic_event_name, " -> ", specific_event_name, " already entered -- skipping.\n")
        next
      }

      # if we said merge_adjacent_pivot_events, we need to check whether we have already stored a pivot event adjacent to this one
      if(merge_adjacent_pivot_events & ! non_pivot_event){
        # get the vector of adjacent_events from the event_list spec
        generic_adjacent_event_names <- event_list[[generic_event_name]]$adjacent_events
        if(!is.null(generic_adjacent_event_names)){
          # regularize the names
          specific_adjacent_event_names <- generic_adjacent_event_names %>%
            map_chr(turn_generic_to_specific) %>%  # i_jk -> c_ba
            map_chr(regularize_candidate_order_after_underscore)  # c_ba -> cab
          cat(generic_event_name, " can be merged with ", generic_adjacent_event_names %>% paste(collapse = ", "), ". After transformation, looking to merge ", specific_event_name, " with ", specific_adjacent_event_names %>% paste(collapse = ""), ".\n", sep = "")
          cat("names(out) is ", names(out) %>% paste(collapse = ", "), ".\n")
          in_out <- which(specific_adjacent_event_names %in% names(out))[1]
          if(!is.na(in_out)){
            cat("Found one! ", specific_adjacent_event_names[in_out], ".\n", sep = "")
            # store existing result from this adjacent event
            out[[specific_event_name]] <- out[[specific_adjacent_event_names[in_out]]]
            # but not the P matrix -- that should be from this event
            out[[specific_event_name]]$W_mat <- this_P[cand_order, ballot_order]
            # store time
            time_diff <- Sys.time() - this_time_start
            units(time_diff) <- "secs"
            out[[specific_event_name]]$seconds_elapsed <- as.double(time_diff)
            next  # on to the next event
          }
        }
      }

      cat(generic_event_name, " -> ", specific_event_name, ".\n", sep = "")

      if(!from_sims){

        # get the S array, either from storage or from conditions
        if(generic_event_name %in% names(S_list)){
          this_S <- S_list[[generic_event_name]]
        }else{
          this_S <- S_array_from_inequalities_and_conditions(event_list[[generic_event_name]]$conditions, rows_to_alter = event_list[[generic_event_name]]$rows_to_alter, drop_dimension = drop_dimension, limits = limits)  # qhull options, epsilon
          S_list[[generic_event_name]] <- this_S
        }

        # compute the integral
        if(distribution == "dirichlet"){
          out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[ballot_order], maxEvals = maxEvals, tol = tol, ...)
          # when we allowed drop_dimension for compound pivot events, we would get a point (e.g. at a_bc for k = 3). adaptIntegrateSimplex didn't know what to do, so we had to specify "give me the density at the point. now this shouldn't happen.
          #         out[[specific_event_name]] <- list(integral = gtools::ddirichlet(as.vector(this_S), alpha[ballot_order])/sqrt(length(alpha)), functionEvaluations = 1, message = "OK")
        }else if(distribution == "logisticnormal"){
          out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[ballot_order], sigma = sigma[ballot_order, ballot_order], maxEvals = maxEvals, tol = tol, ...)
        }
      }else{
        cat(".")
        # get the conditions
        this_im <- event_list[[generic_event_name]]$conditions
        C_mat <- this_im[,-ncol(this_im)]
        rta <- event_list[[generic_event_name]]$rows_to_alter
        # permute the sims -- could permute the conditions, but I can't see an advantage. (thought we could skip some computations, but there should be no repeats)
        these_sims <- sims[,ballot_order]
        # apply the conditions
        XA_prime <- these_sims %*% t(C_mat)
        if(length(rta) == 0){
          proportion_in_window <- mean(apply(XA_prime, 1, all_positive))
        }else{
          proportion_in_window <- mean(apply(XA_prime[,-rta, drop = F], 1, all_positive) & apply(XA_prime[, rta, drop = F], 1, all_in_window, window = limits))
        }
        this_window_factor <- ifelse(length(rta) > 0, sqrt(2)*(limits[2] - limits[1]), 1)
        out[[specific_event_name]] <- list(integral = proportion_in_window/this_window_factor) # unscaled -- scaling factor applied below.
      }

      # when we drop a dimension we may need to scale the integral by the area we have lost
      if(drop_dimension & !is.null(event_list[[generic_event_name]]$scaling_factor)){
        out[[specific_event_name]]$integral <- out[[specific_event_name]]$integral*event_list[[generic_event_name]]$scaling_factor
      }

      # store the permuted P matrix
      out[[specific_event_name]]$W_mat <- this_P[cand_order, ballot_order]

      # on the clock?
      if(store_time){
        time_diff <- Sys.time() - this_time_start
        units(time_diff) <- "secs"
        out[[specific_event_name]]$seconds_elapsed <- as.double(time_diff)
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

regularize_candidate_order_after_underscore <- function(specific_event_name){
  if(str_detect(specific_event_name, "_") & !str_detect(specific_event_name, "\\|")){
    split_name <- specific_event_name %>% str_split("") %>% unlist()
    if(length(split_name) > 3){
      specific_event_name <- paste0(paste(split_name[1:2], collapse = ""), paste(sort(split_name[3:length(split_name)]), collapse = ""))
    }
  }
  specific_event_name
}

all_positive <- function(vec){
  all(vec > 0)
}

all_in_window <- function(vec, window){
  all(vec > window[1] & vec < window[2])
}
