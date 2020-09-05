# I think this is the way forward.

event_probabilities_from_event_list <- function(event_list = NULL, alpha = NULL, mu = NULL, precision = NULL, sigma = NULL, sims = NULL, draw_sims = F, num_sims = 100000, cand_names = NULL, ordinal = T, drop_dimension = F, merge_adjacent_pivot_events = F, skip_non_pivot_events = F, store_time = T, ...){

  if(is.null(event_list)){stop("You have to pass an event_list.")}

  # storage for output
  out <- list()

  # start the clock
  time_start <- Sys.time()

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
    distribution <-  "sims"
  }

  # get 1/n from the conditions -- is this risky?
  conditions_1 <- event_list[[1]]$conditions
  base_limit <- conditions_1[1, ncol(conditions_1)]
  if(class(base_limit) != "numeric"){
    stop(paste0("Was expecting last column of first row of first condition to be 1/n; instead it is ", base_limit, "."))
  }
  if(base_limit <= 0){
    stop("I assumed that the last column of the first condition on the first event is 1/n, but it is", base_limit, ".")
  }

  # set limits (width of integration region) based on arguments
  if(merge_adjacent_pivot_events & drop_dimension){
    limits <- c(0, 0) # only first limit used -- a line at the boundary
  } else if(merge_adjacent_pivot_events){
    limits <- c(-base_limit/2, base_limit/2)  # a thin strip straddling the boundary
  } else if(drop_dimension){
    limits <- rep(base_limit/2, 2)   # only first limit used -- a line just short of the boundary
  } else{
    limits <- c(0, base_limit)     # a thin strip short of the boundary
  }

  if(ordinal){
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
    if(length(mu) != 6){
      stop("For ordinal methods I expect parameters for 6 ballot types. You provided parameters indicating more or fewer than 6 ballot types.")
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
  colnames(candidate_permutation_mat) <- letters[9:(9 + length(cand_names) - 1)]

  # storage so we don't make the same S_array more than once
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

      this_P <- event_list[[names(event_list)[el_index]]]$P

      # we use the P matrix to determine if this is a pivot event
      non_pivot_event <- this_P %>% apply(1, mean) %>% max() == 1
      if(skip_non_pivot_events & non_pivot_event){next}

      # start the clock
      this_time_start <- Sys.time()

      # we need to convert e.g. i_j to c_a
      generic_event_name <- names(event_list)[el_index]

      # two functions I use a couple times here. uses local variables in definition.
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

      # we may have already done this one (after regularization)
      if(specific_event_name %in% names(out)){
        cat(generic_event_name, " -> ", specific_event_name, " already entered -- skipping.\n")
        next
      }

      # check whether we have already stored a pivot event adjacent to this one
      if(merge_adjacent_pivot_events & ! non_pivot_event){
        generic_adjacent_event_names <- event_list[[generic_event_name]]$adjacent_events
        if(!is.null(generic_adjacent_event_names)){
          specific_adjacent_event_names <- generic_adjacent_event_names %>%
            map_chr(turn_generic_to_specific) %>%
            map_chr(regularize_candidate_order_after_underscore)
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

      # get the S array, either from storage or from conditions
      if(generic_event_name %in% names(S_list)){
        this_S <- S_list[[generic_event_name]]
      }else{
        this_S <- S_array_from_inequalities_and_conditions(event_list[[generic_event_name]]$conditions, rows_to_alter = event_list[[generic_event_name]]$rows_to_alter, drop_dimension = drop_dimension, limits = limits)  # qhull options, epsilon
        S_list[[generic_event_name]] <- this_S
      }

      # compute the integral
      if(distribution == "dirichlet"){
        if(class(this_S) == "matrix" && ncol(this_S) == 1){
          out[[specific_event_name]] <- list(integral = gtools::ddirichlet(as.vector(this_S), alpha[ballot_order])/sqrt(length(alpha)), functionEvaluations = 1, message = "OK")
        }else{
          out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[ballot_order], ...)
        }
      }else if(distribution == "logisticnormal"){
        if(class(this_S) == "matrix" && ncol(this_S) == 1){
          out[[specific_event_name]] <- list(integral = dlogisticnormal(as.vector(this_S), mu[ballot_order], sigma[ballot_order, ballot_order])/sqrt(length(mu)), functionEvaluations = 1, message = "OK")
        }else{
          out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[ballot_order], sigma = sigma[ballot_order, ballot_order], ...)
        }
      }

      # when we drop a dimension we need to scale the integral by the area we have lost
      if(drop_dimension){
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

