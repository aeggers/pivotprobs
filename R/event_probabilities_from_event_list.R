#' Compute probability of election events
#'
#' \code{election_event_probs()} takes an \code{election} argument
#' (created by e.g. \code{plurality_election()} or \code{irv_election()}) and
#' returns a list containing, for each election event (i.e. each class of
#' election outcomes),
#' \itemize{
#' \item a scalar \code{integral} indicating the probability of the event,
#' \item a matrix \code{P}  indicating which candidate elected (rows) as a
#' function of which action the voter takes (columns), and
#' \item other objects depending on the method used to compute the event probabilities.}
#'
#' There are four methods for computing the probability of a given election event:
#' \itemize{
#' \item "sc" (or "SC" or "SimplicialCubature") partitions the election event
#' into disjoint simplices and uses the SimplicialCubature to integrate the
#' belief function over these simplices using the adaptive algorithm of
#' Genz and Cools. Slow, but can be sped up  by skipping non_pivot_events
#' (\code{skip_non_pivot_events=T}), dropping a dimension
#' (\code{drop_dimension=T}), merging adjacent pivot events
#' (\code{merge_adjacent_pivot_events=T}), and/or (when relevant) skipping compound
#' pivot events (\code{skip_compound_pivot_events=T}). Currently implemented only for
#' Dirichlet and Logistic Normal beliefs, but with some adjustment to the syntax
#' could be extended to handle other distributions (get in touch!).
#' \item "mc" (or "MC" or "Monte Carlo") counts the proportion of simulated elections
#' (supplied by the user or generated from user-supplied parameters) that satisfy
#' election event conditions. Also slow (especially with large \code{num_sums}),
#' but like "sc" can be sped up. Can generate simulations from Dirichlet or Logistic Normal belief distributions, or can be applied to any matrix of simulated elections
#' supplied by the user.
#' \item "ev" (or "EV" or "Eggers-Vivyan") uses numerical integration of the
#' Dirichlet distribution to estimate pivot event probabilities in plurality
#' elections, making an independence assumption when there are
#' more than four candidates. Fast but overstates the probability of low-probability
#' events and limited to Dirichlet beliefs.
#' \item "en" (or "EN" or "Eggers-Nowacki") uses numerical integration of the
#' Dirichlet distribution to estimate pivot event probabilities in three-candidate
#' IRV elections up to an arbitrary degree of precision. Fast but limited to Dirichlet beliefs.
#' }
#'
#' @param election A list with elements \code{n} (electorate size), \code{ordinal}
#'  (logical flag indicating whether it is an ordinal system or not), and \code{events}, a list of election events each with elements \code{win_conditions}, \code{tie_condition_rows}, and \code{P}. See ?event_list_functions.
#' @param method Method for estimating event probabilities: one of "sc"
#'  (SimplicialCubature), "mc" (Monte Carlo), "ev" (Eggers-Vivyan), or
#'  "en" (Eggers-Nowacki). See below for details.
#' @param mc_method Method for Monte Carlo estimation of pivot
#' event probabilities: one of
#' "rectangular", "density", "naive_density".
#' @param alpha Optional vector of parameters for Dirichlet distribution (one for
#' each ballot type). In a plurality election with k candidates, should be of
#' length k; in an ordinal system with 3 candidates, should be length 6.
#' @param mu Optional vector of location parameters for Dirichlet distribution
#' or logistic normal distribution.
#' @param precision Optional scalar precision parameter for Dirichlet distribution.
#' @param sigma Optional covariance matrix for logistic normal distribution.
#' @param sims Optional matrix of simulated elections, one column per ballot type.
#' @param num_sims  Optional number of simulated elections, if not
#' provided in \code{sims}.
#' @param sim_window Vote share discrepancy within which two candidates will be
#' "nearly tied" for method "mc". Larger \code{sim_window} means more bias,
#' smaller means more variance. Set to .01 by default.
#' @param cand_names Optional vector of candidate names. If not supplied, will
#' be "a", "b", "c", etc.
#' @param drop_dimension Method "sc" can be sped up for pivot events by
#' integrating the belief function on facets where one candidate is exactly
#' 1/2n behind the other (or they are exactly tied, if \code{drop_dimension=T})
#' rather than in the hypervolume where two candidates are nearly tied.
#' @param merge_adjacent_pivot_events If \code{F}, method "sc" computes the probability of e.g. candidate \code{a} being just behind \code{b} and vice
#' versa; if \code{T}, method "sc" saves time by computing the probability of candidates \code{a} and \code{b} being approximately or exactly tied (depending on \code{drop_dimension}).
#' @param skip_non_pivot_events Set to \code{T} to avoid computing the probability
#' of events where one vote cannot make a difference, e.g. where candidate \code{a}
#' wins by more than one vote. Relevant for methods "sc" and "mc".
#' @param skip_compound_pivot_events Set to \code{T} to avoid computing the probability of (near) three-way ties in plurality.
#' @param force_condition_based_mc Set to \code{T} to force Monte Carlo simulations to
#' use the conditions in the election object rather than a faster bespoke method.
#' @param store_time By default we store the time each computation takes, in seconds.
#' @param maxEvals The maximum number of function evaluations to compute an
#' integral via method "sc" (passed to
#' \code{SimplicialCubature::adaptIntegrateSimplex()}). If this is exceeded, you still get a result, but it will not be as precise as requested. You
#' can check whether the limit kicked in by looking at the
#' `message` attribute of each event probability. To prevent the limit
#' from binding, the default is a very large value; the computation
#' may therefore take a very long time.
#' @param tol The relative error requested in computing an
#' integral via method "sc" (passed to
#' \code{SimplicialCubature::adaptIntegrateSimplex()}). For faster imprecise computation (e.g. for testing), set to e.g. .2.
#' @param eval_points_for_1d_integral Sets the \code{subdivisions} argument to \code{integrate()} when
#' \code{SimplicialCubature::adaptIntegrateSimplex()} will be performed on a line.
#' @param ev_increments Increments for numerical integration in method "ev".
#' @param en_increments_1st_round Increments for numerical integration of
#' first-round pivot events in method "en".
#' @param en_increments_2nd_round Increments for numerical integration of
#' second-round pivot events in method "en".
#' @param bw_divisor If ks::hpi determines optimal bandwidth
#' for density estimation is h, we use h/bw_divisor when
#' \code{mc_method="density"} or
#' \code{mc_method="naive_density"}. Allows for undersmoothing.
#'
#' @return A list containing one list per election event. Each of these lower-level
#' lists contains
 #'\itemize{
#' \item \code{integral} the event probability
#' \item \code{P} the election probability matrix, indicating which candidate (rows)
#' is elected at this event depending on which action (columns) the voter takes.
#' \item \code{seconds_elapsed} the time to compute this event probability.
#' }
#' For method "sc" these event-specific lists include other output from \code{SimplicialCubature::adaptIntegrateSimplex()}, such as \code{functionEvaluations} and \code{message}.
#'
#' @examples
#' alpha3 <- c(.4, .35, .25)*85
#' sc_out <- plurality_election(k = 3) %>%
#'   election_event_probs(method = "sc", alpha = alpha3, tol = .1)
#' sc_out[["a_b"]]$integral
#' sc_out[["a_b"]]$seconds_elapsed
#'
#' mc_out <- plurality_election(k = 3) %>%
#'   election_event_probs(method = "mc", alpha = alpha3, num_sims = 500000)
#' mc_out[["a_b"]]$integral
#' mc_out[["a_b"]]$seconds_elapsed
#'
#' mc_out <- plurality_election(k = 3) %>%
#'   election_event_probs(method = "ev", alpha = alpha3)
#' mc_out[["a_b"]]$integral
#' mc_out[["a_b"]]$seconds_elapsed
#' @export
election_event_probs <- function(election,
                                 method = "sc",  # sc, mc, ev, en
                                 mc_method = "density", # "density", "naive_density", "rectangular"
                                 # distribution parameters
                                 alpha = NULL, # dirichlet
                                 mu = NULL,  # dirichlet or logisticnormal
                                 precision = NULL, # dirichlet
                                 sigma = NULL, # logisticnormal
                                 sims = NULL, # provide MC simulations
                                 num_sims = 100000, # for drawing
                                 sim_window = .01,
                                 cand_names = NULL, # optional
                                 drop_dimension = F,
                                 merge_adjacent_pivot_events = F,
                                 skip_non_pivot_events = T,
                                 skip_compound_pivot_events = T,
                                 minimum_volume = 0,
                                 force_condition_based_mc = F,
                                 store_time = T,
                                 maxEvals = .Machine$integer.max, tol = .01, eval_points_for_1d_integral = 100, # SimplicialCubature arguments
                                 ev_increments = 50,
                                 en_increments_1st_round = 30,
                                 en_increments_2nd_round = 100,
                                 bw_divisor = 2,
                                 ... # other arguments to SimplicialCubature::adaptIntegrateSimplex()
                                 ){

  # start the clock
  time_start <- Sys.time()

  sc_method_names <- c("sc", "SC", "SimplicialCubature")
  mc_method_names <- c("mc", "MC", "Monte Carlo")
  ev_method_names <- c("ev", "EV", "Eggers-Vivyan")
  en_method_names <- c("en", "EN", "Eggers-Nowacki")

  if(! method %in% c(sc_method_names, mc_method_names, ev_method_names, en_method_names)){
    stop("I don't recognize the supplied method.")
  }

  ## get the meta parametesr from the election object
  if(class(election) != "list"){stop("Expection `election` (first argument) to be a list.")}
  # ballot type
  ordinal <- election$ordinal
  # electorate size
  if(is.null(ordinal)){stop("election$ordinal must be specified so that we know whether this is an ordinal voting method or not.")}
  n <- election$n
  if(is.null(n)){stop("election$n must be specified so that we know the electorate size (which affects limits of integration and normalization factors.")}

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
    mu <- apply(sims, 2, mean)
  }

  # if the user passes sims, we override the method provided.
  if(!is.null(sims)){
    if(!method %in% mc_method_names){
      warning("You provided `sims` but specified a method or than simulation. Setting method to 'mc' (Monte Carlo).")
      method <- "mc"
    }
  }else if(method %in% ev_method_names){
    if(ordinal){stop("The Eggers-Vivyan method can only be used for plurality elections. Your election argument says it is an ordinal method.")}
    if(distribution != "dirichlet"){stop("The Eggers-Vivyan method can only be used with the Dirichlet distribution. The parameters you supplied do not indicate a Dirichlet distribution.")}
    skip_compound_pivot_events <- T
    skip_non_pivot_events <- T
    merge_adjacent_pivot_events <- T # TODO: consider relaxing this.
  }else if(method %in% en_method_names){
    if(!ordinal){stop("The Eggers-Nowacki method can only be used for IRV elections. Your election argument says it is a non-ordinal method.")}
    skip_compound_pivot_events <- T
    skip_non_pivot_events <- T
    merge_adjacent_pivot_events <- T # TODO: consider relaxing this.
  }

  # for simulations with mc_method = "rectangular", we need to merge_adjacent_pivot_events. we also set drop_dimension = F so that we don't apply the scaling factor.
  # you can think of it as dropping a dimension, but then the scaling factors cancel out.
  # better to think of it like this: the ratio of the probability of being within sim_window to the pivot probability is the ratio of sim_window to 1/n. so pivot probability is probability of being within sim_window over (sim_window times n)
  if(method %in% mc_method_names){
    if(is.null(sims)){
      # we need to draw simulations.
      if(!is.numeric(num_sims)){stop("If draw_sims is TRUE, we need a numeric num_sims.")}
      if(!num_sims > 0){stop("num_sims needs to be a positive number.")}
      if(num_sims < 10000){warning(paste0("num_sims is only ", num_sims, "; probably want more."))}
      # now, generate the simulations
      if(distribution == "dirichlet"){
        sims <- gtools::rdirichlet(n = num_sims, alpha = alpha)
      }else if(distribution == "logisticnofmal"){
        sims <- rlogisticnormal(n = num_sims, mu = mu, sigma = sigma)
      }
    }
    if(!force_condition_based_mc & election$system %in% c("plurality", "positional", "irv", "kemeny_young")){
      # we use a "standalone" method for simulations, because these are faster
      if(election$system == "plurality"){
        return(plurality_event_probs_from_sims(sims = sims, n = n, window = sim_window, cand_names = cand_names, method = mc_method, drop = drop_dimension, merge = merge_adjacent_pivot_events, bw_divisor = bw_divisor, skip_non_pivot_events = skip_non_pivot_events))
      }else if(election$system == "positional"){
        return(positional_event_probs_from_sims(sims = sims, n = n, window = sim_window, cand_names = cand_names, method = mc_method, merge = merge_adjacent_pivot_events, s = election$s, bw_divisor = bw_divisor))
      }else if(election$system == "irv"){
        return(irv_event_probs_from_sims(sims = sims, n = n, window = sim_window, cand_names = cand_names, method = mc_method, merge = merge_adjacent_pivot_events, s = election$s, bw_divisor = bw_divisor))
      }else if(election$system == "kemeny_young"){
        return(condorcet_event_probs_from_sims(sims = sims, n = n, window = sim_window, cand_names = cand_names, kemeny = T, method = mc_method, merge = merge_adjacent_pivot_events, bw_divisor = bw_divisor))
      }else{
        # if we are using the election conditions -- legacy method, basically
        merge_adjacent_pivot_events <- T
        drop_dimension <- F
      }
    }
  }

  # for the limits we define below, we want the value of one vote in terms of vote share, i.e. 1/n.
  base_limit <- 1/n

  # set limits (width of integration region) based on method and arguments
  if(method %in% mc_method_names){
    # legacy: if we are using the election conditions to do Monte Carlo, we are using the rectangular method, merging adjacent, and just counting cases where the key margin is between -1/(2*sim_window) and 1/(2*sim_window)
    limits <- (sim_window/2)*c(-1, 1)
  }else if(merge_adjacent_pivot_events & drop_dimension){
    limits <- c(0, 0) # only first limit used -- a line at the boundary
  }else if(merge_adjacent_pivot_events){
    limits <- (base_limit/2)*c(-1, 1)  # a thin strip straddling the boundary: -1/2n, 1/2n
  }else if(drop_dimension){
    limits <- rep(base_limit/2, 2)   # only first limit used -- a line just short of the boundary
  }else{
    limits <- c(0, base_limit)     # a thin strip short of the boundary -- 0 to 1/n
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
  }else{  # single-ballot system
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

  ## Special case: if we're going to be doing a line integral, we need maxEvals to not be so high
  # Would be better if condition were based on S rather than election method, but S could in principle be an array including lines and other things.
  if(election$system == "plurality" && election$k == 3 & drop_dimension == T){
    # cat("Setting maxEvals to lower number to avoid passing subdivisions=NA to integrate().")
    maxEvals <- eval_points_for_1d_integral*21L
  }

  # now we permute through possible orderings of candidates
  for(p_row in 1:nrow(candidate_permutation_mat)){
    cand_vector <- candidate_permutation_mat[p_row,]  # e.g. c("c", "a", "b")

    cand_param_indexes <- rank(cand_vector) # this allows us to reorder the parameters correctly for this permutation, assuming parameters ranked by candidate name. e.g. "cab" -> (3, 1, 2)
    cand_order <- order(cand_vector) # this allows us to reorganize the P matrix, e.g. cab -> c(2, 3, 1)
    if(ordinal){
      obv <- ordered_ballot_vector_from_candidate_vector(cand_vector)  # e.g. cab, cba, acb, abc, bca, bac
      # these are the indices for reordering the parameters such that e.g. c becomes i, a becomes j, b becomes k.
      ballot_param_indexes <- obv %>% names() %>% as.numeric()   # (5, 6, 2, 1, 4, 3)
      # these are the indices for reordering the ballots in the P matrix
      ballot_order <- c(order(obv), length(obv) + 1)
    }else{
      ballot_param_indexes <- cand_param_indexes  # (3, 1, 2)
      ballot_order <- c(cand_order, length(cand_order) + 1)                  # (2, 3, 1, 4)
    }

    if(class(election$events) != "list"){
      stop("election$events must be a list of election events.")
    }

    event_list <- election$events

    for(el_index in 1:length(names(election$events))){

      generic_event_name <- names(event_list)[el_index]

      # start the clock on this election event
      this_time_start <- Sys.time()

      this_event <- event_list[[generic_event_name]]
      this_P <- this_event$P

      # we use the P matrix to determine if this is a pivot event
      pivot_event <- this_P %>% apply(1, mean) %>% max() < 1
      if(skip_non_pivot_events & !pivot_event){next}

      # we can specify skip_compound_pivot_events, i.e. those with more than one "row to alter"; we MUST skip compound pivot events when drop_dimension is TRUE.
      if(length(this_event$tie_condition_rows) > 1 & (skip_compound_pivot_events | drop_dimension )){next}

      # if no scaling factor is provided we set it to 1
      scaling_factor <- this_event$scaling_factor
      if(is.null(scaling_factor)){scaling_factor <- 1}

      # we need to convert the generic_event_name (e.g. i_j) to a specific name (e.g. c_a) based on the permutation
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
        regularize_candidate_order_after_underscore()  # c_ba => c_ab, but c_ba|cb stays same

      # Now we check to see if we need to store this one.
      store <- TRUE
      if(specific_event_name %in% names(out)){
        # There are two reasons this occurs.
        # 1) In Eggers-Nowacki we get all the first-round pivot events for a given pair at once (because more efficient that way), so when we get to i_j|ik we will have already stored it.
        # 2) In plurality we have compound pivot events that are identical but not until after regularization. e.g. this is i_jk, which translates to c_ba before regularization but c_ab after regularization and we had already done c_ab.
        store <- FALSE
      } else if(merge_adjacent_pivot_events & pivot_event){
        # if we said merge_adjacent_pivot_events, we need to check whether we have already stored a pivot event adjacent to this one
        # get the vector of adjacent_events from the event_list spec
        generic_adjacent_event_names <- this_event$adjacent_events
        if(!is.null(generic_adjacent_event_names)){
          # regularize the names
          specific_adjacent_event_names <- generic_adjacent_event_names %>%
            purrr::map_chr(turn_generic_to_specific) %>%  # i_jk -> c_ba
            purrr::map_chr(regularize_candidate_order_after_underscore)  # c_ba -> cab
          # Now see if any of the adjacent pivot events have already been measured
          in_out <- which(specific_adjacent_event_names %in% names(out))[1]
          if(!is.na(in_out)){
            # if so, copy from this adjacent event
            store <- FALSE
            out[[specific_event_name]] <- out[[specific_adjacent_event_names[in_out]]]
            # but wipe the P matrix and time -- these should be from this iteration.
            out[[specific_event_name]]$P <- NULL
            out[[specific_event_name]]$seconds_elapsed <- NULL
            # we fill them in below
          }
        }
      }

      # Now we store the integral if store == T
      if(store){
        if(method %in% sc_method_names){
          # simplicial cubature to integrate a function

          # get the S array, either from storage or from conditions
          if(generic_event_name %in% names(S_list)){
            this_S <- S_list[[generic_event_name]]
          }else{
            this_S <- S_array_from_inequalities_and_conditions(this_event$win_conditions, rows_to_alter = this_event$tie_condition_rows, drop_dimension = drop_dimension, limits = limits)  # qhull options, epsilon
            if(!is.null(this_S) & !is.null(minimum_volume) & minimum_volume > 0){
              # some matrices in the S array produce an integral of zero
              # we can filter these out here. I thought this would speed thingds up but it didn't. I am leaving in here in case we can figure out a better way later.
              # I considered also having the filtering built into the election_list
              # methods because e.g. it's the same for every Condorcet election,
              # but we have so much flexibility with s in positional and IRV
              # methods that it didn't look very elegant/general.
              # this is a general solution that doesn't seem to solve anything.
              indices_to_keep <- c()
              for(k in 1:(dim(this_S)[3])){
                this_val <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S[,,k], alpha = rep(1, 6), maxEvals = 100000, tol = .2)
                if(this_val$integral > minimum_volume){
                  cat(".")
                  indices_to_keep <- c(indices_to_keep, k)
                }else{
                  cat("-")
                }
              }
              cat("\n")
              # subset the S array
              this_S <- this_S[,,indices_to_keep]
            }
            S_list[[generic_event_name]] <- this_S
          }

          if(is.null(this_S)){
            # the conditions aren't met anywhere,
            # e.g. if first round of IRV is Borda count, event i_j|ij cannot occur.
            out[[specific_event_name]] <- list(integral = 0)
          }else{
            # compute the integral
            if(distribution == "dirichlet"){
              # for error diagnosis, output the arguments we pass to SimplicialCubature
              # if(election$system == "irv"){cat("Saving arguments for ", generic_event_name, ".\n"); save(generic_event_name, alpha, ballot_param_indexes, maxEvals, this_S, tol, file = paste0("~/Dropbox/research/strategic_voting/pivotprobs_paper/data/this_S_etc_n_", generic_event_name, "_", election$n, ".RData"))}
              out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = dirichlet_for_integration, S = this_S, alpha = alpha[ballot_param_indexes], maxEvals = maxEvals, tol = tol, ...)
              # when we allowed drop_dimension for compound pivot events, we would get a point (e.g. at a_bc for k = 3). adaptIntegrateSimplex didn't know what to do, so we had to specify "give me the density at the point":
              #         out[[specific_event_name]] <- list(integral = gtools::ddirichlet(as.vector(this_S), alpha[ballot_order])/sqrt(length(alpha)), functionEvaluations = 1, message = "OK")
              # now we don't need this because we shouldn't be trying to do compound pivot events with drop_dimension
            }else if(distribution == "logisticnormal"){
              out[[specific_event_name]] <- SimplicialCubature::adaptIntegrateSimplex(f = logisticnormal_for_integration, S = this_S, mu = mu[ballot_param_indexes], sigma = sigma[ballot_param_indexes, ballot_param_indexes], maxEvals = maxEvals, tol = tol, ...)
            }
          }

          # when we drop a dimension we need to scale the integral by the width of the interval where a single vote can make a difference
          if(drop_dimension & pivot_event){
            out[[specific_event_name]]$integral <- out[[specific_event_name]]$integral*(1/(n*scaling_factor))
            # note if we deal with compound events we probably have to square this.
          }
        }else if(method %in% ev_method_names){
          # Eggers Vivyan is easy: tie for first between first two candidates (based on alpha vector). So no information from the pivot event list necessary.
          out[[specific_event_name]] <- list(integral = eggers_vivyan_probability_of_tie_for_first(alpha = alpha[ballot_param_indexes], increments = ev_increments)$estimate*(1/(n*scaling_factor)))
        }else if(method %in% en_method_names){
          # Eggers Nowacki for IRV elections
          # detect second round or first round from name
          first_round <- specific_event_name %>% str_detect("\\|")
          if(first_round & this_event$scaling_factor == sqrt(6)){stop("Name suggests this is a first-round pivot event in IRV, but the scaling factor is consistent with a second-round pivot event.")}
          # if second round
          if(!first_round){
            # do normally
            the_estimate <- eggers_nowacki_second_round_pivot_probability(alpha = alpha[ballot_param_indexes], increments = en_increments_2nd_round)
            out[[specific_event_name]] <- list(integral = the_estimate*(1/n)) # scaling factor is in the code
          }else{
            # if first round, we get all three first round pivot events at once
            the_estimates <- eggers_nowacki_first_round_pivot_probabilities(alpha = alpha[ballot_param_indexes], increments = en_increments_1st_round)
            # this returns i_j|ij, i_j|kj, i_j|ik
            for(this_generic_event_name in names(the_estimates)){
              this_specific_event_name <- this_generic_event_name %>%
                turn_generic_to_specific()
              out[[this_specific_event_name]] <- list(integral = the_estimates[[this_generic_event_name]]*(1/n)) # scaling factor is in the code
            }
          }
        }else if(method %in% mc_method_names){
          # Monte Carlo simulation
          # this code is elegant but slow so not used.
          # get the event conditions as stated in the event_list
          this_im <- this_event$win_conditions
          C_mat <- this_im[,-ncol(this_im)] # take off the vector of constants (1/n) -- we are checking Ax >= 0. Could change that -- would make a tiny difference.
          rta <- this_event$tie_condition_rows
          # permute the sims -- could permute the conditions, but I can't see an advantage. (thought we could skip some computations, but there should be no repeats)
          these_sims <- sims[,ballot_param_indexes]
          # apply the conditions
          XA_prime <- these_sims %*% t(C_mat)
          if(length(rta) == 0){ # non pivot event
            proportion_in_window <- mean(apply(XA_prime, 1, all_positive)) # do all the conditions hold?
          }else{
            proportion_in_window <- mean(
              apply(XA_prime[,-rta, drop = F], 1, all_positive) & # all raw conditions hold
                apply(XA_prime[, rta, drop = F], 1, all_in_window, window = limits)) # and the rows_to_alter conditions area near equality
          }
          mc_scaling_factor <- ifelse(pivot_event, (sim_window*n)^length(rta), 1)
          out[[specific_event_name]] <- list(integral = proportion_in_window/mc_scaling_factor) # this is our estimate of the integral.
        }
      }

      # always store the permuted P matrix
      out[[specific_event_name]]$P <- this_P[cand_order, ballot_order]

      # on the clock?
      if(store_time){
        time_diff <- Sys.time() - this_time_start
        units(time_diff) <- "secs"
        out[[specific_event_name]]$seconds_elapsed <- as.double(time_diff)
      }
    }
  }

  # previously stored total time, but makes the output ugly, i.e. a list of events and then total time.
  # if(store_time){
  #   time_diff <- Sys.time() - time_start
  #   units(time_diff) <- "secs"
  #   out[["total"]]$seconds_elapsed <- as.double(time_diff)
  # }

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
