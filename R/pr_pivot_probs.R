
dhondt_pivot_probs_cands_1_and_2_dirichlet <- function(alpha_vec, M, increments = 50){

  piv_prob_list = list()
  for(x in 1:M){
    for(y in 1:M){
      if(x + y > M + 1){next}
      this.piv.prob = 0
      min.a = x/(M + 2)
      max.a = x/(M + 1)
      start.point = c(min.a, min.a*y/x, 1 - min.a - min.a*y/x)
      end.point = c(max.a, max.a*y/x, 1 - max.a - max.a*y/x)
      distance.covered = sqrt(sum((start.point - end.point)^2))
      as = seq(min.a, max.a, length = increments)
      for(i in 1:(length(as) - 1)){
        a = (as[i] + as[i+1])/2
        this.piv.prob = this.piv.prob + gtools::ddirichlet(x = c(a, a*y/x, 1 - a - a*y/x), alpha = alpha_vec)/sqrt(3)
      }
      piv_prob_list[[paste0(x, ".", y)]] = this.piv.prob*(distance.covered/increments)
    }
  }
  piv_prob_list
}

dhondt_pivot_probs_dirichlet <- function(alpha_vec, M, increments = 50, cand_names = c("a", "b", "c"), n = 10000){
  stopifnot(length(alpha_vec) == 3)
  stopifnot(length(cand_names) == 3)
  ppdd.1.2 <- dhondt_pivot_probs_cands_1_and_2_dirichlet(alpha_vec, M, increments)
  ppdd.1.3 <- dhondt_pivot_probs_cands_1_and_2_dirichlet(alpha_vec[c(1,3,2)], M, increments)
  ppdd.2.3 <- dhondt_pivot_probs_cands_1_and_2_dirichlet(alpha_vec[c(2,3,1)], M, increments)
  out <- list()
  out[[paste0(cand_names[1], cand_names[2])]] <- list("integral" = sum(unlist(ppdd.1.2))/n,
                                                      "raw" = ppdd.1.2)
  out[[paste0(cand_names[1], cand_names[3])]] <- list("integral" = sum(unlist(ppdd.1.3))/n,
                                                      "raw" = ppdd.1.3)
  out[[paste0(cand_names[2], cand_names[3])]] <- list("integral" = sum(unlist(ppdd.2.3))/n,
                                                      "raw" = ppdd.2.3)
  out
}

divisor_gridpoints_ab <- function(M, divisors = "dhondt"){

  if(length(divisors) > 1){
    divs <- divisors
  }else if(divisors == "dhondt"){
    divs <- 1:(M+1)
  }else if(divisors == "sainte-lague"){
    divs <- 1 + 2*(0:M)
  }else if(divisors == "koudeika"){
    divs <- c(sqrt(2), 2:(M+1))
  }else{
    stop("We don't recognize that kind of divisor.")
  }

  # for a-b ties where c is getting 0 seats, we need a divisor that corresponds to getting 0 seats. we set this at 1e-10.
  divs <- c(1e-10, divs)

  # "pot_seats" refers to the number of seats the candidate could win given one more vote; a_pot_seats = 1 and b_pot_seats = 1 means one of them is going to win one seat, but not both
  expand_grid(a_pot_seats = 1:M, b_pot_seats = 1:M, c_boost = c(0,1)) %>%
    filter(a_pot_seats + b_pot_seats <= M + 1) %>%
    mutate(
      a_divisor = divs[1 + a_pot_seats],
      b_divisor = divs[1 + b_pot_seats],
      c_seats = M - a_pot_seats - b_pot_seats + 1,
      c_divisor = divs[1 + c_seats + c_boost],
      a = 1/(1 + (b_divisor + c_divisor)/a_divisor),
      b = a*b_divisor/a_divisor,
      c = 1 - a - b
    ) %>%
    select(a, b, c, a_pot_seats, b_pot_seats, c_seats)
}

divisor_gridpoints <- function(M, divisors = "dhondt"){

  dgp <- divisor_gridpoints_ab(M, divisors = divisors) %>%
    mutate(a_low = a_pot_seats - 1,
           b_low = b_pot_seats - 1,
           a_high = a_pot_seats,
           b_high = b_pot_seats) %>%
    select(a_low, a_high, b_low, b_high, c_lowhigh = c_seats, a, b, c)

  bind_rows(
    # ab
    dgp %>%
      mutate(
        event_1 = pmap_chr(select(., a_low, b_high, c_lowhigh), str_c, sep = "."),
        event_2 = pmap_chr(select(., a_high, b_low, c_lowhigh), str_c, sep = "."),
        event_class = "ab"
      ),
    # ac
    dgp %>%
      mutate(
        event_1 = pmap_chr(select(., a_low, c_lowhigh, b_high), str_c, sep = "."),
        event_2 = pmap_chr(select(., a_high, c_lowhigh, b_low), str_c, sep = "."),
        event_class = "ac"
      ) %>%
      rename(c = b,
             b = c),
    # bc
    dgp %>%
      mutate(
        event_1 = pmap_chr(select(., c_lowhigh, a_low, b_high), str_c, sep = "."),
        event_2 = pmap_chr(select(., c_lowhigh, a_high, b_low), str_c, sep = "."),
        event_class = "bc"
      ) %>%
      rename(
        a = c,
        b = a,
        c = b
      )
  ) %>%
    mutate(event = paste0(event_1, "_", event_2)) %>%
    select(event, event_class, a, b, c)

}



# # illustration
# divisor_gridpoints(M = 4, divisors = "dhondt") %>%
#   mutate(a = a + .5*b,
#          b = sqrt(3/4)*b) %>%
#   ggplot(aes(x = a, y = b)) +
#   geom_line(aes(group = event, col = event_class)) +
#   expand_limits(x = c(0,1), y = c(0,1)) +
#   coord_fixed() +
#   theme_void()
#
# divisor_gridpoints(M = 4, divisors = "sainte-lague") %>%
#   mutate(a = a + .5*b,
#          b = sqrt(3/4)*b) %>%
#   ggplot(aes(x = a, y = b)) +
#   geom_line(aes(group = event, col = event_class)) +
#   expand_limits(x = c(0,1), y = c(0,1)) +
#   coord_fixed() +
#   theme_void()
#
# divisor_gridpoints(M = 4, divisors = "koudeika") %>%
#   mutate(a = a + .5*b,
#          b = sqrt(3/4)*b) %>%
#   ggplot(aes(x = a, y = b)) +
#   geom_line(aes(group = event, col = event_class)) +
#   expand_limits(x = c(0,1), y = c(0,1)) +
#   coord_fixed() +
#   theme_void()
#
# # c(1, 4, 7, 11))
# divisor_gridpoints(M = 13, divisors = seq(1, 200, by = 1)) %>%
#   mutate(a = a + .5*b,
#          b = sqrt(3/4)*b) %>%
#   ggplot(aes(x = a, y = b)) +
#   geom_line(aes(group = event, col = event_class)) +
#   expand_limits(x = c(0,1), y = c(0,1)) +
#   coord_fixed() +
#   theme_void()



# Next to do:
# add rows for the non-pivot events for both plotting and integration.
# add that to event_list_functions so that we can get the P matrix including non pivot events.

outcome_names_from_M <- function(M){
  expand_grid(a = 0:M, b = 0:M) %>%
    filter(a + b <= M) %>%
    mutate(c = M - a - b) %>%
    mutate(outcome = pmap_chr(select(., a, b, c), str_c, sep = ".")) %>%
    pull(outcome)
}

P_matrix_from_PR_event_name_and_M <- function(event_name, M, cand_names = c("a", "b", "c"), reverse = F){
  # an event_name is like "0.1.11_1.0.11" -- it combines two outcome names
  outcome_names <- outcome_names_from_M(M)
  mat <- matrix(0, nrow = length(outcome_names), ncol = 4)
  rownames(mat) <- outcome_names
  colnames(mat) <- c(cand_names, "0")
  outcomes <- str_split(event_name, "_") %>% unlist()
  if(reverse){outcomes <- outcomes[c(2,1)]}
  outcome_seat_counts <- outcomes %>% str_split("\\.") %>% lapply(as.integer)
  outcome_seat_diff <- outcome_seat_counts[[1]] - outcome_seat_counts[[2]]
  mat[outcomes[1], -which(outcome_seat_diff < 0)] <- 1
  mat[outcomes[2], which(outcome_seat_diff < 0)] <- 1
  mat
}

