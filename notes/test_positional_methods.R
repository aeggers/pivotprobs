## test positional methods

# testing grid-based vs sampling.
# I think show for different grid sizes

N <- 100
n <- 100000
v_mat <- gtools::rdirichlet(N, alpha = c(5, 4, 3, 3, 2, 4))
s <- 20

alpha_vec <- as.numeric(v_mat[1,]*s)

bind_rows(
  positional_pivot_probs_grid_based(alpha_vec, increment = .05) %>% unlist(),
  positional_pivot_probs_grid_based(alpha_vec, increment = .035) %>% unlist(),
  positional_pivot_probs_grid_based(alpha_vec, increment = .02)
) -> gridders

gridders <- gridders %>% mutate(increment = c(.05, .035, .02), type = "grid_based")


sims <- gtools::rdirichlet(1000000, alpha_vec)
sim_based <- bind_rows(
  positional_piv_probs_simulation(sims, s = .5, tol = .03),
  positional_piv_probs_simulation(sims, s = .5, tol = .01)
) %>% mutate(increment = c(.03, .01), type = "sim")

bind_rows(gridders, sim_based) %>%
  pivot_longer(cols = ab:bc, names_to = "pp") -> for_plot

for_plot %>%
  ggplot(aes(x = increment, y = value, col = pp, linetype = type)) +
  geom_line(aes(group = interaction(pp, type))) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

# we have something that basically works.
# TODO: test more thoroughly. try to speed up?

# what testing approach?
# the problem is that we don't actually know the truth for any case.
# for the analytic approaches before, we showed that the answer gets closer and closer to the simulation answer as we increase the number of simulations.
# so maybe I want to do 3 sims and 3 grid approaches. then I could show the RMSE between the grid approach for a fixed increment and the sim approach for numbers of simulations



## update: in this version we do different grid sizes and with and without boosting low alphas. and we compare to one large simulation draw.
rmses <- tibble(j = NULL, increment = NULL, rmse = NULL)
Ns <- c(1000000)
J <- 100
set.seed(12345)
v_mat <- gtools::rdirichlet(J, alpha = c(3,2,2,2,1,3)/1.5)
ss <- runif(J, min = 10, max = 50)
for(j in 1:J){
  cat(j, ": ", sep = "")
  v_vec <- as.numeric(v_mat[j,])
  s <- ss[j]
  alpha_vec <- v_vec*s

  # simulations
  cat(" s ")
  df <- tibble(n = Ns)
  df %>%
    mutate(sims = pmap(., gtools::rdirichlet, alpha = alpha_vec),
           sim_result = map(sims, positional_piv_probs_simulation, s = .5, tol = .015)) %>%
    select(-sims) %>%
    unnest_longer(col = sim_result, values_to = "value", indices_to = "pp") -> sim_results

  # grid-based results
  cat(" g ")
  df <- expand_grid(increment = c(.025, .05), boost_low_alphas = c(T, F)) %>%
    mutate(gb_result = pmap(., positional_pivot_probs_grid_based, alpha_vec = alpha_vec)) %>%
    unnest_longer(col = gb_result, values_to = "value", indices_to = "pp") -> gb_results

  cat(" r")
  gb_results %>%
    left_join(sim_results %>% select(-n), by = "pp", suffix = c("_grid", "_sim")) %>%
    group_by(increment, boost_low_alphas) %>%
    summarize(rmse = sqrt(sum((value_grid - value_sim)^2))) -> these_rmses

  cat(".\n")

  rmses <- bind_rows(rmses, these_rmses %>% mutate(j = j))
}

# now to explore this? we want to compare rmses across four approaches.

rmses %>%
  ggplot(aes(x = factor(increment), y = rmse, group = interaction(increment, boost_low_alphas), fill = boost_low_alphas)) +
  geom_violin()

rmses %>% group_by(increment, boost_low_alphas) %>% summarize(mean(rmse))

# so looking at this I see that there is a lot of discrepancy between the simulation and the grid answer.
# things I can do:
  ## check the geometry to make sure that it's not excluding points.
ptg <- positional_tie_grid()

# how do I filter the data so that only two dimensions vary?
# if I fix three dimensions, any two of the remaining three form a line.
ptg %>% filter(round(ab, 4) == .2375 & round(ba, 4) == .1875 & round(ca, 4) == .0125) -> ptg2

ptg2 %>% ggplot(aes(x = ac, y = bc)) + geom_point(alpha = .1) + coord_fixed()
ptg2 %>% ggplot(aes(x = ac, y = cb)) + geom_point(alpha = .1) + coord_fixed()

ptg %>% filter(round(ab, 4) == .2375 & round(ba, 4) == .1875) -> ptg2

ptg2 %>% ggplot(aes(x = ac, y = bc, col = ca)) + geom_point(alpha = .5) + coord_fixed()
ptg2 %>% ggplot(aes(x = ac, y = cb)) + geom_point(alpha = .1) + coord_fixed()


  ## check bad cases

rmses %>%
  filter(increment == .025 & boost_low_alphas == T & rmse > 1) %>% pull(j) -> bad_js

rmses %>%
  filter(increment == .025 & boost_low_alphas == T & rmse < .05) %>% pull(j) -> good_js

apply(v_mat[bad_js, ]*ss[bad_js], 1, min)
apply(v_mat[good_js, ]*ss[good_js], 1, min)
# so tiny alphas is an issue, but not the only issue.


# points I need to make:
  ## the RMSE is low and gets lower as simulations increase.
  ## grid approach is more stable than the simulation when pivot probs are low -- how to show that?
  # I guess this is the idea: one massive simulation to get the truth, with 10M draws; then
  ## is there something better than just boosting up the tiny alphas? should I scale all of them?



rmses %>%
  ggplot(aes(x = increment, y = rmse)) +
  geom_violin(aes(group = increment))

# another view
rmses %>%
  ggplot(aes(x = increment, y = rmse, group = j)) +
  geom_line()

# sim_results %>%
#   ggplot(aes(x = n, y = value, col = pp)) +
#   geom_line(aes(group = pp)) +
#   geom_point() +
#   scale_x_log10() +
#   scale_y_log10()

# and then to compare the results via e.g. RMSE:
# I think I know how to do the correct normalization. and I think the problem will depend on s. let's show that, using s = 1/3?

# stopping now without

