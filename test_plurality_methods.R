### for plurality: test grid based and analytical approx against simulation
N <- 100
n <- 100000
v_mat <- gtools::rdirichlet(N, alpha = c(5, 4, 3, 2, 1))
s <- 20

results
for(i in 1:nrow(v_mat)){
  alpha_vec <- as.vector(v_mat[i,]*s)
  names(alpha_vec) <- letters[1:5]
  # simulation based
  samp <- gtools::rdirichlet(n, alpha_vec)
  colnames(samp) <- names(alpha_vec)
  pps_sim_based <- pivotal.probabilities(samp, tol = .02)

  # grid based
  pps_grid_based <- plurality_pivot_probs_grid_based(alpha_vec, increment = .02)

  # analytical approximation
  pps_analytical <- pivotal.probabilities.analytical(alpha_vec, n = 1)

}

bind_rows(
  plurality_pivot_probs_grid_based(alpha_vec, increment = .05) %>% unlist(),
  plurality_pivot_probs_grid_based(alpha_vec, increment = .025) %>% unlist(),
  plurality_pivot_probs_grid_based(alpha_vec, increment = .01) %>% unlist(),
  plurality_pivot_probs_grid_based(alpha_vec, increment = .005) %>% unlist()
) -> gridders

gridders %>% mutate(increment = c(.05, .025, .01, .005)) -> gridders

big_sim <- gtools::rdirichlet(5000000, alpha_vec)
colnames(big_sim) <- names(alpha_vec)

bind_rows(
  pivotal.probabilities(big_sim, tol = .03) %>% unlist(),
  pivotal.probabilities(big_sim, tol = .02) %>% unlist(),
  pivotal.probabilities(big_sim, tol = .01) %>% unlist()
) %>% mutate(
  increment = c(.03, .02, .01)
) -> sims

bind_rows(
  pivotal.probabilities.analytical(alpha_vec, n = 1) %>% unlist(),
  pivotal.probabilities.analytical(alpha_vec, increments = 100, n = 1) %>% unlist()
) %>% mutate(
  increment = c(.005, .01)
) -> analyticals


bind_rows(gridders %>% mutate(type = "grid"),
          sims %>% mutate(type = "sims"),
          analyticals %>% mutate(type = "analytical")) %>%
  pivot_longer(cols = ab:de, names_to = "pp") -> for_plot

for_plot %>%
  ggplot(aes(x = increment, y = value, col = pp, linetype = type)) +
  geom_line(aes(group = interaction(pp, type))) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

# so in a few cases, I find that they all closely agree, but often the gridded approach gets closer to the simulation "truth".
# TODO: Can test this across cases as in the AV paper.
# TODO: also check on speed. maybe use those metrics they use in numerical integration methods.


