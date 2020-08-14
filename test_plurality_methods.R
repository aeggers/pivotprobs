### for plurality: test grid based and analytical approx against simulation
set.seed(123)
J <- 100
v_mat <- gtools::rdirichlet(J, alpha = c(5, 4, 3, 2))
ss <- runif(J, min = 10, max = 50)

rmses <- tibble(j = NULL, increment = NULL, rmse = NULL)

# v_mat <- gtools::rdirichlet(J, alpha = c(3,2,2,2,1,3))
N_for_sim_analysis <- 1000000
Ns <- c(N_for_sim_analysis)


for(j in 1:J){
  cat(j, ": ", sep = "")
  v_vec <- as.numeric(v_mat[j,])
  s <- ss[j]
  alpha_vec <- v_vec*s

  # simulation-based answers
  cat(" s ")
  df <- tibble(N = Ns)
  df %>%
    mutate(sim_result = pmap(., plurality_pivot_probs_simulation_based, alpha = alpha_vec)) -> df2

  df2 %>%
    unnest_longer(col = sim_result, values_to = "value", indices_to = "pp") -> sim_results

  # grid-based results
  cat(" g ")
  df <- tibble(increment = (2:5)/100) %>%
    mutate(gb_result = pmap(., plurality_pivot_probs_grid_based, alpha_vec = alpha_vec, cand_names = letters[1:length(alpha_vec)])) %>%
    unnest_longer(col = gb_result, values_to = "value", indices_to = "pp") -> gb_results

  # plotting?
  expand_grid(
    sim_results %>% filter(N == N_for_sim_analysis),
    increment = c(.01, .05)) %>% mutate(type = "sim") %>%
  bind_rows(gb_results %>%
    mutate(type = "grid")) %>%
  ggplot(aes(x = increment, y = value, col = pp, linetype = type)) +
    geom_line(aes(group = interaction(pp, type))) +
    scale_y_log10() +
    scale_x_log10()

  cat(" r")
  gb_results %>%
    left_join(sim_results %>% filter(N == N_for_sim_analysis) %>% select(-N), by = "pp", suffix = c("_grid", "_sim")) -> for_rmses

  for_rmses %>%
    group_by(increment) %>%
    summarize(rmse = sqrt(sum((value_grid - value_sim)^2)/n()),
              rmse_norm = sqrt(sum((value_grid/sum(value_grid) - value_sim/sum(value_sim))^2)/n())) -> these_rmses

  cat(".\n")

  rmses <- bind_rows(rmses, these_rmses %>% mutate(j = j))
}

rmses %>%
  ggplot(aes(x = increment, y = rmse_norm)) +
  geom_violin(aes(group = increment)) -> p1

# another view
rmses %>%
  ggplot(aes(x = increment, y = rmse_norm, group = j)) +
  geom_line() -> p2

library(patchwork)
p1 + p2

rmses %>% group_by(increment) %>% summarize(mean(rmse), mean(rmse_norm))

# OK now a different analysis:
# compare analytical and fine-grid grid-based to simulations with bigger and bigger size -- should get closer and closer.

Ns <- c(10000, 100000, 1000000)

rmses <- tibble(N = NULL, j = NULL, rmse_grid = NULL, rmse_analytical = NULL)
for(j in 1:J){

  cat(j, ": ", sep = "")
  v_vec <- as.numeric(v_mat[j,])
  s <- ss[j]
  alpha_vec <- v_vec*s

  # simulation-based answers
  cat(" s ")
  df <- tibble(N = Ns)
  df %>%
    mutate(sim_result = pmap(., plurality_pivot_probs_simulation_based, alpha = alpha_vec)) -> df2

  df2 %>%
    unnest_longer(col = sim_result, values_to = "value", indices_to = "pp") -> sim_results

  # grid-based results
  cat(" g ")
  df <- tibble(increment = 1/100) %>%
    mutate(gb_result = pmap(., plurality_pivot_probs_grid_based, alpha_vec = alpha_vec, cand_names = letters[1:length(alpha_vec)])) %>%
    unnest_longer(col = gb_result, values_to = "value", indices_to = "pp") -> gb_results

  # analytical results
  cat(" a ")
  alpha.vec <- alpha_vec
  names(alpha.vec) <- letters[1:length(alpha.vec)]
  dfa <- tibble(n = 1) %>%
    mutate(a_result = pmap(., pivotal.probabilities.analytical, alpha.vec = alpha.vec)) %>%
    unnest_longer(col = a_result, values_to = "value", indices_to = "pp") -> a_results

  cat(" r")
  sim_results %>%
    left_join(gb_results %>% select(-increment), by = "pp", suffix = c("_sim", "_grid")) %>%
    left_join(a_results %>% select(-n), by = "pp") %>%
    rename(value_analytical = value) -> for_rmses

  for_rmses %>%
    group_by(N) %>%
    summarize(rmse_grid = sqrt(sum((value_grid - value_sim)^2)/n()),
              rmse_analytical = sqrt(sum((value_analytical - value_sim)^2)/n())) -> these_rmses

  cat(".\n")

  rmses <- bind_rows(rmses, these_rmses %>% mutate(j = j))
}

rmses %>%
  pivot_longer(cols = c(rmse_grid, rmse_analytical), names_to = "type", values_to = "rmse", names_prefix = "rmse_") -> rmse_longer

rmse_longer %>%
  ggplot(aes(x = N, y = rmse, fill = type, group = interaction(type, N))) +
  geom_violin() +
  scale_x_log10()

# another approach

rmses %>%
  ggplot(aes(x = rmse_analytical, y = rmse_grid)) +
  geom_point() +
  # expand_limits(x = 0, y = 0) +
  scale_x_log10() +
  scale_y_log10() +
  expand_limits(x = .001, y = .001) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(. ~ N)

# or what about this -- seems good.
rmses %>%
  ggplot(aes(x = rmse_analytical, y = rmse_grid, group = factor(j))) +
  geom_point(aes(col = factor(N))) +
  geom_line(alpha = .25) +
  # expand_limits(x = 0, y = 0) +
  scale_x_log10() +
  scale_y_log10() +
  expand_limits(x = .001, y = .001) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0)


# yet another approach
rmse_longer %>%
  ggplot(aes(x = N, y = rmse, group = j)) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(. ~ type)

rmse_longer %>%
  ggplot(aes(x = N, y = rmse, group = j)) +
  geom_line() +
  scale_x_log10() +
  facet_wrap(. ~ type)


## Impressions from this: the grid approach seems to perform well generally: as the number of simulations gets large, the simulation result tends to get closer and closer to the grid result with .01 increments.
# It does not perform as well as the analytical result, which uses an approximation -- that approach goes faster and has lower mean RMSE.

rmses %>% group_by(N) %>% summarize(mean(rmse_grid), mean(rmse_analytical))

# If I have designed it correctly, the grid RMSE might eventually get better than the analytical version, because it does not include the independence assumption.

# actually let's check that for one case?

v_vec <- c(.35, .3, .27, .07)
s <- 30
alpha_vec <- v_vec*s
names(alpha_vec) <- letters[1:4]

df <- tibble(increment = c(1/100, 1/500, 1/1000)) %>%
  mutate(gb_result = pmap(., plurality_pivot_probs_grid_based, alpha_vec = alpha_vec, cand_names = letters[1:length(alpha_vec)])) %>%
  unnest_longer(col = gb_result, values_to = "value", indices_to = "pp") -> gb_results

alpha.vec <- alpha_vec
names(alpha.vec) <- letters[1:length(alpha.vec)]
dfa <- tibble(n = 1) %>%
  mutate(a_result = pmap(., pivotal.probabilities.analytical, alpha.vec = alpha.vec)) %>%
  unnest_longer(col = a_result, values_to = "value", indices_to = "pp") -> a_results

df <- tibble(N = c(50000, 500000, 5000000))
df %>%
  mutate(sim_result = pmap(., plurality_pivot_probs_simulation_based, alpha = alpha_vec)) %>%
  unnest_longer(col = sim_result, values_to = "value", indices_to = "pp") -> sim_results

# so now I want to know: is there a point at which the grid based is closer than the analytical?

a_results %>%
  select(-n) %>%
  mutate(type = "analytical") %>%
  bind_rows(gb_results %>% filter(increment == .01) %>% select(-increment) %>% mutate(type = "grid")) %>%
  left_join(sim_results %>% rename(value_sim = value), by = "pp") -> for_plot

for_plot %>%
  ggplot(aes(x = value_sim, y = value, col = pp, pch = type)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(. ~ N)
# It's hard to tell from that.



# let's look at the bad cases and good cases
rmses %>% filter(N > 500000 & rmse_grid > .05) %>% pull(j) -> bad_js
rmses %>% filter(N > 500000 & rmse_grid < .005) %>% pull(j) -> good_js
v_mat[bad_js, ]
ss[bad_js]
v_mat[good_js, ]
ss[good_js]

# The issue seems to relate to small values in the v_vec -- so I think this is again about getting close to the edge.

# let's look at a case.
# j = 2 is an example with a very small alpha value

j <- 2
v_vec <- v_mat[j,]
s <- ss[j]
(alpha_vec <- v_vec*s)
names(alpha_vec) <- letters[1:4]

# and we get a very different value for ab via analytical and grid_based
pps_a <- pivotal.probabilities.analytical(alpha_vec, n = 1)
pps_g <- plurality_pivot_probs_grid_based(alpha_vec, increment = .01)
tibble(pp = names(pps_a), analytical = pps_a %>% unlist(), grid_based = pps_g %>% unlist()) -> comparison

comparison %>% ggplot(aes(x = analytical, y = grid_based, col = pp)) +
  geom_point() +
  coord_fixed() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0)

# only the pivot probs not involving d (the one with the low alpha) are wrong. interesting. what if we just boosted d up to 1 for those cases?

# very different.

# now just confirm I can get the ab case hand:
ptg <- plurality_tie_grid(increment = .01, k = 4)
# now just compute the ab piv prob
ptg %>%
  mutate(dens = gtools::ddirichlet(as.matrix(.), alpha_vec)) -> ptg2
sum(ptg2$dens)*(.01^2)
# there it is.

# now I think the reason is that, with such a low alpha in the mix, the edge pieces of the grid get huge densities and the others get very little.
# here is a simple example of numerical integration going wrong with Dirichlet:
incr <- .001
xs <- seq(from = incr/2, to = 1 - incr/2, by = incr)
sum(dbeta(xs, 1, .1))*incr # should be 1, but isn't.
sum(dbeta(xs, 1, 2))*incr # is 1
# so when we are computing the a and b pivot probability, if d's alpha is really small we end up with a pivot prob that is too small.


# compare the distribution with this s and a higher s:
ptg %>%
  mutate(dens = gtools::ddirichlet(as.matrix(.), alpha_vec),
         dens2 = gtools::ddirichlet(as.matrix(.), alpha_vec*5)) -> for_plot

for_plot %>% summarize(sum(dens), sum(dens2))

for_plot %>%
  pivot_longer(cols = c(dens, dens2)) %>%
  ggplot(aes(x = value, col = name, fill = name)) +
  scale_x_log10() +
  geom_density(alpha = .25)

# and we can see this issue here

ptg2 %>% filter(dens > 100)
# we get big spiky densities for grid points where the component with the small alpha is supposed to get a very small result.
# but what can I do about this?

# what about boosting alpha up to 1 when it is excluded -- does that work?
alpha_vec2 <- alpha_vec
alpha_vec2[alpha_vec2 < 1] <- 1
ptg %>%
  mutate(dens = gtools::ddirichlet(as.matrix(.), alpha_vec),
         dens2 = gtools::ddirichlet(as.matrix(.), alpha_vec2)) -> check_fudge

sum(check_fudge$dens)*(.01^2)
sum(check_fudge$dens2)*(.01^2)
# well that works for one case -- but what could go wrong?
# seems like that improves things, though it's unsatisfyingly arbitrary.

# overall, the grid approach works (with this caveat about low alphas).
# now how about turning to the positional case? do we want to boost alphas there too?
# maybe set up a thing showing whether it works better or worse with boosted alphas.




## There is a general pattern that the grid based approach doesn't do very well when the alpha values are really small and some values are small. We were getting grid points very close to the edge of the simplex. If the corresponding alpha component is close to zero, then the density at this grid point will be huge. The grid points close to the edge were happening in the "last" column, 1 - everything else. To cover it up we require all points to be effective_increment/2 away from the edge.

# based on 100 runs:
# non-monotonic -- smallest at .03, largest at .04!
# this is probably about how the grid points lie in better or worse places in the space.
# nothing obviously different about the v_vecs that make bigger RMSEs.
# I'd like to know:
  ## is it the same with other values of s?
  ## do the sequences look like I expect?
  ## are the relative values good, but just not the absolute? the relative values are better, but not perfect. this could be because the increments affect how much volume is cut off.

# So is this good enough or not?
# Tends to be better with smaller increment, but not always. I think that's because the different increments imply cutting off different data.
# especially when there's a lot of density near the edge, a smaller grid size means getting more of that distortion.


# a bit of fun visualization

ptg <- plurality_tie_grid(increment = .001, k = 4)
colnames(ptg) <- c("a", "b", "c", "d")
ptg <- as_tibble(ptg)

# missing points:
# a = .31, b = .31, c = .07, d = .31
# I suspect it's because d is considered to be winning that one.
# this looks like a floating point issue.

ptg %>%
  mutate(x = a + .5*b,
         y = sqrt(3/4)*b,
         z = sqrt(3/4)*d
         # dens = gtools::ddirichlet(., alpha = c(.4, .35, .2, .05)*20),
         # col_ntile = dplyr::ntile(dens, n = 20)
         ) -> for_plot

#cols <- rainbow(20)
#rgl::plot3d(x = for_plot$x, y = for_plot$y, z = for_plot$z, col = cols[for_plot$col_ntile], type = "p")
rgl::plot3d(x = for_plot$x, y = for_plot$y, z = for_plot$z, xlim = c(1/3, 1), ylim = c(0,.5), zlim = c(0,.5), cex = 2, type = "p")






## earlier and less systematic work below.

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


