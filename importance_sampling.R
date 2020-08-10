## importance sampling tutorial

## set up a simple problem
## what is the probability of a tie when one candidate is very highly favored?

v_vec <- c(.8, .2)
s <- 20
N <- 1000000
alpha_vec <- v_vec*s

# analytical
answer <- pbeta(.505, alpha_vec[1], alpha_vec[2]) - pbeta(.5, alpha_vec[1], alpha_vec[2])

sims <- rbeta(N, alpha_vec[1], alpha_vec[2])

library(tidyverse)

# direct monte carlo
tibble(
  sim = sims,
  ab_tie = as.integer(sims > .5 & sims < .505),
  cum_nume = cumsum(ab_tie),
  ones = 1,
  cum_ones = cumsum(ones),
  cum_denom = cum_ones,
  cum_piv_prob = cum_nume/cum_denom) -> sim_version

# importance sampling version
q_alpha_vec <- c(.5, .5)*20 # was 20, so that we got some weird samples.

is_sims <- rbeta(N, q_alpha_vec[1], q_alpha_vec[2])

tibble(
  sim = is_sims,
  lr = dbeta(is_sims, alpha_vec[1], alpha_vec[2])/dbeta(is_sims, q_alpha_vec[1], q_alpha_vec[2]),
  ab_tie = as.integer(is_sims > .5 & is_sims < .505),
  cum_nume = cumsum(ab_tie*lr),
  ones = 1,
  cum_ones = cumsum(ones),
  cum_denom = cumsum(lr),
  cum_piv_prob = cum_nume/cum_denom) -> is_version

bind_rows(
  sim_version %>% mutate(method = "direct"),
  is_version %>% mutate(method = "IS")
) -> combined

combined %>%
  ggplot(aes(x = cum_ones, y = cum_piv_prob, col = method)) +
  geom_line(aes(group = method)) -> p

p +
  geom_hline(yintercept = answer, col = "red")

# interesting: they all give the right answer, but direct works better here.
# seems that IS version suffers from big drops Why would this happen? Probably we get a draw that would be really unlikely under the biased distribution and likely under the true distribution, so the denominator jumps way up.
# see what happened in this version

is_version %>%
  mutate(change_prob = cum_piv_prob - lag(cum_piv_prob)) %>%
  filter(cum_ones > 5000) %>%
  arrange(change_prob) %>%
  select(sim, lr, ab_tie, cum_ones, cum_piv_prob, change_prob) %>%
  head()

# so this suggests to me that I need to make the biased distribution more targeted to the piv prob

is_version %>%
  ggplot(aes(x = sim)) +
  geom_density()

# and actually, this is a trivial example but why don't we use a uniform distribution so that all the samples fit?

is_sims <- runif(N, min = .5, .505)

tibble(
  sim = is_sims,
  lr = dbeta(is_sims, alpha_vec[1], alpha_vec[2])/(1/.005),
  ab_tie = as.integer(is_sims > .5 & is_sims < .505),
  cum_nume = cumsum(ab_tie*lr),
  ones = 1,
  cum_ones = cumsum(ones),
  cum_denom = cumsum(lr),
  cum_piv_prob = cum_nume/cum_denom) -> is_version

bind_rows(
  sim_version %>% mutate(method = "direct"),
  is_version %>% mutate(method = "IS")
) -> combined

combined %>%
  filter(cum_ones < 20000 & cum_ones > 1000) %>%
  ggplot(aes(x = cum_ones, y = cum_piv_prob, col = method)) +
  geom_line(aes(group = method)) -> p2

p2 +
  geom_hline(yintercept = answer, col = "green", linetype = 2)

# this looks good except for the constant:
  ## no big jumps
  ## converging to something. what is the ratio converging to?

combined %>%
  pivot_wider(id_cols = cum_ones, names_from = method, values_from = cum_piv_prob) %>%
  mutate(ratio = IS/direct) %>%
  filter(cum_ones < 20000 & cum_ones > 1000) %>%
  ggplot(aes(x = cum_ones, y = ratio)) +
  geom_line()


# OK then how about just making a huge grid?

incr <- .01
the_seq <- seq(incr/2, 1 - incr/2, by = incr)
the_grid <- expand_grid(a = the_seq, b = the_seq) %>%
  filter(a + b <= 1) %>%
  mutate(c = 1 - a - b) %>%
  filter(a > 0 & b > 0 & c > 0)

the_grid %>%
  mutate(a = a + .5*b, b = sqrt(3/4)*b) %>%
  ggplot(aes(x = a, y = b)) +
  coord_fixed() +
  theme_void() +
  geom_point(size = .25, alpha = .1)
# takes a long time to plot

# this is the slow part, but only nec for the normalization
the_grid %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha = c(5,5,5))) -> the_grid_plus

# I think the simplex has an area of 1/2 in total.
normalizer <- .5/nrow(the_grid)

normalizer*sum(the_grid_plus$ddir_555)  ## and that's close to correct.

# so if we want to get the probability of someone winning
the_grid_plus %>%
  filter(a > b & a > c) %>%
  summarize(sum(ddir_555)*normalizer)

# fine.

# what about a piv prob?

# the number depends a lot on the tol chosen.
# ideally it wouldn't. maybe it wouldn't if we had more points.
# but this does not bode well for higher dimensions.

# the problem here was that I had not yet recognized that I could make a grid on the pivot probability hyperplane.
ab_tie_given_tol <- function(tol = .01){
  the_grid_plus %>%
    filter(a > c & b > c & abs(a - b) < tol/2) %>%
    summarize(sum(ddir_555)*normalizer)/tol
}

tibble(
  tol = seq(from = .001, to = .03, by = .00069),
  ab_tie_pp = map(tol, ab_tie_given_tol)
) %>% unnest(ab_tie_pp) -> out

# what's the actual answer?
incr <- .001
tibble(a = seq(1/3 + incr/2, 1/2, by = incr)) %>%
  mutate(b = a,
         c = 1 - 2*b) %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha = c(5,5,5))) %>%
  summarize(sum(ddir_555))*.001 -> answer

sims <- gtools::rdirichlet(1000000, alpha = c(5,5,5))

mean(sims[,1] > sims[,3] & sims[,2] > sims[,3] & abs(sims[,1] - sims[,2]) < tol/2)/tol -> sim_answer

out %>% ggplot(aes(x = tol, y = `sum(ddir_555) * normalizer`)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = as.numeric(answer), col = "red") +
  geom_hline(yintercept = as.numeric(sim_answer), col = "blue")

# why do we get the sawtooth pattern? because of the regularity of the grid and the integral being taken. the result wouldn't depend on the tol in this pattern if, instead of a grid, we had a big Dirichlet uniform sample.

## Can I show this first?

N <- 500000
the_grid <- gtools::rdirichlet(N, alpha = c(1,1,1))
colnames(the_grid) = c("a", "b", "c")
the_grid <- as_tibble(the_grid)

# this is very slow, and unnecessary. but it shows the randomness.
the_grid %>%
  mutate(a = a + .5*b, b = sqrt(3/4)*b) %>%
  ggplot(aes(x = a, y = b)) +
  coord_fixed() +
  theme_void() +
  geom_point(size = .25, alpha = .1)


# this is the slow part, but only nec for the normalization
the_grid %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha = c(5,5,5))) -> the_grid_plus

# I think the simplex has an area of 1/2 in total.
normalizer <- .5/nrow(the_grid)

normalizer*sum(the_grid_plus$ddir_555)  ## and that's close to correct.

# so if we want to get the probability of someone winning
the_grid_plus %>%
  filter(a > b & a > c) %>%
  summarize(sum(ddir_555)*normalizer)

tibble(
  tol = seq(from = .001, to = .03, by = .001),
  ab_tie_pp = map(tol, ab_tie_given_tol)
) %>% unnest(ab_tie_pp) -> out_samp

bind_rows(out_samp %>% mutate(type = "sample"), out %>% mutate(type = "grid")) %>%
  ggplot(aes(x = tol, y = `sum(ddir_555) * normalizer`, col = type)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = as.numeric(answer), col = "red") +
  geom_hline(yintercept = as.numeric(sim_answer), col = "blue")


# what we're doing so far is summing the dirichlet density over grid points within the tol window.
# we're saying the area is the standard grid size, and the density for each area is the density at the point.
# the problem is that the number of grid points in the key area depends too much (and in arbitrary ways) on the definition of the tol window.

# I have a vague idea that I could define points right on this surface and integrate that way.

# let's consider a simple case: tie for first between two of four candidates.
incr <- .005
expand_grid(
  a = seq(.25 + incr/2, .5 - incr/2, by = incr),
  c = seq(0 + incr/2, 1/3 - incr/2, by = incr)
) %>%
  mutate(b = a, d = 1 - a - b - c) %>%
  filter(d > 0 & d < 1 & a > c & a > d) %>%
  select(a, b, c, d) -> grid2

# plot this:
grid2 %>% mutate(x = a + .5*b, y = sqrt(3/4)*b, z = sqrt(3/4)*d) -> for_plot

rgl::plot3d(x = for_plot$x, y = for_plot$y, z = for_plot$z, type = "p")

# so we're in 2d, so incr^2 is the area.

alpha_vec <- c(.4, .35, .2, .05)*40
grid2 %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha_vec)) %>%
  summarize(area = sum(ddir_555)*incr^2) -> out_ab

grid2 %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha_vec[c(2,3,1,4)])) %>%
  summarize(area = sum(ddir_555)*incr^2) -> out_bc

grid2 %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha_vec[c(1,3,2,4)])) %>%
  summarize(area = sum(ddir_555)*incr^2) -> out_ac

grid2 %>%
  mutate(ddir_555 = gtools::ddirichlet(as.matrix(.), alpha_vec[c(1,4,2,3)])) %>%
  summarize(area = sum(ddir_555)*incr^2) -> out_ad

# let's compare to the simulation version
sims4 <- gtools::rdirichlet(10000000, alpha_vec)
tol <- .01
ab <- mean(abs(sims4[,1] - sims4[,2]) < tol/2 & sims4[,1] > sims4[,3] & sims4[,1] > sims4[,4])/tol
bc <- mean(abs(sims4[,2] - sims4[,3]) < tol/2 & sims4[,2] > sims4[,1] & sims4[,2] > sims4[,4])/tol
ac <- mean(abs(sims4[,1] - sims4[,3]) < tol/2 & sims4[,1] > sims4[,2] & sims4[,1] > sims4[,4])/tol
ad <- mean(abs(sims4[,1] - sims4[,4]) < tol/2 & sims4[,1] > sims4[,2] & sims4[,1] > sims4[,3])/tol
analytical <- unlist(c(out_ab, out_ac, out_ad, out_bc))
names(analytical) <- c("ab", "ac", "ad", "bc")
simulation_based <- c(ab, ac, ad, bc)
names(simulation_based) <- c("ab", "ac", "ad", "bc")
analytical/simulation_based
# so ab and ac agree closely, bc less so, and ad not very much.
# oh but hold up -- with more simulations it gets closer and closer. ad is very unlikely.
# OK so I think that DOES work. cool.

# now we see if we can do it for borda count.

# borda count. we have 4 degrees of freedom.
# this takes a long time with small incr, but only needs to happen once.
# OTOH we do need to compute exactly below, and that gets slow if incr is small
incr <- .01
s <- .5
expand_grid(
  ab = seq(0 + incr/2, 1 - incr/2, by = incr),
  ba = seq(0 + incr/2, 1 - incr/2, by = incr),
  ca = seq(0 + incr/2, 1 - incr/2, by = incr),
  cb = seq(0 + incr/2, 1 - incr/2, by = incr)
) %>%
  mutate(
    ac = .5*(1 - s*ba + (s-2)*ab - (1 + s)*ca - (1 - s)*cb),
    bc = 1 - ab - ac - ba - ca - cb,
    score_a = ab + ac + s*(ba + ca),
    score_c = ca + cb + s*(ac + bc)
  ) %>%
  filter(bc > 0 & bc < 1 & ac > 0 & ac < 1 & score_c < score_a) %>% # filtering on ac not necessary I think.
  select(ab, ac, ba, bc, ca, cb) -> grid_bc

# check that's it's correct -- yes it is now correct (was not before).

grid_bc %>% mutate(
  score_a = ab + ac + s*(ba + ca),
  score_b = ba + bc + s*(ab + cb),
  score_c = ca + cb + s*(ac + bc),
  score_diff_ab = score_a - score_b,
  score_diff_ac = score_a - score_c
) %>% ggplot(aes(x = score_diff_ab, y = score_diff_ac)) + geom_point(alpha = .1)
# okay that works. (was not getting )

# 28570 rows with .025 incr. 1.14M rows with .01 incr.
# so to get a piv prob:
v_vec <- c(12, 4, 7, 6, 3, 10)
v_vec <- v_vec/sum(v_vec)
s <- 20
alpha_vec <- v_vec*s

sims <- gtools::rdirichlet(10000000, alpha = alpha_vec)
# how often do a and b tie for first in Borda count?
pps <- positional_piv_probs_simulation(sims)

grid_bc %>%
  mutate(ddir_xxx = gtools::ddirichlet(as.matrix(.), alpha_vec)) %>%
  summarize(area = sum(ddir_xxx)*incr^4) -> out_ab

grid_bc %>%
  mutate(ddir_xxx = gtools::ddirichlet(as.matrix(.), alpha_vec[c(2,1,5,6,3,4)])) %>%
  summarize(area = sum(ddir_xxx)*incr^4) -> out_ac

grid_bc %>%
  mutate(ddir_xxx = gtools::ddirichlet(as.matrix(.), alpha_vec[c(4,3,6,5,1,2)])) %>%
  summarize(area = sum(ddir_xxx)*incr^4) -> out_bc

pps_analytical = unlist(c(out_ab, out_ac, out_bc))
names(pps_analytical) = c("ab", "ac", "bc")

pps_analytical/unlist(pps)

# so they are all close to twice as high, but not exactly.
# sources of error:
# sample size in the simulation.
# grid size in numerical version.
# so, with smaller increments/finer grids they get closer to 2.
# but this grid is very fine-grained and it takes a long time to compute.
# are the absolute levels wrong because I am getting the geometry wrong?

# are the hypercubes really size incr^4? for this to be true, we would need the grid points to be a distance of incr from each other in each dimension. is that the case?

library(patchwork)
# you can switch cb and bc and looks similar.
grid_bcx <- grid_bc %>% filter(ab == .225 & ba == .225)
p1 <- grid_bcx %>% ggplot(aes(x = ac, y = bc)) + geom_point(alpha = .1)
# interesting -- ac and bc are also in a grid; it's just slanty.
p2 <- grid_bcx %>% ggplot(aes(x = ab, y = ba)) + geom_point(alpha = .1)
p3 <- grid_bcx %>% ggplot(aes(x = ca, y = cb)) + geom_point(alpha = .1)
p1 + p2 + p3

# it appears that the increment is not as large in the bc/ac grid: about 50 dots in a .05x.05 space, vs 25 for the cb/ca grid. it seems to be in the bc dimension only that we pack it in.
# I can dimly see how this might happen. it only affects the ratios of course.

## So next steps: assess the method compared to direct MC. trading off speed and accuracy. need to consider edge cases, where I think the numerical approach will do better than direct Monte Carlo.

# need to stop on this.




# would be useful to have an adaptIntegrate thing here -- do I now understand enough about the geometry that I can do this?



# but I do think that the step sizes not being uniform is an issue
# can we deal with that?

# ac and bc take on lots of values.

# if we simplify: suppose it was a grid of a and b, and then c takes on various values. then in summing up we could go through the grid by values of a and b; the volumes would be incr^2 \times vertical (c) discrep. so summing up all of these hypercubes would take a long time

# suppose we have three shares, we have a gridded out, a condition on b, and c is just 1 - a - b. So this is a line on the unit simplex. if I can get the density at each point on the grid, then knowing it's a line I need to also get the distance between points, which depends on the a increments *and* the angle of the line.

# In my 4D case this was pretty straightforward because it yielded a perfectly regular grid, though maybe there's an issue with normalizing. In 6D I feel completely lost. Leave it there for now.

# In my sleep (Sunday night) I thought I had this resolved. One thought was that it is enough to work out the volumes once, so it's okay if the process is slow.

# The idea I had was to arrange the data so that I basically have gridded slices through a volume. But

grid_bc %>% arrange(ac, bc, ab, ba, ca, cb) %>% head()

library(patchwork)
p1 <- grid_bc %>% ggplot(aes(x = ac, y = bc)) + geom_point(alpha = .1)
# interesting -- ac and bc are also in a grid; it's just slanty.
p2 <- grid_bc %>% ggplot(aes(x = ab, y = ba)) + geom_point(alpha = .1)
p3 <- grid_bc %>% ggplot(aes(x = ca, y = cb)) + geom_point(alpha = .1)

p1 + p2 + p3

# you can switch cb and bc and looks similar.
grid_bcx <- grid_bc %>% filter(ab == .2625 & ba == .2625)
p1 <- grid_bcx %>% ggplot(aes(x = ac, y = cb)) + geom_point(alpha = .1)
# interesting -- ac and bc are also in a grid; it's just slanty.
p2 <- grid_bcx %>% ggplot(aes(x = ab, y = ba)) + geom_point(alpha = .1)
p3 <- grid_bcx %>% ggplot(aes(x = ca, y = bc)) + geom_point(alpha = .1)

p1 + p2 + p3


# so that's very promising.
# think of it as: for each

# apparently unrelated answers and I don't know why.
# a key reason I think is that this is not actually on a grid.
# these points are on the surface of interest, but they could be any distance apart.
grid_bc %>% filter(ba > .1 & ba < .11) %>% ggplot(aes(x = ab, y = ac)) + geom_point(alpha = .1)
grid_bc %>% ggplot(aes(x = ab, y = ac)) + geom_point(alpha = .1)

grid_bc %>% ggplot(aes(x = ab, y = ca)) + geom_point(alpha = .1)

grid_bc %>% filter(ab == .2 - incr/2, ca == .1 - incr/2, ba == .2 - incr/2) %>% ggplot(aes(x = cb, y = ac)) + geom_point(alpha = .1)

ggplot(aes(x = ab, y = ca)) + geom_point(alpha = .1)
