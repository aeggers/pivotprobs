
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pivotprobs

<!-- badges: start -->

<!-- badges: end -->

`pivotprobs` is a package for computing the probability of pivot events
(election results where a single ballot can determine the winner) and
other election events.

It aims to offer both general methods that can handle any voting system
and belief distribution and specialized methods that handle important
cases efficiently.

## Installation

You can install development version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("aeggers/pivotprobs")
```

## Example

Lots of dependencies – you also need `rcdd` and `geometry`.

``` r
library(pivotprobs)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

alpha3 <- c(.4, .35, .25)*85
plurality3 <- plurality_election(k = 3)

sc_out <- plurality3 %>%
   election_event_probabilities(method = "sc", alpha = alpha3, tol = .1)
sc_out[["a_b"]]$integral
#> [1] 0.003483393
sc_out[["a_b"]]$seconds_elapsed
#> [1] 0.04988289
sc_out[["a_c"]]$integral
#> [1] 0.0003622554
sc_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.1887867

mc_out <- plurality3 %>%
   election_event_probabilities(method = "mc", alpha = alpha3, num_sims = 100000)
mc_out[["a_b"]]$integral
#> [1] 0.003511
mc_out[["a_b"]]$seconds_elapsed
#> [1] 0.2849569
mc_out[["a_c"]]$integral
#> [1] 0.000332
mc_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.8242304

ev_out <- plurality3 %>%
  election_event_probabilities(method = "ev", alpha = alpha3)
ev_out[["a_b"]]$integral
#> [1] 0.003490655
ev_out[["a_b"]]$seconds_elapsed
#> [1] 0.007524014
ev_out[["a_c"]]$integral
#> [1] 0.0003546949
ev_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.02365303
```
