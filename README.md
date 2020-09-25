
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

## Setup

``` r
library(pivotprobs)
library(tidyverse)
library(kableExtra) # for presentation
```

## Basic syntax for getting a P-matrix from Dirichlet parameters

Suppose we just want a P-matrix given Dirichlet parameters and an
electorate size.

### 3-candidate plurality

We will use the Eggers-Vivyan method, which is biased for elections with
more than three candidates but fine in the three-candidate case (and
very fast).

``` r
# define some Dirichlet parameters 
alpha3 <- c(.4, .35, .25)*85
electorate_size <- 15000

plurality_election(k = 3, n = electorate_size) %>% 
  election_event_probs(method = "ev", alpha = alpha3) %>% 
  combine_P_matrices()
#>              [,1]         [,2]         [,3]
#> [1,] 5.127133e-04 2.364633e-05 2.327103e-04
#> [2,] 1.359716e-05 4.926149e-04 2.327103e-04
#> [3,] 1.359716e-05 2.364633e-05 7.448696e-05
```

### 3-candidate IRV

We will use the Eggers-Nowacki method, which is validated and very fast.

``` r
# define some Dirichlet parameters 
alpha6 <- c(.3, .05, .2, .15, .1, .2)*85
electorate_size <- 15000

irv_election(n = electorate_size) %>% 
  election_event_probs(method = "en", alpha = alpha6) %>% 
  combine_P_matrices()
#>              [,1]         [,2]         [,3]         [,4]         [,5]
#> [1,] 4.756231e-04 4.756231e-04 9.725688e-05 3.237750e-06 4.384982e-04
#> [2,] 1.109281e-04 1.077429e-04 5.765073e-04 5.765073e-04 8.111857e-05
#> [3,] 8.725784e-05 9.044304e-05 4.481390e-08 9.406394e-05 1.541923e-04
#>              [,6]
#> [1,] 0.0002387696
#> [2,] 0.0002808471
#> [3,] 0.0001541923
```

## More detailed illustrations

### Election object

We always start with an election object, which is just a list containing
information about the voting system and the electorate size.

``` r
plurality_5 <- plurality_election(n = 5000, k = 5) # 5-candidate plurality
borda <- positional_election(n = 5000, s = .5) # 3-candidate Borda count
irv <- irv_election(n = 5000) # 3-candidate IRV
irv_borda <- irv_election(n = 5000, s = .5) # 3-candidate IRV with first round Borda count 
ky <- kemeny_young_election(n = 5000) # kemeny-young election 
```

### Computing event probabilities

For each of these, we can compute election event probabilities by
several methods using the `election_event_probs()` function.

(Pick up here.)

``` r
### 3 candidate plurality case ### 
# define the election object 
plurality3 <- plurality_election(k = 3, n = 15000) 
# this is the necessary first argument to election_event_probs().
# it is just a list containing the key attributes of the voting system.
# those attributes include "election events", each of which is associated 
# with a set of conditions under which it happens and a matrix showing how 
# a single ballot can affect the outcome at that event. 
# you must also specify the electorate, n.

# a bit of inspection 
plurality3$events$i_j$P
#>      [,1] [,2] [,3]
#> [1,]    1    0    1
#> [2,]    0    1    0
#> [3,]    0    0    0
plurality3$n
#> [1] 15000

# define some Dirichlet parameters 
alpha3 <- c(.4, .35, .25)*85

# compute pivot event probabilities using the Eggers-Vivyan method 
plurality3 %>% 
  election_event_probs(method = "ev", alpha = alpha3) -> pps

# inspect pivot probabilities 
tibble(event = names(pps), prob = pps %>% map("integral")) %>% 
  kbl(caption = "Pivot event probabilities computed by Eggers-Vivyan method") %>% 
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Pivot event probabilities computed by Eggers-Vivyan method

</caption>

<thead>

<tr>

<th style="text-align:left;">

event

</th>

<th style="text-align:left;">

prob

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

a\_b

</td>

<td style="text-align:left;">

0.000232710301369787

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c

</td>

<td style="text-align:left;">

2.36463262588568e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a

</td>

<td style="text-align:left;">

0.000232710301369787

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c

</td>

<td style="text-align:left;">

1.35971556490655e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a

</td>

<td style="text-align:left;">

2.36463262588568e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b

</td>

<td style="text-align:left;">

1.35971556490655e-05

</td>

</tr>

</tbody>

</table>

``` r

# the result includes P matrices for each event
pps$a_b$P 
#>      [,1] [,2] [,3]
#> [1,]    1    0    1
#> [2,]    0    1    0
#> [3,]    0    0    0
pps$a_c$P
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
#> [2,]    0    0    0
#> [3,]    0    0    1

# a convenience method to combine these event-specific P matrices according to event probabilities to get the full P matrix:  
pps %>% combine_P_matrices()
#>              [,1]         [,2]         [,3]
#> [1,] 5.127133e-04 2.364633e-05 2.327103e-04
#> [2,] 1.359716e-05 4.926149e-04 2.327103e-04
#> [3,] 1.359716e-05 2.364633e-05 7.448696e-05
```

Three-candidate IRV:

``` r
### 3 candidate IRV case ### 
# define the election object 
irv <- irv_election(n = 15000) 

# the Dirichlet parameters 
alpha6 <- c(.3, .05, .2, .15, .1, .2)*85

# compute pivot event probabilities using the Eggers-Nowacki method 
irv %>% 
  election_event_probs(method = "en", alpha = alpha6) -> pps

# inspect pivot probabilities 
tibble(event = names(pps), prob = pps %>% map("integral")) %>% 
  kbl(caption = "Pivot event probabilities computed by Eggers-Nowacki method") %>% 
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Pivot event probabilities computed by Eggers-Nowacki method

</caption>

<thead>

<tr>

<th style="text-align:left;">

event

</th>

<th style="text-align:left;">

prob

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

a\_b

</td>

<td style="text-align:left;">

9.98642894042279e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_b|ab

</td>

<td style="text-align:left;">

3.91826056180126e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_b|cb

</td>

<td style="text-align:left;">

3.91835133743872e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_b|ac

</td>

<td style="text-align:left;">

2.22804805450754e-08

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c

</td>

<td style="text-align:left;">

4.70095636226301e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c|ac

</td>

<td style="text-align:left;">

7.20246814001112e-12

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c|bc

</td>

<td style="text-align:left;">

2.45732301218642e-10

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c|ab

</td>

<td style="text-align:left;">

1.37622549172375e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a

</td>

<td style="text-align:left;">

9.98642894042279e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a|ba

</td>

<td style="text-align:left;">

3.91826056180126e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a|bc

</td>

<td style="text-align:left;">

3.91835133743872e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a|ca

</td>

<td style="text-align:left;">

2.22804805450754e-08

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c

</td>

<td style="text-align:left;">

1.59259776505929e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c|bc

</td>

<td style="text-align:left;">

7.96005580659376e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c|ac

</td>

<td style="text-align:left;">

9.30758803536643e-07

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c|ba

</td>

<td style="text-align:left;">

9.97823499940729e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a

</td>

<td style="text-align:left;">

4.70095636226301e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a|ca

</td>

<td style="text-align:left;">

7.20246814001112e-12

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a|cb

</td>

<td style="text-align:left;">

2.45732301218642e-10

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a|ba

</td>

<td style="text-align:left;">

1.37622549172375e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b

</td>

<td style="text-align:left;">

1.59259776505929e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b|cb

</td>

<td style="text-align:left;">

7.96005580659376e-06

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b|ca

</td>

<td style="text-align:left;">

9.30758803536643e-07

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b|ab

</td>

<td style="text-align:left;">

9.97823499940729e-05

</td>

</tr>

</tbody>

</table>

``` r

# make P matrix 
pps %>% combine_P_matrices()
#>              [,1]         [,2]         [,3]         [,4]         [,5]
#> [1,] 4.756231e-04 4.756231e-04 9.725688e-05 3.237750e-06 4.384982e-04
#> [2,] 1.109281e-04 1.077429e-04 5.765073e-04 5.765073e-04 8.111857e-05
#> [3,] 8.725784e-05 9.044304e-05 4.481390e-08 9.406394e-05 1.541923e-04
#>              [,6]
#> [1,] 0.0002387696
#> [2,] 0.0002808471
#> [3,] 0.0001541923
```

Other ways to compute pivot probabilities for the plurality-Dirichlet
case:

``` r
# SimplicialCubature: adaptive integration of election distribution
sc_out <- plurality3 %>%
   election_event_probs(method = "sc", alpha = alpha3, tol = .01)
tibble(event = names(sc_out), prob = sc_out %>% map("integral")) %>% 
  kbl(caption = "Election events in plurality (including non-pivot events)") %>% 
  kable_classic(full_width = F)
```

<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>

<caption>

Election events in plurality (including non-pivot events)

</caption>

<thead>

<tr>

<th style="text-align:left;">

event

</th>

<th style="text-align:left;">

prob

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

a\_

</td>

<td style="text-align:left;">

0.691496481272913

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_b

</td>

<td style="text-align:left;">

0.000232294171740642

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_bc

</td>

<td style="text-align:left;">

2.15985792933378e-08

</td>

</tr>

<tr>

<td style="text-align:left;">

a\_c

</td>

<td style="text-align:left;">

2.36393583466928e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_

</td>

<td style="text-align:left;">

0.288378500983843

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_a

</td>

<td style="text-align:left;">

0.000232206578211485

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_ac

</td>

<td style="text-align:left;">

2.1589403033054e-08

</td>

</tr>

<tr>

<td style="text-align:left;">

b\_c

</td>

<td style="text-align:left;">

1.36208326079811e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_

</td>

<td style="text-align:left;">

0.0194004433070118

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_a

</td>

<td style="text-align:left;">

2.36101763803129e-05

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_ab

</td>

<td style="text-align:left;">

2.1571058309774e-08

</td>

</tr>

<tr>

<td style="text-align:left;">

c\_b

</td>

<td style="text-align:left;">

1.35984212178784e-05

</td>

</tr>

</tbody>

</table>

``` r

#sc_out[["a_b"]]$integral
#sc_out[["a_b"]]$seconds_elapsed
#sc_out[["a_c"]]$integral
sc_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.386107

mc_out <- plurality3 %>%
   election_event_probs(method = "mc", alpha = alpha3, num_sims = 100000)
mc_out[["a_b"]]$integral
#> [1] 0.000233
mc_out[["a_b"]]$seconds_elapsed
#> [1] 0.2706101
mc_out[["a_c"]]$integral
#> [1] 2.453333e-05
mc_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.8887336

ev_out <- plurality3 %>%
  election_event_probs(method = "ev", alpha = alpha3)
ev_out[["a_b"]]$integral
#> [1] 0.0002327103
ev_out[["a_b"]]$seconds_elapsed
#> [1] 0.01159
ev_out[["a_c"]]$integral
#> [1] 2.364633e-05
ev_out %>% map("seconds_elapsed") %>% unlist() %>% sum()
#> [1] 0.04245305
```

We can Other methods:

``` r
# borda count 
positional_election(n = 5000, s = .5) %>% 
  election_event_probs(method = "mc", alpha = alpha6, num_sims = 500000) %>% 
  combine_P_matrices()
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] 0.06904854 0.06923930 0.06866348 0.06846918 0.06904500 0.06865994
#> [2,] 0.92901322 0.92880176 0.92940686 0.92958904 0.92898394 0.92937758
#> [3,] 0.00274976 0.00277046 0.00274118 0.00275330 0.00278258 0.00277400

# IRV with borda count first round
irv_election(n = 5000, s = .5) %>% 
  election_event_probs(method = "mc", alpha = alpha6, num_sims = 500000) %>% 
  combine_P_matrices()
#>            [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] 0.17569966 0.17569284 0.17475070 0.17473836 0.17568050 0.17473154
#> [2,] 0.82257064 0.82254894 0.82352198 0.82353186 0.82255882 0.82351016
#> [3,] 0.00272946 0.00275798 0.00272708 0.00272954 0.00276044 0.00275806

# kemeny-young (condorcet) 
#kemeny_young_election(n = 5000) %>% 
#  election_event_probs(method = "sc", alpha = alpha6, tol = .1) %>% 
#  combine_P_matrices()
```
