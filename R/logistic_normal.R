## logistic normal

# this is not correct -- see http://people.csail.mit.edu/tomasz/papers/huang_hln_tech_report_2006.pdf
rlogisticnormal <- function(n, mean, sigma){
  #draw from multivariate normal
  draws <- mvtnorm::rmvnorm(n, mean = mean, sigma = sigma)
  # logit transformation
  # softmax -- I'm not sure if this is correct
  ed <- exp(draws)
  ed/apply(ed, 1, sum)

}

# follows from above.
dlogisticnormal <- function(x, mean, sigma){

}
