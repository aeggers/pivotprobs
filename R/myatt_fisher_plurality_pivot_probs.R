incomplete_beta <- function(x,a,b){
  pbeta(x,a,b)*beta(a,b)
}

one_myatt_fisher_plurality_pivot_prob <- function(alpha){
  incomplete_beta_part <- incomplete_beta(1/3, alpha[3] + 1, sum(alpha[1:2]) + 1)
  (2^alpha[3])*incomplete_beta_part/(2^(sum(alpha) + 1))
}

myatt_fisher_plurality_pivot_probs <- function(alpha, cand_names = NULL, sep = ""){
  if(is.null(cand_names)){cand_names <- letters[1:3]}
  out <- list()
  for(i in 1:2){
    for(j in (i+1):3){
      out[[paste0(cand_names[i], sep, cand_names[j])]] <- one_myatt_fisher_plurality_pivot_prob(c(alpha[c(i,j)], alpha[-c(i,j)]))
    }
  }
  out
}

# this was the code previously tested
# posterior.pivotal.probability = function(pie.vec, s, which.ones = c(1,1,0), prior = 1){
#   pie = pie.vec[which(which.ones == 0)]
#   inc.beta.part = incomplete.beta(1/3, pie*s + prior, (1 - pie)*s + prior)
#   2^(pie*s)*inc.beta.part/2^(s + prior)
# }
