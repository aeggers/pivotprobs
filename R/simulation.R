#' Compute pivot probabilities from a matrix of simulated election results

ab_plurality_tie_for_first_from_sims <- function(sims, tol = .01){
  row_max <- apply(sims, 1, max)
  mean((sims[,1] == row_max | sims[,2] == row_max) & abs(sims[,1] - sims[,2]) < tol/2)/(tol*sqrt(2)) # Normalization here because the width of a channel where a and b is separated by tol holding fixed others' votes is sqrt(2)*tol
}

plurality_tie_probs_from_sims <- function(sims, tol = .01, cand_names = NULL, sep = ""){

  out <- list()
  if(is.null(cand_names)){cand_names <- letters[1:ncol(sims)]}
  for(i in 1:(ncol(sims)-1)){
    for(j in (i+1):ncol(sims)){
      out[[paste0(cand_names[i], sep, cand_names[j])]] = ab_plurality_tie_for_first_from_sims(cbind(sims[,c(i,j)], sims[,-c(i,j)]))
    }
  }
  out

}
