# grid-based -- intended to supersede grid_based_pivot_probs.R

grid_based_plurality_pivot_probs <- function(increment = .01, alpha = NULL, mu = NULL, precision = NULL, sigma = NULL, cand_names = NULL, sep = ""){

  if(!is.null(alpha) | (!is.null(mu) & !is.null(precision)) & is.null(sigma)){
    # this is Dirichlet
    distribution = "dirichlet"
    if(is.null(alpha)){alpha <- mu*precision}
    k <- length(alpha)
  }else if(!is.null(mu) & !is.null(sigma)){
    # this is logistic normal
    distribution = "logisticnormal"
    k <- length(mu)
  }

  if(is.null(cand_names)){cand_names <- letters[1:k]}

  this_ptg_name <- paste0("k_", k, "_increment_", increment)
  if(this_ptg_name %in% names(TIE_GRID_STORE)){
    ptg <- TIE_GRID_STORE[[this_ptg_name]]
  }else{
    TIE_GRID_STORE[[this_ptg_name]] <<- ptg <- plurality_tie_grid(increment = increment, k = k)
  }

  grid_hypercube_volume <- plurality_grid_hypervolume_from_ptg(ptg)

  out <- list()
  for(i in 1:(length(cand_names)-1)){
    for(j in (i+1):length(cand_names)){
      pp_name <- paste0(cand_names[i], sep, cand_names[j])
      these_indices <- c(i,j,(1:k)[-c(i,j)])
      if(distribution == "dirichlet"){
        out[[pp_name]] <- grid_hypercube_volume*sum(gtools::ddirichlet(as.matrix(ptg), alpha = alpha[these_indices]))/sqrt(k)
      }else if(distribution == "logisticnormal"){
        out[[pp_name]] <- grid_hypercube_volume*sum(dlogisticnormal(as.matrix(ptg), mu = mu[these_indices], sigma = sigma[these_indices, these_indices]))/sqrt(k)
      }
    }
  }
  out
}

#' @export
plurality_tie_grid <- function(k = 4, increment = .025){
  # this one is base R
  stopifnot(k >= 3)
  # beware using a small increment and high k because that will take a long time.
  increments <- round(1/increment, 0)
  effective_increment <- 1/increments
  standard_seq <- seq(0 + effective_increment/2, 1 - effective_increment/2, length = increments)

  to_expand <- as.data.frame(matrix(standard_seq, ncol = k - 2, nrow = length(standard_seq), byrow = F))

  the_grid <- expand.grid(to_expand)
  the_grid$V0 <- the_grid$V1
  the_grid$last <- 1 - apply(the_grid, 1, sum)
  the_grid <- the_grid[which(the_grid$last >= effective_increment/2),]
  # slight fudge to deal with floating point stuff, which was making slightly differently shaped grids depending on the increment.
  the_grid <- the_grid[which(the_grid$V0 > apply(the_grid, 1, max) - .00001),]
  the_grid[,c("V0", "V1", names(the_grid)[which(!names(the_grid) %in% c("V0", "V1", "last"))], "last")]

}

#sum_of_grid_point_dirichlet_densities <- function(tie_grid, alpha){
#  sum(gtools::ddirichlet(as.matrix(tie_grid), alpha)/sqrt(length(alpha)))
#}

TIE_GRID_STORE <- list()



plurality_grid_hypervolume_from_ptg <- function(ptg){
  # we get the product of the first k-2 unique distances -- could produce an error if a diagonal distance is less than an orthogonal distance, or if two distances are judged to be different even though they are the same
  # k-2 because the hyperplane is in k-2 dimensions.
  k <- ncol(ptg)
  i <- floor(nrow(ptg)/2) # picking a point from the middle of the grid matrix
  p <- as.numeric(ptg[i, ])
  m <- as.matrix(ptg[-i,])
  sdpmp <- unique(sort(round(distance_from_point_to_matrix_of_points(p, m), 10))) # rounding to avoid floating point issues
  # for reporting until no longer need:
  cat("Judging these to be the dimensions of the hypervolume: ", paste(sdpmp[1:(k-2)], collapse = (", ")), ".\n")
  prod(sdpmp[1:(k-2)]) # the product of the first k-2 unique distances
}

distance_between_points <- function(p1, p2){
  sqrt(sum((p1 - p2)^2))
}

distance_from_point_to_matrix_of_points <- function(p, m){
  stopifnot(length(p) == ncol(m))
  apply(m, 1, distance_between_points, p2 = p)
}


