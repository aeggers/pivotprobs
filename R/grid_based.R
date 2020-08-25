# grid-based -- intended to supersede grid_based_pivot_probs.R

plurality_tie_grid <- function(increment = .025, k = 4){
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

grid_based_plurality_pivot_probs <- function(distribution = "dirichlet", normalizing_factor = 1, ...){

  args <- list(...)

  if(!distribution %in% c("dirichlet", "mvnorm")){
    stop("For `grid_based_plurality_pivot_probs(), the `distribution` argument must be `dirichlet` or `mvnorm`.")
  }
  if(distribution == "dirichlet"){
    if(is.null(args$alpha)){
      if(is.null(args$mu)){
        stop("For `grid_based_plurality_pivot_probs()` and `distribution` 'dirichlet', you must pass either `mu` or `alpha`.")
      }else{
        if(is.null(args$precision)){stop("For `grid_based_plurality_pivot_probs()` and `distribution` 'dirichlet', if you do not pass `alpha` you must pass `precision`.")}
        if(args$precision <= 0){stop("For `grid_based_plurality_pivot_probs()` and `distribution` 'dirichlet', `precision` must be non-negative.")}
        args$alpha <- args$mu*args$precision
      }
    }else{
      if(min(args$alpha) <= 0){
        stop("For `grid_based_plurality_pivot_probs()` and `distribution` 'dirichlet', all elements of `alpha` must be positive.")
      }
    }
    k <- length(args$alpha)

  }else if(distribution == "mvnorm"){
    if(is.null(args$mu) | is.null(args$sigma)){
      stop("For `grid_based_plurality_pivot_probs()` and `distribution` 'dirichlet' and `distribution` 'mvnorm', you must pass `mu` and `sigma`.")
    }
    k <- length(args$mu)
  }

  if(is.null(args$cand_names)){args$cand_names <- letters[1:k]}
  if(is.null(args$sep)){args$sep = ""}
  if(is.null(args$increment)){args$increment = .025}

  this_ptg_name <- paste0("k_", k, "_increment_", args$increment)
  if(this_ptg_name %in% names(TIE_GRID_STORE)){
    ptg <- TIE_GRID_STORE[[this_ptg_name]]
  }else{
    TIE_GRID_STORE[[this_ptg_name]] <<- ptg <- plurality_tie_grid(increment = args$increment, k = k)
  }

  grid_hypercube_volume <- plurality_grid_hypervolume_from_ptg(ptg)

  out <- list()
  for(i in 1:(length(args$cand_names)-1)){
    for(j in (i+1):length(args$cand_names)){
      pp_name <- paste0(args$cand_names[i], args$sep, args$cand_names[j])
      indices <- c(i,j,(1:k)[-c(i,j)])
      if(distribution == "dirichlet"){
        out[[pp_name]] <- (1/normalizing_factor)*grid_hypercube_volume*sum(gtools::ddirichlet(as.matrix(ptg), alpha = args$alpha[indices]))/sqrt(length(alpha))
      }else if(distribution == "mvnorm"){
        out[[pp_name]] <- (1/normalizing_factor)*grid_hypercube_volume*sum(mvtnorm::dmvnorm(as.matrix(ptg), mean = args$mu[indices], sigma = args$sigma[indices, indices]))
      }else{
        stop("The `distribution` argument must be 'dirichlet' or 'mvnorm'.")
      }
    }
  }
  out
}

plurality_grid_hypervolume_from_ptg <- function(ptg){
  # we get the product of the first k-2 unique distances -- could produce an error if a diagonal distance is less than an orthogonal distance
  # k-2 because the hyperplane is in k-2 dimensions.
  k <- ncol(ptg)
  i <- floor(nrow(ptg)/2) # picking a point from the middle of the grid matrix
  p <- as.numeric(ptg[i, ])
  m <- as.matrix(ptg[-i,])
  sdpmp <- unique(sort(distance_from_point_to_matrix_of_points(p, m)))
  prod(sdpmp[1:(k-2)]) # the product of the first k-2 unique distances
}

distance_between_points <- function(p1, p2){
  sqrt(sum((p1 - p2)^2))
}

distance_from_point_to_matrix_of_points <- function(p, m){
  stopifnot(length(p) == ncol(m))
  apply(m, 1, distance_between_points, p2 = p)

}


