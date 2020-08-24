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

one_pivot_prob_from_tie_grid_and_alpha_unnormalized <- function(tie_grid, alpha){
  sum(gtools::ddirichlet(as.matrix(tie_grid), alpha))
}

TIE_GRID_STORE <- list()

grid_based_plurality_pivot_probs <- function(distribution = "dirichlet", increment = .025, cand_names = NULL, sep = "", ...){

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

  if(is.null(cand_names)){cand_names <- letters[1:k]}

  this_ptg_name <- paste0("k_", k, "_increment_", increment)
  if(this_ptg_name %in% names(TIE_GRID_STORE)){
    ptg <- TIE_GRID_STORE[[this_ptg_name]]
  }else{
    TIE_GRID_STORE[[this_ptg_name]] <<- ptg <- plurality_tie_grid(increment, k)
  }

  effective_increment <- 1/round(1/increment, 0)
  out <- list()
  for(i in 1:(length(cand_names)-1)){
    for(j in (i+1):length(cand_names)){
      pp_name <- paste0(cand_names[i], sep, cand_names[j])
      indices <- c(i,j,(1:k)[-c(i,j)])
      normalizer <- effective_increment^(k-2)
      if(distribution == "dirichlet"){
        out[[pp_name]] <- sum(gtools::ddirichlet(as.matrix(tie_grid), alpha = args$alpha[indices]))*normalizer
      }else if(distribution == "mvrnorm"){
        out[[pp_name]] <- sum(mvtnorm::mvnorm(as.matrix(tie_grid), mu = args$mu[indices], sigma = args$sigma[indices, indices]))*normalizer
      }
    }
  }
  out
}


# I thought I might need to normalize, but it appears I don't?
# gr <- expand_grid(x = standard_seq, y = standard_seq) %>% filter(x + y < 1) %>% mutate(z = 1 - x - y)
# sum(gtools::ddirichlet(as.matrix(gr), alpha = rep(1, 3)))*increment^2
# I guess the point is that each of these points is actually a triangle? could I get the area of the volumes in the grid?
# TODO: look into this, perhaps once I have some results


