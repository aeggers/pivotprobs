#'
#' @export
plurality_pivot_probs <- function(method = "eggers-vivyan", ...){

  args <- list(...)

  if(method == "eggers-vivyan"){

    # must pass either mu (length k) and precision (positive) or alpha (length k, positive)
    if(is.null(args$alpha)){
      if(is.null(args$mu)){
        stop("For method eggers-vivyan, you must pass either `mu` or `alpha`.")
      }else{
        if(is.null(args$precision)){stop("For method eggers-vivyan, if you do not pass `alpha` you must pass `precision`.")}
        if(args$precision <= 0){stop("For method eggers-vivyan, `precision` must be non-negative.")}
        args$alpha <- args$mu*args$precision
      }
    }else{
      if(min(args$alpha) <= 0){
        stop("For method eggers-vivyan, all elements of `alpha` must be positive.")
      }
    }

    eggers_vivyan_plurality_pivot_probs(args$alpha, cand_names = args$cand_names, sep = args$sep)

  }else if(method == "myatt-fisher"){

    # Computes exact relative pivot probabilities for three-candidate races using the method of Fisher and Myatt (2017 WP)
    # User must pass either mu (length 3) and precision (positive) or alpha (length 3, positive)
    if(is.null(args$alpha)){
      if(is.null(args$mu)){
        stop("For method myatt-fisher, you must pass either `mu` or `alpha`.")
      }else{
        if(length(args$mu) != 3){
          stop("For method myatt-fisher, `mu` must be length 3.")
        }
        if(is.null(args$precision)){stop("For method myatt-fisher, if you do not pass `alpha` you must pass `precision`.")}
        if(args$precision <= 0){stop("For method myatt-fisher, `precision` must be non-negative.")}
        args$alpha <- args$mu*args$precision
      }
    }else{
      if(length(args$alpha) !=3 ){
        stop("For method myatt-fisher, `alpha` must be length 3.")
      }else if(min(args$alpha) <= 0){
        stop("For method myatt-fisher, all elements of `alpha` must be positive.")
      }
    }

    myatt_fisher_plurality_pivot_probs(args$alpha, cand_names = args$cand_names, sep = args$sep)

  }else if(method == "simplicial-cubatature"){

    # Pivot probability calculation via numerical integration of the relevant facets.
    # Must pass a function name and parameters

    if(is.null(args$distribution)){
      stop("For method `simplicial-cubature` you must pass a `distribution` argument.")
    }
    if(!args$distribution %in% c("dirichlet", "mvnorm")){
      stop("For method `simplicial-cubature`, the `distribution` argument must be `dirichlet` or `mvnorm`.")
    }
    if(args$distribution == "dirichlet"){
      if(is.null(args$alpha)){
        if(is.null(args$mu)){
          stop("For `method` 'simplicial-cubature' and `distribution` 'dirichlet', you must pass either `mu` or `alpha`.")
        }else{
          if(is.null(args$precision)){stop("For `method` 'simplicial-cubature' and `distribution` 'dirichlet', if you do not pass `alpha` you must pass `precision`.")}
          if(args$precision <= 0){stop("For `method` 'simplicial-cubature' and `distribution` 'dirichlet', `precision` must be non-negative.")}
          args$alpha <- args$mu*args$precision
        }
      }else{
        if(min(args$alpha) <= 0){
          stop("For `method` 'simplicial-cubature' and `distribution` 'dirichlet', all elements of `alpha` must be positive.")
        }
      }

      simplicial_cubature_plurality_pivot_probs_dirichlet(alpha = args$alpha, ...)

    }else if(args$distribution == "mvnorm"){
      if(is.null(args$mu) | is.null(args$sigma)){
        stop("For `method` 'simplicial-cubature' and `integrand_function` 'mvnorm', you must pass `mu` and `sigma`.")
      }

      simplicial_cubature_plurality_pivot_probs_mvnorm(mu = args$mu, sigma = args$sigma, ...)
    }

  }else if(method == "simulation"){
    # you can pass the matrix of simulations, or
    if(is.null(args$sims)){
      if(is.null(args$distribution)){
        stop("For `method` 'simulation' you must pass either `sims` or a `distribution` argument.")
      }
      if(! args$distribution %in% c("dirichlet", "mvnorm")){
        stop("For `method` 'simulation' the `distribution` argument must be 'dirichlet' or 'mvnorm'.")
      }
      if(is.null(args$n)){
        stop("For `method` 'simulation', you must pass `n`.")
      }
      if(args$distribution == "dirichlet"){
        if(is.null(args$alpha)){
          if(is.null(args$mu)){
            stop("For `method` 'simulation' and `distribution` 'dirichlet', you must pass either `mu` or `alpha`.")
          }else{
            if(is.null(args$precision)){stop("For `method` 'simulation' and `distribution` 'dirichlet', if you do not pass `alpha` you must pass `precision`.")}
            if(args$precision <= 0){stop("For `method` 'simulation' and `distribution` 'dirichlet', `precision` must be non-negative.")}
            args$alpha <- args$mu*args$precision
          }
        }else{
          if(min(args$alpha) <= 0){
            stop("For `method` 'simulation' and `distribution` 'dirichlet', all elements of `alpha` must be positive.")
          }
        }
        args$sims <- gtools::rdirichlet(n, args$alpha)
      }else if(args$distribution == "mvnorm"){
        if(is.null(args$mu) | is.null(args$sigma)){
          stop("For `method` 'simulation' and `distribution` 'mvnorm', you must pass `mu` and `sigma`.")
        }
        args$sims <- mvtnorm::mvnorm(n, args$mu, args$sigma)
      }

      plurality_tie_probs_from_sims(args$sims, ...) # tol, cand_names, sep

    }
  }else if(method == "grid-based"){
    # I put the error reporting in the method in this case.
    grid_based_plurality_pivot_probs(distribution = args$distribution, increment = args$increment, ...)
  }

}




