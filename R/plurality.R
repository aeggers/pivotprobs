#' Compute pivot probabilities for a plurality election
#'
#' @param method One of
#' \itemize{
#' \item \code{"eggers-vivyan"}, which is based on the Dirichlet
#' distribution but makes an independence assumption for efficiency -- in brief, it assumes that the probability of $a$ and $b$ tying for first at vote share $x$ is equal to the probability of $a$ and $b$ receiving $x$ times the probability of $c$ receiving no more than $x$ times the probability of $d$ receiving no more than $x$ times the probability of $e$ receiving no more than $x$, etc. (In reality, these events are not independent.)
#' \item \code{"myatt-fisher"}, which computes exact (but not normalized) pivot probabilities for the three-candidate case assuming a Dirichlet distribution over election outcomes.
#' \item \code{"simplicial-cubature"}, which uses an adaptive algorithm to integrate the distribution of election results (either Dirichlet or logistic normal) over the relevant simplices.
#' \item \code{"simulation"}, which takes a matrix of simulated election results from any distribution, or arguments to generate such a matrix from a Dirichlet or logistic normal distribution, and computes pivot probabilities.
#' \item \code{"grid-based"}, which constructs a grid of points and computes pivot probabilities using either Dirichlet or logistic normal parameters.
#' }
#' @param ... Other arguments including
#' \itemize{
#' \item \code{distribution}, either "dirichlet" or "logisticnormal", relevant unless \code{method = "simulation"} and user passes \code{sims}; if \code{distribution} is not supplied, code tries to infer \code{distribution} from parameters passed (\code{alpha} implies "dirichlet", \code{sigma} implies "logisticnormal"),
#' \item \code{alpha}, a$k$-length vector of parameters for the Dirichlet
#'  distribution,
#'  \item \code{mu}, a $k$-length vector of expected vote shares for the Dirichlet distribution or a $k-1$ length vector of means for thre logistic normal distribution,
#'  \item \code{sigma}, a $k-1$-by-$k-1$ matrix of variances and covariances for the logistic normal distribution,
#'  \item \code{sims} to be used if \code{method = "simulation"} (not required if you pass parameters for the simulations to be generated from a Dirichlet or logistic normal distribution),
#'  \item \code{n} the number of simulations, if \code{method = "simulation"} and not passing \code{sims},
#'  \item \code{increment}, the space between grid-points if using \code{method = "grid-based"}
#'  }
#'
#'
#' @export
#'
plurality_pivot_probs <- function(method = "eggers-vivyan", ...){
  if(method == "eggers-vivyan"){
    eggers_vivyan_plurality_pivot_probs(...)
  }else if(method == "myatt-fisher"){
    myatt_fisher_plurality_pivot_probs(...)
  }else if(method == "simulation"){
    simulation_based_plurality_pivot_probs(...)
  }else if(method == "grid-based"){
    grid_based_plurality_pivot_probs(...)
  }else if(method == "simplicial-cubature"){
    simplicial_cubature_based_plurality_pivot_probs(...)
  }else{
    stop("Unknown method ", method, " for computing plurality pivot probs.")
  }
}

# this is a bad old method that was complicated to debug
plurality_pivot_probs2 <- function(method = "eggers-vivyan", normalize = T, ...){

  args <- list(...)

  # guessing the distribution from the parameters passed
  if(is.null(args$sims) & is.null(args$distribution)){
    if((!is.null(args$alpha) | !is.null(args$precision)) & !is.null(args$sigma)){
      stop("You did not pass a `distribution` argument and you passed both (`alpha` or `precision`) and `sigma` -- unable to determine which distribution you intended. Either state the distribution or pass unambiguous arguments.")
    }
    if(!is.null(args$alpha) | !is.null(args$precision)){
      args$distribution <- "dirichlet"
    }else if(!is.null(args$sigma)){
      args$distribution <- "logisticnormal"
    }else{
      stop("You did not pass a `distribution` argument and you did not pass arguments that make clear what distribution you intended.")
    }
  }

  normalizing_factor <- 1 # previously normalized for multivariate normal

  if(method == "eggers-vivyan"){

    # must pass either mu (length k) and precision (positive) or alpha (length k, positive)
    if(is.null(args$alpha)){
      if(is.null(args$mu) | is.null(args$precision)){
        stop("For method eggers-vivyan, you must pass either `mu` and `precision` or `alpha`.")
      }
      if(length(args$precision) > 1){stop("Precision argument must be a scalar.")}
      if(args$precision <= 0){stop("For method eggers-vivyan, `precision` must be non-negative.")}
      args$alpha <- args$mu*args$precision
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

  }else if(method == "simplicial-cubature"){

    # Pivot probability calculation via numerical integration of the relevant facets.
    # Must pass a function name and parameters

    if(!args$distribution %in% c("dirichlet", "logisticnormal")){
      stop("For method `simplicial-cubature`, the `distribution` argument must be `dirichlet` or `logisticnormal`.")
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

      # take out arguments that don't apply here
      args$distribution <- NULL
      args$mu <- NULL
      args$precision <- NULL
      sc_out <- exec("simplicial_cubature_plurality_pivot_probs_dirichlet", !!!args)
     #  sc_out <- simplicial_cubature_plurality_pivot_probs_dirichlet(alpha = args$alpha, ...)

    }else if(args$distribution == "logisticnormal"){
      if(is.null(args$mu) | is.null(args$sigma)){
        stop("For `method` 'simplicial-cubature' and `integrand_function` 'logisticnormal', you must pass `mu` and `sigma`.")
      }

      args$distribution = NULL
      args$normalizing_factor = normalizing_factor
      sc_out <- exec("simplicial_cubature_plurality_pivot_probs_logisticnormal", !!!args) # converts args into explicit arguments.
    }

    sc_out

  }else if(method == "simulation"){
    # you can pass the matrix of simulations, or
    if(is.null(args$sims)){
      if(is.null(args$distribution)){
        stop("For `method` 'simulation' you must pass either `sims` or a `distribution` argument.")
      }
      if(! args$distribution %in% c("dirichlet", "logisticnormal")){
        stop("For `method` 'simulation' the `distribution` argument must be 'dirichlet' or 'logisticnormal'.")
      }
      if(is.null(args$nsims)){
        stop("For `method` 'simulation', you must pass `nsims`. (args was ", paste(args, collapse = "\n"), ".\n")
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
        args$sims <- gtools::rdirichlet(args$nsims, args$alpha)
      }else if(args$distribution == "logisticnormal"){
        if(is.null(args$mu) | is.null(args$sigma)){
          stop("For `method` 'simulation' and `distribution` 'logisticnormal', you must pass `mu` and `sigma`.")
        }
        args$sims <- rlogisticnormal(args$nims, args$mu, args$sigma)
      }

    }

    plurality_tie_probs_from_sims(args$sims, tol = args$tol, cand_names = args$cand_names, sep = args$sep) # explicitly naming arguments here.

  }else if(method == "grid-based"){
    # I put the error reporting in the method in this case.
    args$normalizing_factor = normalizing_factor
    exec("grid_based_plurality_pivot_probs", !!!args) # (distribution = args$distribution, increment = args$increment, ...)
  }else{
    stop("Unknown method ", args$method, ".\n")
  }

}




