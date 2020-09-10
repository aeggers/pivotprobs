#' @export
positional_pivot_probs_froms_sims <- function(sims, tol = .01, n = 1000, s = .5, cand_names = NULL, sep = ""){

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  if(is.null(cand_names)){
    cand_names <- names(sims)
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
  }

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  normalization_factor <- n*tol # used to be tol/sqrt(3)
  # this is direct Monte Carlo -- nothing fancy
  out <- list(
    mean(score_a > score_b & score_a - score_b < tol & score_b > score_c)/normalization_factor,
    mean(score_a > score_c & score_a - score_c < tol & score_c > score_b)/normalization_factor,
    mean(score_b > score_c & score_b - score_c < tol & score_c > score_a)/normalization_factor
  )

  names(out) <- c(paste0(cand_names[1], sep, cand_names[2]), paste0(cand_names[1], sep, cand_names[3]), paste0(cand_names[2], sep, cand_names[3]))

  out

}

positional_intersection_point <- function(pvec3, s){
  # given a vector of second preference rations (pAB, pBA, pCA)
  # returns the intersection point, which may not be in the simplex.
  pAB <- pvec3[1]
  pAC <- 1 - pvec3[1]
  pBA <- pvec3[2]
  pBC <- 1 - pvec3[2]
  pCA <- pvec3[3]
  pCB <- 1 - pvec3[3]

  d1 = 1 - s*(pAB - pCB + pCA)
  s1 = (1 - s*(pBA + pCB - pCA))
  i1 = s*(pCB - pCA)

  d2 = 2 - s*(pAC + pCA)
  s2 = (s*(pBC + pCA - pBA) - 1)
  i2 = (1 - s*pCA)

  # then substitute one into the other to solve for v_B
  v_b.star = ((i1/d1)*d2 - i2)/(s2 - (d2*s1/d1))
  v_a.star = (v_b.star*s1 + i1)/d1
  c(v_a.star, v_b.star, 1 - v_b.star - v_a.star)

}

dirichlet_function <- function(x, alpha_vec){
  gtools::ddirichlet(as.numeric(x), alpha_vec)
}


piv_probs_from_p_vec <- function(p_vec, positional_s, fp_alpha){

  full_p_vec <- c(p_vec[1], 1 - p_vec[1], p_vec[2], 1 - p_vec[2], p_vec[3], 1 - p_vec[3])
  tpms <- votevizr::ternary_point_mats_from_p_vec_and_s_v3(full_p_vec, positional_s)
  pip <- positional_intersection_point(p_vec, positional_s)

  ## If the intersection point is not on the simplex, then one of the tie lines is length 0, and the others are as is
  if(min(pip) < 0){
    tpms[[which(pip < 0)]] = matrix(0, nrow = 2, ncol = 3)
  }else{
    # otherwise the second point in each one should be the intersection point
    for(j in 1:3){
      tpms[[j]][2,] <- pip
    }
  }

  # now we have endpoints of each of the lines we need to integrate along

  # use SimplicialCubature to do integration
  unlist(list(
    "AB" = SimplicialCubature::adaptIntegrateSimplex(dirichlet_function, S = t(tpms[[1]]), alpha_vec = fp_alpha)$integral,
    "AC" = SimplicialCubature::adaptIntegrateSimplex(dirichlet_function, S = t(tpms[[2]]), alpha_vec = fp_alpha)$integral,
    "BC" = SimplicialCubature::adaptIntegrateSimplex(dirichlet_function, S = t(tpms[[3]]), alpha_vec = fp_alpha)$integral
  ))

}

piv_probs_positional_2_step <- function(v_vec, s, positional_s, N = 100){

  # reformulate the Dirichlet parameters
  alpha_vec <- v_vec*s
  fp_alpha <- c(sum(alpha_vec[1:2]), sum(alpha_vec[3:4]), sum(alpha_vec[5:6]))

  # draw second preference shares from beta distribution
  p_vecs <- data.frame(
    p_ab = rbeta(N, alpha_vec[1], alpha_vec[2]),
    p_ba = rbeta(N, alpha_vec[3], alpha_vec[4]),
    p_ca = rbeta(N, alpha_vec[5], alpha_vec[6])
  )

  # compute pivotal probabilities for each set of second preference shares
  piv_probs <- apply(p_vecs, 1, piv_probs_from_p_vec, positional_s, fp_alpha)

  piv_probs

}






## importance sampling test

importance_sample_positional <- function(N, positional_s, v_vec, s, v_vec_is = rep(1, 6)/6, s_is = 6, tol = .01){

  # draw simulations from is distribution
  sims <- draw_dirichlet_sims(N, v_vec_is, s_is)

  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  piv_df <- data.frame(
    ab = as.integer(score_a > score_b & score_a - score_b < tol & score_b > score_c),
    ac = as.integer(score_a > score_c & score_a - score_c < tol & score_c > score_b),
    bc = as.integer(score_b > score_c & score_b - score_c < tol & score_c > score_a),
    dens_is = gtools::ddirichlet(sims, alpha = v_vec_is*s_is),
    dens_true = gtools::ddirichlet(sims, alpha = v_vec*s)
  )

  piv_df$dens_ratio <- piv_df$dens_true/piv_df$dens_is

  list(
    "AB" = (piv_df$ab %*% piv_df$dens_ratio)/(tol*sum(piv_df$dens_ratio)),
    "AC" = (piv_df$ac %*% piv_df$dens_ratio)/(tol*sum(piv_df$dens_ratio)),
    "BC" = (piv_df$bc %*% piv_df$dens_ratio)/(tol*sum(piv_df$dens_ratio))
  )

}
