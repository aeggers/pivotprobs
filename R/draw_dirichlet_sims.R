## get dirichlet sims

draw_dirichlet_sims <- function(N, v_vec, s){

  if(!sum(v_vec) == 1){
    v_vec <- v_vec/sum(v_vec)
  }

  gtools::rdirichlet(N, v_vec*s)
}


