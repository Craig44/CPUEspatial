#' simulate_data
#' 
#' @details 
#' A function for simulating data from TMB model (specifically CPUEspatail obj). The idea is for the simulated data to be used for measuring consistentency with observed data using an idea analgous to prediction distributions from the Bayesian framework
#' This will feed into the DHARMa R package see ?DHARMa for more information
#' @param obj TMB object which has been optimised and checkef that it successfully converged
#' @param n_sims number of simulated data sets 
#' @param fixed_effect integer[0,1] if you want to use plug-in (=0) method or Pseudo Bayes (=1) where you simulate values from the MVN joint distribution
#' @param random_effect integer[0,1] specifying how random effects are drawn. 0 = Bayes estiamtes, 1 = simulate a new field
#' @param calibrate boolean on whether we should try and calibrate the P-value.
#' @export
#' @return: A matrix of simulated data, with rows being each data record and cols being simulations i.e. n_sims
simulate_data <- function(obj, n_sims = 100, fixed_effect = 0, random_effect = 0, calibrate = FALSE) {
  ## get MLE and Bayes estimates
  MLE_pars = obj$env$last.par.best
  ## set simulate random effects or not
  obj$env$data$simulate_GF = random_effect
  sim_data = NULL
  if(fixed_effect == 1) {
    sd_rep = sdreport(obj, getJointPrecision = T)
    # cols are the pars
    sim_pars = rmvnorm_prec(MLE_pars, sd_rep$jointPrecision, n_sims)
    sim_data =sapply(1:n_sims, FUN = function(x){
      #print(x); 
      this_sim = obj$simulate(sim_pars[,x])
      ## you can look at other information here if you want
      this_sim$y_i
    });
  } else {
    sim_data = sapply(1:n_sims, FUN = function(x){
      #print(x); 
      this_sim = obj$simulate(MLE_pars)
      ## you can look at other information here if you want
      this_sim$y_i
    });
  }
  return(sim_data)
}
